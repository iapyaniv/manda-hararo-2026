from __future__ import annotations

import json
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import patch

import numpy as np

from mh26.distributed import (
    _strike_vectors,
    assemble,
    build_geometry,
    geometry_from_json,
    geometry_sha256,
    geometry_to_json,
    load_geometry_source,
)
from mh26.gunw import block_average, corrected_los, los_unit_neu, valid_sample_mask
from mh26.params import CENTER_E, CENTER_N, FARFIELD_RADIUS, WAVELENGTH_L
from mh26.s1_check import coherence_fraction
from mh26.stress import fault_patch_receivers
from mh26.uncertainty import referenced_noise_draw
from mh26.wright2006 import DAT_DEFAULT, build_contract, coverage_report


class ValidityMaskTests(unittest.TestCase):
    def test_decodes_documented_digits_and_provisional_high_flags(self) -> None:
        values = np.array([11, 12, 111, 124, 0, 1, 10, 100, 255,
                           0x0100000C, 0x01000000], dtype=np.uint32)
        expected = np.array([True, True, True, True, False, False, False,
                             False, False, True, False])
        np.testing.assert_array_equal(valid_sample_mask(values), expected)

    def test_corrected_los_applies_every_required_screen(self) -> None:
        shape = (1, 7)
        data = {
            "unwrappedPhase": np.ones(shape),
            "ionospherePhaseScreen": np.zeros(shape),
            "coherence": np.array([[0.9, 0.9, 0.2, 0.9, 0.9, 0.9, 0.9]]),
            "connectedComponents": np.array([[1, 1, 1, 0, 65535, 1, 1]], dtype=np.uint16),
            "mask": np.array([[11, 10, 11, 11, 11, 0x0100000C, 255]], dtype=np.uint32),
        }
        los, _ = corrected_los(data)
        finite = np.isfinite(los)
        np.testing.assert_array_equal(
            finite, [[True, False, False, False, False, True, False]]
        )
        self.assertAlmostEqual(los[0, 0], -WAVELENGTH_L / (4 * np.pi))


class NumericUtilityTests(unittest.TestCase):
    def test_empty_blocks_do_not_emit_runtime_warning(self) -> None:
        array = np.full((4, 4), np.nan)
        array[:2, :2] = 2.0
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            points, values, (_, _, grid) = block_average(
                array, np.arange(4.0), np.arange(4.0), block=2
            )
        self.assertEqual(caught, [])
        self.assertEqual(len(points), 1)
        np.testing.assert_allclose(values, [2.0])
        self.assertTrue(np.isnan(grid[1, 1]))

    def test_uint8_coherence_is_scaled_to_fraction(self) -> None:
        scaled = coherence_fraction(np.array([0, 64, 255], dtype=np.uint8))
        np.testing.assert_allclose(scaled, [0, 64 / 255, 1])
        unchanged = coherence_fraction(np.array([0.0, 0.25, 1.0]))
        np.testing.assert_allclose(unchanged, [0.0, 0.25, 1.0])

    def test_los_interpolation_returns_unit_neu_vectors(self) -> None:
        data = {
            "cube_heightAboveEllipsoid": np.array([0.0, 1000.0]),
            "cube_yCoordinates": np.array([0.0, 1000.0]),
            "cube_xCoordinates": np.array([0.0, 1000.0]),
            "cube_losUnitVectorX": np.full((2, 2, 2), 0.6),
            "cube_losUnitVectorY": np.full((2, 2, 2), 0.2),
        }
        vectors = los_unit_neu(data, np.array([[500.0, 500.0]]), height=500.0)
        np.testing.assert_allclose(vectors[0, :2], [0.2, 0.6])
        self.assertAlmostEqual(np.linalg.norm(vectors[0]), 1.0)


class GeometryConventionTests(unittest.TestCase):
    MODEL = np.array([683225.36, 1350366.91, 152.58, 17150.0, 1494.0,
                      4989.0, 2.812])

    def test_strike_frame_round_trip_matches_map_transform(self) -> None:
        # Regression: an opposite perpendicular sign in the map transform
        # swapped the east and west ridges in the published figure.
        from mh26.faults import strike_frame
        e_c, n_c, strike = self.MODEL[0], self.MODEL[1], self.MODEL[2]
        sr = np.radians(strike)
        s, p = 3_000.0, 1_650.0
        e = e_c + s * np.sin(sr) - p * np.cos(sr)
        n = n_c + s * np.cos(sr) + p * np.sin(sr)
        along, perp = strike_frame(np.array([e]), np.array([n]), self.MODEL)
        self.assertAlmostEqual(along[0], s, places=5)
        self.assertAlmostEqual(perp[0], p, places=5)

    def test_fault_columns_terminate_at_the_dike_plane(self) -> None:
        # Regression: ridge-anchored columns once continued through the dike
        # onto the far side of the rift.
        offsets = {"east": np.full(14, 500.0), "west": np.full(14, -700.0)}
        g = build_geometry(self.MODEL, fault_trace_offsets=offsets)
        _, perp = _strike_vectors(self.MODEL[2])
        dip_r = np.radians(g["fault_dip"])
        for side, sign in (("east", 1.0), ("west", -1.0)):
            sel = g["kinds"] == f"fault_{side}"
            pat = g["patches"][sel]
            off = ((pat[:, 0] - self.MODEL[1]) * perp[0]
                   + (pat[:, 1] - self.MODEL[0]) * perp[1])
            self.assertTrue((sign * off > 0).all())
            deepest = sign * off - pat[:, 8] * np.cos(dip_r)
            self.assertGreater(float(deepest.min()), -1.0)
            self.assertEqual(int(sel.sum()),
                             int(g["fault_masks"][side].sum()))
        full = build_geometry(self.MODEL)
        for side in ("east", "west"):
            self.assertTrue(full["fault_masks"][side].all())

    def test_fault_patch_centres_lie_on_the_dipping_planes(self) -> None:
        # Patch centres must lie on the dipping plane — z/tan(dip) toward
        # the axis from the surface trace — not vertically below the trace.
        offsets = {"east": np.linspace(1_200.0, 3_000.0, 14),
                   "west": np.full(14, -1_400.0)}
        g = build_geometry(self.MODEL, fault_trace_offsets=offsets)
        _, perp = _strike_vectors(self.MODEL[2])
        dip_r = np.radians(g["fault_dip"])
        for side, sign in (("east", 1.0), ("west", -1.0)):
            cols = np.argwhere(g["fault_masks"][side])[:, 1]
            sel = g["kinds"] == f"fault_{side}"
            pat = g["patches"][sel]
            off = ((pat[:, 0] - self.MODEL[1]) * perp[0]
                   + (pat[:, 1] - self.MODEL[0]) * perp[1])
            expected = (np.asarray(offsets[side])[cols]
                        - sign * pat[:, 2] / np.tan(dip_r))
            np.testing.assert_allclose(off, expected, atol=1e-6)

    def test_saved_geometry_round_trip_preserves_every_patch(self) -> None:
        geometry = build_geometry(self.MODEL)
        restored = geometry_from_json(geometry_to_json(geometry))
        self.assertEqual(geometry_sha256(restored), geometry_sha256(geometry))
        np.testing.assert_array_equal(restored["patches"], geometry["patches"])
        for side in ("east", "west"):
            self.assertEqual(len(restored["fault_column_offsets_m"][side]), 14)

    def test_complete_geometry_file_preserves_source_patches_and_identity(self) -> None:
        geometry = build_geometry(self.MODEL, fault_offset=1_500.0,
                                  fault_dip=55.0)
        wrapper = {
            "geometry": geometry_to_json(geometry),
            "geometry_sha256": geometry_sha256(geometry),
        }
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "old_geometry.json"
            source.write_text(json.dumps(wrapper))
            restored, provenance = load_geometry_source(source)
            wrapper["geometry_sha256"] = "0" * 64
            source.write_text(json.dumps(wrapper))
            with self.assertRaisesRegex(ValueError, "hash does not match"):
                load_geometry_source(source)
        np.testing.assert_array_equal(restored["patches"], geometry["patches"])
        self.assertEqual(provenance["geometry_sha256"],
                         geometry_sha256(geometry))
        self.assertEqual(provenance["source_kind"], "distributed_result")
        self.assertTrue(provenance["sha256"])

    def test_signed_cross_rift_increment_follows_the_ene_frame(self) -> None:
        # Regression: an opposite perpendicular sign in a map transform once
        # swapped the east and west ridges. The signed phase increment must
        # follow the _strike_vectors ENE-positive frame, and its aliasing
        # mark must trigger exactly beyond λ/4 of LOS per pixel.
        from mh26.summary_figures import (
            aliasing_flag,
            signed_cross_rift_increment,
        )
        strike = self.MODEL[2]
        _, perp = _strike_vectors(strike)
        x = 683_000.0 + 30.0 * np.arange(40)
        y = 1_350_000.0 - 30.0 * np.arange(40)
        X, Y = np.meshgrid(x, y)
        rate = 0.05  # rad per metre, gradient pointing ENE
        phase = rate * (perp[1] * X + perp[0] * Y)
        signed = signed_cross_rift_increment(
            phase, np.ones_like(phase), x, y, strike
        )
        np.testing.assert_allclose(signed[1:-1, 1:-1], rate * 30.0, atol=1e-9)

        los = 0.004 * (perp[1] * X + perp[0] * Y)  # 12 cm per pixel > λ/4
        flagged = aliasing_flag(los, y, strike, WAVELENGTH_L / 4.0)
        self.assertTrue(flagged[1:-1, 1:-1].all())
        calm = aliasing_flag(0.1 * los, y, strike, WAVELENGTH_L / 4.0)
        self.assertFalse(calm[1:-1, 1:-1].any())

    def test_dike_frame_rotation_matches_a_direct_rotated_solve(self) -> None:
        # The rotated formal sigmas need the east–north cross-covariance
        # (large for near-polar looks); dropping it fails silently, so the
        # rotation is checked against a solve done in the rotated basis.
        from mh26.decompose3d import dike_frame, solve_cells_3d
        rng = np.random.default_rng(7)
        shape = (2, 3)
        cells = []
        for k, tilt in enumerate((0.15, -0.15, 0.25, -0.05)):
            unit = np.array([tilt, 0.6 + 0.05 * k, 0.75])
            unit /= np.linalg.norm(unit)
            cells.append({"los": rng.normal(size=shape),
                          "n": np.full(shape, unit[0]),
                          "e": np.full(shape, unit[1]),
                          "u": np.full(shape, unit[2]),
                          "sigma": 0.5 + 0.2 * k})
        fields = solve_cells_3d(cells)
        strike = self.MODEL[2]
        frame = dike_frame(fields, strike)
        along, perp = _strike_vectors(strike)
        weights = np.array([1.0 / c["sigma"] for c in cells])
        for r in range(shape[0]):
            for c in range(shape[1]):
                A = np.column_stack([
                    [cl["e"][r, c] * along[1] + cl["n"][r, c] * along[0]
                     for cl in cells],
                    [cl["e"][r, c] * perp[1] + cl["n"][r, c] * perp[0]
                     for cl in cells],
                    [cl["u"][r, c] for cl in cells],
                ])
                d = np.array([cl["los"][r, c] for cl in cells])
                cov = np.linalg.inv((A * weights[:, None]).T @ (A * weights[:, None]))
                m = cov @ ((A * weights[:, None]).T @ (d * weights))
                sigma = np.sqrt(np.diag(cov))
                for i, key in enumerate(("par", "perp", "u")):
                    self.assertAlmostEqual(frame[key][r, c], m[i], places=9)
                    self.assertAlmostEqual(frame[f"sigma_{key}"][r, c],
                                           sigma[i], places=9)
        # Strike 90° makes along = east and perpendicular-left = north.
        identity = dike_frame(fields, 90.0)
        np.testing.assert_allclose(identity["par"], fields["e"], atol=1e-12)
        np.testing.assert_allclose(identity["perp"], fields["n"], atol=1e-12)
        np.testing.assert_allclose(identity["sigma_perp"], fields["sigma_n"],
                                   atol=1e-12)

    def test_idealized_geometry_has_stress_receivers(self) -> None:
        geometry = build_geometry(self.MODEL)
        receivers = fault_patch_receivers(geometry, 65.0)
        for side in ("east", "west"):
            self.assertGreater(len(receivers[side]["points"]), 0)
            saved = geometry["patches"][geometry["kinds"] == f"fault_{side}"]
            np.testing.assert_allclose(receivers[side]["points"], saved[:, :3],
                                       atol=1e-9)


class Wright2006ContractTests(unittest.TestCase):
    def test_dat_segments_chain_and_contract_round_trips(self) -> None:
        geometry = build_contract(DAT_DEFAULT, extend_se=2, with_faults=False)
        prov = geometry["wright2006_provenance"]
        # Sub-metre junctions are what validate the inferred column meaning;
        # growth here means the file format or the parser changed.
        self.assertLessEqual(max(prov["junction_residuals_m"]), 0.5)
        restored = geometry_from_json(geometry_to_json(geometry))
        self.assertEqual(geometry_sha256(restored), geometry_sha256(geometry))
        # The scalar frame must be the chord of the s axis: origin at s=0,
        # azimuth toward increasing s, within the recorded chord error.
        along, _ = _strike_vectors(geometry["strike"])
        e_c, n_c = geometry["trace_centre"]
        ns = geometry["dike_shape"][1]
        top_row = geometry["patches"][:ns]
        reconstructed = (np.array([n_c, e_c])
                         + np.outer(geometry["dike_s_centres"], along))
        errors = np.linalg.norm(top_row[:, :2] - reconstructed, axis=1)
        self.assertLessEqual(errors.max(),
                             prov["straight_frame_max_error_m"] + 0.1)
        self.assertLess(errors.max(), 6_000.0)

    def test_coverage_projection_detects_an_uncovered_trace(self) -> None:
        geometry = build_contract(DAT_DEFAULT, extend_se=0, with_faults=False)
        coverage = coverage_report(geometry)
        # The 2026 uniform trace overruns the unextended 2005 plane SE end
        # by ~2.8 km; a projection that clamps to the plane cannot see that.
        self.assertLess(coverage["margin_se_km"], -2.0)
        self.assertGreater(coverage["margin_nw_km"], 20.0)


class ReferenceAndCacheTests(unittest.TestCase):
    MODEL = GeometryConventionTests.MODEL

    def test_noise_draw_obeys_far_field_plane_constraints(self) -> None:
        angles = np.linspace(0, 2 * np.pi, 24, endpoint=False)
        radius = FARFIELD_RADIUS + 5_000.0
        far_points = np.column_stack([
            CENTER_N + radius * np.sin(angles),
            CENTER_E + radius * np.cos(angles),
        ])
        near_points = np.array([[CENTER_N, CENTER_E],
                                [CENTER_N + 1_000.0, CENTER_E]])
        points = np.vstack([far_points, near_points])
        draw = referenced_noise_draw(
            points, np.eye(len(points)), np.random.default_rng(4)
        )
        design = np.column_stack([
            np.ones(len(far_points)),
            far_points[:, 1] - CENTER_E,
            far_points[:, 0] - CENTER_N,
        ])
        np.testing.assert_allclose(design.T @ draw[:len(far_points)], 0.0,
                                   atol=1e-8)

    def test_greens_cache_is_reused_only_for_the_same_geometry(self) -> None:
        first = build_geometry(self.MODEL)
        second = build_geometry(self.MODEL, fault_offset=1_500.0)
        dataset = {
            "pts": np.zeros((2, 2)),
            "units": np.zeros((2, 3)),
            "vals": np.zeros(2),
            "sigma": 1.0,
        }
        calls = []

        def fake_greens(geometry, points, units):
            calls.append(geometry["fault_offset"])
            return np.full((len(points), len(geometry["kinds"])), len(calls))

        with patch("mh26.distributed.greens_los", side_effect=fake_greens):
            assemble([dataset], first)
            assemble([dataset], first)
            matrix, _, _ = assemble([dataset], second)
        self.assertEqual(calls, [2_000.0, 1_500.0])
        np.testing.assert_array_equal(matrix[:, :-1], 2.0)


if __name__ == "__main__":
    unittest.main()
