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
                # U∥ is positive toward strike + 180° (near-north), so the
                # rotated-basis solve uses the negated along vector.
                A = np.column_stack([
                    [-(cl["e"][r, c] * along[1] + cl["n"][r, c] * along[0])
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
        # Strike 90° makes along = east (so U∥ = −east) and
        # perpendicular-left = north.
        identity = dike_frame(fields, 90.0)
        np.testing.assert_allclose(identity["par"], -fields["e"], atol=1e-12)
        np.testing.assert_allclose(identity["perp"], fields["n"], atol=1e-12)
        np.testing.assert_allclose(identity["sigma_perp"], fields["sigma_n"],
                                   atol=1e-12)

        # The figure once drew +U∥ opposite to the solved component. Check
        # the rendered shaft in map coordinates, including its label.
        import matplotlib.pyplot as plt
        from matplotlib.text import Annotation
        from mh26.summary_figures import (
            X_MIN, X_MAX, Y_MIN, Y_MAX, draw_frame_axes,
        )
        fig, ax = plt.subplots(figsize=(4, 5))
        try:
            ax.set(xlim=(X_MIN, X_MAX), ylim=(Y_MIN, Y_MAX))
            draw_frame_axes(ax, strike)
            fig.canvas.draw()
            arrows = [text for text in ax.texts
                      if isinstance(text, Annotation)]
            self.assertEqual(len(arrows), 2)
            for label, component in (("U⊥", "perp"), ("U∥", "par")):
                position = next(text.get_position() for text in ax.texts
                                if text.get_text() == label)
                arrow = min(arrows, key=lambda a:
                            np.linalg.norm(np.asarray(a.xy) - position))
                shaft = arrow.arrow_patch.get_path().vertices[:3]
                display = arrow.arrow_patch.get_transform().transform(shaft)
                mapped = ax.transData.inverted().transform(display)
                direction = mapped[-1] - mapped[0]
                direction /= np.linalg.norm(direction)
                projected = dike_frame(
                    {**fields, "e": direction[0], "n": direction[1]}, strike)
                self.assertAlmostEqual(projected[component], 1.0, places=9)
        finally:
            plt.close(fig)

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


class SummaryPatchFrameTests(unittest.TestCase):
    MODEL = GeometryConventionTests.MODEL

    def test_patch_frame_matches_the_straight_frame_on_a_straight_plane(self) -> None:
        from mh26.faults import strike_frame
        from mh26.summary_figures import (
            patch_fault_offset_km,
            patch_profile_frame,
            patch_top_trace,
            straight_profile_frame,
        )
        geometry = build_geometry(self.MODEL, fault_offset=1_500.0)
        magnitudes = np.where(geometry["kinds"] == "dike", 1.0, 0.0)
        frame = patch_profile_frame(geometry, magnitudes)
        reference = straight_profile_frame(geometry)
        self.assertAlmostEqual(frame["strike_deg"] % 360.0,
                               reference["strike_deg"] % 360.0, places=6)
        _, perp_origin = strike_frame(np.array([frame["origin_e"]]),
                                      np.array([frame["origin_n"]]),
                                      self.MODEL)
        self.assertLess(abs(float(perp_origin[0])), 1e-6)
        east, north = patch_top_trace(geometry, "dike")
        finite = np.isfinite(east)
        along, perp = strike_frame(east[finite], north[finite], self.MODEL)
        np.testing.assert_allclose(perp, 0.0, atol=1e-6)
        self.assertAlmostEqual(float(along.min()),
                               float(geometry["dike_s_edges_m"][0]), places=4)
        self.assertAlmostEqual(float(along.max()),
                               float(geometry["dike_s_edges_m"][-1]), places=4)
        self.assertAlmostEqual(patch_fault_offset_km(geometry, frame), 1.5,
                               places=5)

    def test_patch_frame_follows_the_local_segment_on_the_wright_plane(self) -> None:
        # The chord frame is ~4° off the southernmost segment and its origin
        # sits at the plane's NW end; a profile built from it would cross
        # tens of kilometres from the 2026 opening. The patch frame must
        # follow the weighted local segment, flipped into the chord sense.
        from mh26.summary_figures import patch_profile_frame
        geometry = build_contract(DAT_DEFAULT, extend_se=2, with_faults=False)
        kinds = geometry["kinds"]
        magnitudes = np.zeros(kinds.shape)
        ns = geometry["dike_shape"][1]
        dike = np.flatnonzero(kinds == "dike")
        magnitudes[dike[ns - 3:ns]] = 1.0
        frame = patch_profile_frame(geometry, magnitudes)
        segment_strike = float(geometry["patches"][dike[ns - 1], 3])
        self.assertAlmostEqual(frame["strike_deg"],
                               (segment_strike + 180.0) % 360.0, places=6)
        self.assertNotAlmostEqual(frame["strike_deg"],
                                  float(geometry["strike"]) % 360.0, places=1)


class ReferenceAndCacheTests(unittest.TestCase):
    MODEL = GeometryConventionTests.MODEL

    def test_closure_replay_preserves_unknown_settings_and_checks_file_identity(self) -> None:
        # Legacy caches once acquired the current defaults as false
        # provenance. Same-shape input changes also escaped all validation.
        import contextlib
        import io
        import shutil
        from mh26 import closure, validate
        from scripts import build_release

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            results = root / "results"
            data = root / "data" / "gslc"
            source = root / "mh26"
            for path in (results, data, source):
                path.mkdir(parents=True)
            for name in ("closure.py", "gslc.py", "params.py"):
                shutil.copyfile(Path(closure.__file__).parent / name, source / name)
            (results / "okada_fit.json").write_text(
                json.dumps({"values": self.MODEL.tolist()}))
            grid = np.arange(-60.0, 90.0, 30.0)
            phase = np.zeros((5, 5), dtype=np.float32)
            first_tag = next(iter(closure.PRE_GRANULES))
            recorded_config = {"stat_cost_mode": "DEFO", "init_method": "MST"}
            with (patch.multiple(closure, REPO=root, RESULTS=results, DATA_GSLC=data,
                                 CENTER_E=self.MODEL[0], CENTER_N=self.MODEL[1],
                                 FARFIELD_RADIUS=50.0, NEAR_AXIS_PERP_M=20.0,
                                 NEAR_AXIS_ALONG_M=50.0, OFF_AXIS_RADIUS_M=60.0),
                  patch.multiple(validate, REPO=root, RESULTS=results),
                  patch.multiple(build_release, ROOT=root, DATA=root / "data"),
                  patch.object(closure, "_run_snaphu",
                               side_effect=AssertionError("snaphu must not run")),
                  contextlib.redirect_stdout(io.StringIO())):
                for tag in closure.PRE_GRANULES:
                    pre, (ref, sec) = closure.PRE_GRANULES[tag], closure.PAIRS[tag]
                    granules = {"pre": (pre, ref), "co": (ref, sec), "skip": (pre, sec)}
                    for leg in closure.LEGS:
                        np.savez_compressed(
                            closure._leg_file(tag, leg), phase=phase,
                            coherence=np.ones_like(phase),
                            x=self.MODEL[0] + grid, y=self.MODEL[1] + grid[::-1],
                            epsg=closure.EPSG_UTM37N, looks=[6, 3], granules=granules[leg])
                        raw = phase.copy()
                        if leg == "co":
                            raw[2, 2] = 2 * np.pi
                        metadata = ({"config": np.array(json.dumps(recorded_config))}
                                    if tag == first_tag and leg == "pre" else {})
                        np.savez_compressed(data / f"{tag}_loop{leg}_unwraw.npz",
                                            unw=raw, **metadata)
                cache = data / f"{first_tag}_loopco_unwraw.npz"
                original_cache = cache.read_bytes()
                with np.load(cache) as stored:
                    wrong_phase = stored["unw"] + 0.1
                np.savez_compressed(cache, unw=wrong_phase)
                with self.assertRaisesRegex(ValueError, "not congruent"):
                    closure.check(cached_only=True, regenerate=True)
                cache.write_bytes(original_cache)

                ifg = closure._leg_file(first_tag, "co")
                original_ifg = ifg.read_bytes()
                with np.load(ifg) as stored:
                    fields = {name: stored[name] for name in stored.files}
                invalid_fields = {
                    "granules": fields["granules"][::-1],
                    "epsg": np.array(4326), "looks": np.array([3, 6]),
                    "phase": phase[:-1], "coherence": np.ones(5),
                    "x": fields["x"] * 2, "y": fields["y"] + np.nan,
                }
                for name, invalid in invalid_fields.items():
                    with self.subTest(invalid_ifg_field=name):
                        np.savez_compressed(ifg, **{**fields, name: invalid})
                        with self.assertRaisesRegex(ValueError, f"interferogram {name}"):
                            closure.check(cached_only=True, regenerate=True)
                ifg.write_bytes(original_ifg)
                payload = closure.check(cached_only=True, regenerate=True)
                self.assertIsNone(payload[first_tag]["unwrap_configuration"]["co"])
                self.assertEqual(payload[first_tag]["unwrap_configuration"]["pre"],
                                 recorded_config)
                self.assertIsNone(payload[first_tag]["snaphu_execution"]["co"])
                self.assertGreater(payload[first_tag]["near_axis"]
                                   ["fraction_error_ge_1_cycle"], 0)
                verified = validate.validate_closure(require_inputs=True)
                self.assertEqual(set(verified["reproduced_tracks"]),
                                 set(closure.PRE_GRANULES))
                rounded = json.loads(json.dumps(payload))
                rounded[first_tag]["whole_scene"]["rms_wrapped_nonclosure_cycles"] += 1e-13
                with patch.object(closure, "_compute", return_value=rounded):
                    validate.validate_closure(require_inputs=True)
                    rounded[first_tag]["whole_scene"]["n_pixels"] += 1
                    with self.assertRaisesRegex(ValueError, "saved values differ"):
                        validate.validate_closure(require_inputs=True)
                with patch.object(closure, "figure", side_effect=AssertionError(
                        "a partial check must not replace the combined figure")):
                    closure.check(tags=[first_tag], cached_only=True)
                    with self.assertRaisesRegex(ValueError, "requires all"):
                        closure.check(tags=[first_tag], cached_only=True, regenerate=True)

                def snapshot():
                    return {path.relative_to(root): (path.read_bytes(), path.stat().st_mtime_ns)
                            for path in root.rglob("*") if path.is_file()}

                def replay_cli():
                    with patch("sys.argv", ["mh26.closure", "--check", "--cached-only"]):
                        closure.main()

                baseline = snapshot()
                closure.check()
                replay_cli()
                self.assertEqual(snapshot(), baseline)
                with patch.object(closure, "OFF_AXIS_RADIUS_M", 59.0):
                    with self.assertRaisesRegex(ValueError, "settings changed"):
                        replay_cli()
                    with self.assertRaisesRegex(ValueError, "settings changed"):
                        validate.validate_closure(require_inputs=True)
                    with patch.object(build_release, "write_archive") as archive_write:
                        with self.assertRaisesRegex(ValueError, "settings changed"):
                            build_release.build_closure_archive()
                        archive_write.assert_not_called()
                self.assertEqual(snapshot(), baseline)
                source_path = source / "closure.py"
                original_source = source_path.read_bytes()
                source_path.write_bytes(original_source + b"\n# changed source\n")
                changed_baseline = snapshot()
                with self.assertRaisesRegex(ValueError, "content hash mismatch"):
                    replay_cli()
                with self.assertRaisesRegex(ValueError, "content hash mismatch"):
                    validate.validate_closure(require_inputs=True)
                with patch.object(build_release, "write_archive") as archive_write:
                    with self.assertRaisesRegex(ValueError, "content hash mismatch"):
                        build_release.build_closure_archive()
                    archive_write.assert_not_called()
                self.assertEqual(snapshot(), changed_baseline)
                source_path.write_bytes(original_source)

                inventory_path = results / "loop_closure_inventory.json"
                original_inventory = inventory_path.read_bytes()
                for defect in ("schema", "duplicate", "unexpected", "unsafe", "identity"):
                    malformed = json.loads(original_inventory)
                    if defect == "schema":
                        malformed["schema_version"] = None
                    elif defect == "duplicate":
                        malformed["source_files"].append(malformed["source_files"][0])
                    elif defect == "unexpected":
                        malformed["reports"][0]["path"] = "results/unrelated.json"
                    elif defect == "unsafe":
                        malformed["input_files"][0]["path"] = "../outside.npz"
                    else:
                        malformed["input_files"][0]["size_bytes"] = "100"
                    inventory_path.write_text(json.dumps(malformed))
                    with self.subTest(invalid_inventory=defect), \
                            patch.object(build_release, "write_archive") as archive_write:
                        with self.assertRaises(ValueError):
                            build_release.build_closure_archive()
                        archive_write.assert_not_called()
                inventory_path.write_bytes(original_inventory)

                result_path = results / f"loop_closure_{first_tag}.json"
                original = result_path.read_bytes()
                changed = json.loads(original)
                changed["whole_scene"]["fraction_error_ge_1_cycle"] += 0.1
                result_path.write_text(json.dumps(changed))
                with self.assertRaisesRegex(ValueError, "content hash mismatch"):
                    validate.validate_closure()
                result_path.write_bytes(original)

                # These edits preserve shape and phase congruence, so the
                # stored content hashes must detect them.
                for path, key in ((closure._leg_file(first_tag, "pre"), "coherence"),
                                  (data / f"{first_tag}_loopco_unwraw.npz", "unw")):
                    original = path.read_bytes()
                    with np.load(path) as stored:
                        fields = {name: stored[name] for name in stored.files}
                    fields[key] = fields[key] + (-0.1 if key == "coherence" else 2 * np.pi)
                    np.savez_compressed(path, **fields)
                    with self.assertRaisesRegex(ValueError, "content hash mismatch"):
                        validate.validate_closure()
                    with self.assertRaisesRegex(ValueError, "content hash mismatch"):
                        closure.check()
                    with self.assertRaisesRegex(ValueError, "recorded closure input"):
                        closure.check(cached_only=True, regenerate=True)
                    path.write_bytes(original)

                (data / f"{first_tag}_loopco_unwraw.npz").unlink()
                available = validate.validate_closure()
                self.assertNotIn(first_tag, available["reproduced_tracks"])
                self.assertEqual(len(available["missing_inputs"]), 1)
                with self.assertRaises(FileNotFoundError):
                    validate.validate_closure(require_inputs=True)
                with self.assertRaises(FileNotFoundError):
                    closure.check()

    def test_three_epoch_report_cannot_choose_its_delivery_path(self) -> None:
        # A report whose source.directory pointed at a missing directory
        # once skipped the delivery hashes and the recomputation entirely.
        from mh26 import validate
        from mh26.params import REPO

        report = json.loads(
            (REPO / "results" / "uperp_three_epoch.json").read_text())
        report["source"]["directory"] = "data/nonexistent"
        report["epochs"]["sept05"]["volume_km3"] = 999.0
        with self.assertRaisesRegex(ValueError, "tracked\\s+delivery"):
            validate.check_three_epoch_report("uperp_three_epoch.json",
                                              "gunw", report)

    def test_snaphu_records_the_executed_binary_even_with_nonzero_version_status(self) -> None:
        # A Python distribution version once described a different system
        # executable. Capture the command identity at the execution boundary.
        import hashlib
        import subprocess
        from mh26 import gslc

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "snaphu-real"
            executable.write_bytes(b"test executable identity")
            alias = root / "snaphu"
            alias.symlink_to(executable)
            phase = np.arange(6, dtype=np.float32).reshape(2, 3)

            def run(command, **kwargs):
                self.assertEqual(command[0], str(executable))
                if command[1:] == ["--version"]:
                    return subprocess.CompletedProcess(command, 1, "", "snaphu v2.0.6\n")
                config = Path(command[2]).read_text().splitlines()
                output = next(line.split(maxsplit=1)[1] for line in config
                              if line.startswith("OUTFILE "))
                Path(output).write_bytes(phase.tobytes())
                return subprocess.CompletedProcess(command, 0, "", "")

            identity = {}
            with (patch.object(gslc, "DATA_GSLC", root),
                  patch.object(gslc.shutil, "which", return_value=str(alias)),
                  patch.object(gslc.subprocess, "run", side_effect=run) as executed):
                actual = gslc._run_snaphu(phase, np.ones_like(phase),
                                          execution_identity=identity)
            np.testing.assert_array_equal(actual, phase)
            self.assertEqual(executed.call_count, 2)
            self.assertEqual(identity, {
                "resolved_path": str(executable),
                "sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
                "reported_version": "snaphu v2.0.6", "version_exit_code": 1,
            })

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
        pts = np.array([
            [CENTER_N + 40e3, CENTER_E + 40e3],
            [CENTER_N + 40e3, CENTER_E - 40e3],
            [CENTER_N - 40e3, CENTER_E + 40e3],
            [CENTER_N - 40e3, CENTER_E - 40e3],
        ])
        # Orthogonal to the reference-plane basis over these four points, so
        # the values pass through assemble's projection unchanged.
        pattern = np.array([1.0, -1.0, -1.0, 1.0])
        dataset = {
            "pts": pts,
            "units": np.zeros((4, 3)),
            "vals": np.zeros(4),
            "sigma": 1.0,
        }
        calls = []

        def fake_greens(geometry, points, units):
            calls.append(geometry["fault_offset"])
            return len(calls) * pattern[:, None] * np.ones(
                (len(points), len(geometry["kinds"])))

        with patch("mh26.distributed.greens_los", side_effect=fake_greens):
            assemble([dataset], first)
            assemble([dataset], first)
            matrix, _, _ = assemble([dataset], second)
        self.assertEqual(calls, [2_000.0, 1_500.0])
        np.testing.assert_allclose(
            matrix[:, :-1],
            2.0 * pattern[:, None] * np.ones((4, matrix.shape[1] - 1)),
            atol=1e-9)


class ReferenceOperatorTests(unittest.TestCase):
    def test_model_projection_matches_observation_operator(self):
        from mh26.gunw import remove_farfield_plane
        from mh26.params import CENTER_E, CENTER_N, FARFIELD_RADIUS

        rng = np.random.default_rng(7)
        pts = np.column_stack([
            CENTER_N + rng.uniform(-45e3, 45e3, 4000),
            CENTER_E + rng.uniform(-45e3, 45e3, 4000),
        ])
        plane = (0.012 + 2.0e-7 * (pts[:, 1] - CENTER_E)
                 - 1.0e-7 * (pts[:, 0] - CENTER_N))
        out, coef = remove_farfield_plane(pts, plane)
        self.assertLess(float(np.abs(out).max()), 1e-9)
        np.testing.assert_allclose(coef, [0.012, 2.0e-7, -1.0e-7], rtol=1e-6)

        columns = rng.normal(size=(len(pts), 3))
        projected, _ = remove_farfield_plane(pts, columns)
        far = np.hypot(pts[:, 1] - CENTER_E,
                       pts[:, 0] - CENTER_N) > FARFIELD_RADIUS
        basis = np.column_stack([np.ones(len(pts)), pts[:, 1] - CENTER_E,
                                 pts[:, 0] - CENTER_N])
        refit, *_ = np.linalg.lstsq(basis[far], projected[far], rcond=None)
        self.assertLess(float(np.abs(refit[0]).max()), 1e-9)
        self.assertLess(float(np.abs(refit[1:]).max()) * 4.5e4, 1e-9)


if __name__ == "__main__":
    unittest.main()
