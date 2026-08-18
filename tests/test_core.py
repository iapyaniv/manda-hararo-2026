from __future__ import annotations

import unittest
import warnings

import numpy as np

from mh26.gunw import block_average, corrected_los, los_unit_neu, valid_sample_mask
from mh26.params import WAVELENGTH_L
from mh26.s1_check import coherence_fraction


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


if __name__ == "__main__":
    unittest.main()
