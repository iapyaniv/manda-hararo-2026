"""Shared constants for the Manda Hararo July 2026 dike analysis.

Constants shared by multiple processing steps are defined here.
"""

from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATA_NISAR = REPO / "data" / "nisar"
DATA_S1 = REPO / "data" / "s1"
RESULTS = REPO / "results"

# Radar wavelengths, metres.
WAVELENGTH_L = 0.2384     # NISAR L-band, 1257.5 MHz
WAVELENGTH_C = 0.055465   # Sentinel-1 C-band (LiCSAR)

# All NISAR grids in this dataset are UTM zone 37N.
EPSG_UTM37N = 32637

# Geometry-cube interpolation height, metres above ellipsoid; the rift floor is
# approximately 400--600 m above the ellipsoid.
LOS_HEIGHT_AE = 500.0

# Anomaly centre, UTM 37N metres (= 12.18 N, 40.66 E), used for
# far-field referencing and near-field statistics.
CENTER_E = 680_600.0
CENTER_N = 1_347_000.0

# Pixel screening: keep pixels with coherence above this AND a nonzero
# unwrapper connected component.
COHERENCE_MIN = 0.35

# Far-field plane (per-pair offset + ramp) is fitted beyond this radius
# from the centre, metres.
FARFIELD_RADIUS = 28_000.0

# Inversion block averaging: 16 px at 80 m posting = 1.28 km blocks.
BLOCK_PIXELS = 16
GRID_POSTING = 80.0

# Elastic half-space for the Okada dislocation model.
LAME_LAMBDA = SHEAR_MODULUS = 30e9  # Pa

# LiCSAR frame used for the independent Sentinel-1 cross-check.
S1_FRAME = "014A_07688_131313"
S1_PAIR = "20260712_20260718"

# Interferogram crops (data/nisar/<tag>.npz). Event pairs first.
PAIRS = {
    "172_008_A": "asc 07-08→07-20",
    "165_082_D": "desc 07-07→07-19",
    "100_A_008_STRADDLE": "asc 07-15→07-27",
    "100_A_008_CONTROL": "asc 07-03→07-15 (control)",
    "165_082_D_20260719": "desc 07-19→08-12 (post-event, 24 d)",
}

# The two geometries jointly inverted for the dike model.
INVERSION_PAIRS = ("172_008_A", "165_082_D")

# The core archive contains the event, control, continuation, Sentinel-1, and
# GeoTIFF data. The add-on contains the 15 pre-event NISAR crops.
CORE_NISAR_TAGS = (
    "172_008_A",
    "165_082_D",
    "100_A_008_STRADDLE",
    "100_A_008_CONTROL",
    "165_082_D_20260719",
)
PRE_EVENT_NISAR_TAGS = (
    "100_A_008_20260621",
    "165_082_D_20251109",
    "165_082_D_20251121",
    "165_082_D_20251203",
    "165_082_D_20251215",
    "165_082_D_20251227",
    "165_082_D_20260601",
    "165_082_D_20260625",
    "172_008_A_20251110",
    "172_008_A_20251122",
    "172_008_A_20251204",
    "172_008_A_20251216",
    "172_008_A_20251228",
    "172_008_A_20260602",
    "172_008_A_20260626",
)

# Probe points for the pre-event history (UTM 37N): the two co-event flank
# peaks and a southern patch checked as a precursor candidate.
PROBES = {
    "EAST": (687_480.0, 1_348_360.0),
    "WEST": (681_400.0, 1_347_400.0),
    "SOUTH": (688_120.0, 1_338_680.0),
}
