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

# LiCSAR frames with local East/North/Up LOS-vector rasters. 014A is
# ascending; 079D is descending (full two-track set released 2026-08-17).
S1_FRAMES = ("014A_07688_131313", "079D_07694_131313")

# All LiCSAR pairs used by the two-track analysis: (frame, pair, role).
# Pre-event pairs provide an empirical null for the dike-pattern projection.
S1_TRACK_PAIRS = (
    ("014A_07688_131313", "20260712_20260718", "co-event"),
    ("079D_07694_131313", "20260711_20260723", "co-event"),
    ("014A_07688_131313", "20260718_20260730", "post-event"),
    ("014A_07688_131313", "20260730_20260811", "post-event"),
    ("014A_07688_131313", "20260718_20260811", "post-event"),
    ("079D_07694_131313", "20260723_20260804", "post-event"),
    ("079D_07694_131313", "20260804_20260816", "post-event"),
    ("079D_07694_131313", "20260723_20260816", "post-event"),
    ("014A_07688_131313", "20260605_20260612", "pre-event"),
    ("014A_07688_131313", "20260612_20260617", "pre-event"),
    ("014A_07688_131313", "20260617_20260624", "pre-event"),
    ("014A_07688_131313", "20260624_20260629", "pre-event"),
    ("014A_07688_131313", "20260629_20260706", "pre-event"),
    ("014A_07688_131313", "20260706_20260712", "pre-event"),
    ("079D_07694_131313", "20260610_20260622", "pre-event"),
    ("079D_07694_131313", "20260622_20260629", "pre-event"),
    ("079D_07694_131313", "20260629_20260711", "pre-event"),
)

# The two 24-day pairs used for the post-event two-geometry decomposition and
# the co-event pairs decomposed for comparison.
S1_DECOMPOSITION_SETS = {
    "co-event": (
        ("014A_07688_131313", "20260712_20260718"),
        ("079D_07694_131313", "20260711_20260723"),
    ),
    "post-event": (
        ("014A_07688_131313", "20260718_20260811"),
        ("079D_07694_131313", "20260723_20260816"),
    ),
}

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
