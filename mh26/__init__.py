"""Analysis of the July 2026 Manda Hararo (Afar) dike intrusion in NISAR data.

Processing order (run from the repository root; the ordered command list
with the required data archives is in README.md and ANALYSIS_GUIDE.md):
    fetch (or unpack the data archives), okada, los_fields, invert,
    sensitivity, s1_check, history, continuation, s1_tracks,
    distributed, decompose3d, gslc, faults, distributed --mapped,
    stress, uncertainty, figures, export_geotiff, validate.
"""

__version__ = "0.2.0"
