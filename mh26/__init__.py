"""Analysis of the July 2026 Manda Hararo (Afar) dike intrusion in NISAR data.

Pipeline (run from the repository root, see README.md):
    python -m mh26.fetch           # step 0: rebuild data/nisar
    python -m mh26.okada           # sign-convention selftest
    python -m mh26.los_fields      # step 1: LOS panels + E/U decomposition
    python -m mh26.invert          # step 2: Okada dike inversion
    python -m mh26.sensitivity     # step 3: processing-choice model sensitivity
    python -m mh26.s1_check        # step 4: Sentinel-1 cross-check
    python -m mh26.history         # step 5: pre/co/post-event probe history
    python -m mh26.continuation    # step 6: post-event dike-pattern test
    python -m mh26.validate        # data and result validation
    python -m mh26.export_geotiff  # utility: GeoTIFFs for GIS use
"""

__version__ = "0.2.0"
