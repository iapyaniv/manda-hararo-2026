# Manda Hararo 2026

This repository contains reproducible NISAR observations, source modelling,
timing constraints, and a Sentinel-1 comparison for the July 2026 dike
intrusion at the Manda Hararo rift segment in Afar, Ethiopia.

See [`ANALYSIS_GUIDE.md`](ANALYSIS_GUIDE.md) for the contents and processing
order. See [`LIMITATIONS.md`](LIMITATIONS.md) before interpreting the source
parameters or the post-event result. Processing details that affect the
numerical results are in [`PROCESSING_NOTES.md`](PROCESSING_NOTES.md).

Maintainer: Vinay Pai <vinay@yowon.org>.

## Results

- Co-event LOS maxima are +76.6 cm ascending and +58.2 cm descending.
  See `results/los_summary.md`.
- The onset is after 15 July and present by 18 July, based on the NISAR
  control and LiCSAR event pairs.
- The uniform dike has strike 152.6°, length 17.15 km, depth 1.49–4.99 km,
  opening 2.81 m, and RMS 3.05 cm. See `results/okada_fit.md`.
- Processing cases span 17.1–17.5 km length, 1.3–6.0 km depth, and
  2.31–2.92 m opening. Geometric volume is 0.169 km³, with a
  processing-choice range of 0.167–0.188 km³. See
  `results/model_sensitivity.md`.
- The Sentinel-1 ascending comparison gives r = 0.954, amplitude ratio 0.88,
  and RMS 2.3 cm. The descending comparison gives r = 0.917, amplitude ratio
  0.85, and RMS 2.8 cm. See `results/s1_crosscheck.md` and
  `results/s1_tracks.md`.
- The NISAR post-event projection is +0.107 for 19 July–12 August, with a
  processing-choice range of +0.100–+0.123. See `results/continuation.md`.
- Sentinel-1 post-event projections are +0.140 ascending for 18–30 July and
  +0.053 descending for 23 July–4 August. Later windows are at the null
  level. One descriptive fit across all six windows gives a total near 0.15
  of the co-event pattern and a 4.3-day decay time from 18 July. It has no
  formal uncertainty. See `results/s1_tracks.md`.
- The illustrative Sentinel-1-only East/Up decomposition gives +18 cm flank
  divergence. Its two windows are offset by five days. See
  `LIMITATIONS.md` and `results/s1_tracks.md`.

## Model extensions and checks

- The trace-supported distributed fit includes only the east fault. It
  reduces joint RMS from 2.59 to 2.42 cm and improves both 5-km holdout
  splits. The 2.30 cm two-fault idealization remains a geometry sensitivity.
  NISAR→Sentinel-1 is the stronger held-out-sensor transfer direction. See
  `results/distributed_fit_mapped.md` and `LIMITATIONS.md`.
- The joint decomposition, free three-component variant, closure diagnostics,
  and 30-m-input repeat are in `results/decomposition_joint*.md`. Counts now
  distinguish exactly three geometries from four geometries.
- The free three-component solve is also rotated into the dike frame: U⊥
  perpendicular to the fitted strike (ENE positive, the same axis as the
  constrained solve), U∥ along strike (SSE positive), and up. The rotation
  is exact and the rotated formal sigmas carry the east–north
  cross-covariance. U∥ is the weakly observed component: with the 30-m
  inputs its flank probes are +5.7 cm (east) and −4.6 cm (west) against a
  median formal σ of 4.1 cm; the GUNW inputs give +10.3 and −3.9 cm, and
  replacing the fitted strike by the southernmost 2005 segment strike (4°
  different) moves them by about 4 cm. U∥ is therefore reported but not
  resolved. See the dike-frame section of `results/decomposition_joint*.md`.
- The GSLC reprocessing comparison and its near-axis differences are in
  `results/gslc_gunw_compare_*.json`.
- A single-tile SMOOTH/MCF sensitivity changes 6.1% of near-axis pixels by
  more than half a fringe in the ascending model-input pair. It agrees
  closely off axis. This is not an independent unwrap validation. See
  `results/gslc_unwrap_sensitivity_172_008_A.json`.
- Phase-gradient ridges, threshold sensitivity, and an along-strike alignment
  null are in `results/fault_traces.json`. Only the east ridge has consistent
  LOS-step support. The default trace-column fit excludes the unsupported west
  ridge.
- The trace-column fit is a geometry sensitivity, not a fault map. Its source
  path and SHA-256 hash are stored in `results/distributed_fit_mapped.json`.
- The stress report uses the source and receiver geometry stored in that fit.
  It includes spatial stress maps and a depth summary. It does not use the
  uniform rectangle's fitted top. See `results/stress_change.md`.
- The Monte Carlo applies the far-field datum and ramp projection to the model
  and to every covariance draw. See `results/distributed_uncertainty.md`.
- Same-chain rewrapped observations, model fields, and residuals for all four
  geometries are in `results/distributed_fit_fringe_panels_mapped_gslc.png`.
  The separate GUNW-fit-to-GSLC check is in
  `results/fringes_model_comparison.png`.
- Eight draft publication-style figures summarize the trace-supported fit,
  the free E–N–U decomposition, the same solve in the dike frame (U⊥, U∥,
  up), the physical source, the phase-gradient faulting evidence, the
  dike-induced stress on the east-fault receivers, and the fit on the
  fixed 2005 plane (map view and U⊥ cross-section): `results/summary_*.png`
  with matching vector PDF files and the
  identity manifest `results/summary_figures.json`
  (`mh26.summary_figures`; needs the terrain tile described in
  `data/README.md`).
- `mh26.distributed --geometry FILE` accepts a complete saved geometry,
  validates its patch contract, records file and geometry hashes, and uses
  hash-suffixed outputs.
- `mh26.wright2006` builds such a geometry from the 2005 dike-plane source
  model of Wright et al. (2006) in `data/complex_geom2_revised_test.dat`
  (five chained vertical segments, 61 km, 0–10 km depth; see
  `data/README.md`). The 2026 uniform-fit dike is aligned with the
  southernmost segment (strike differs by ~4°; endpoints 0.35–0.88 km off
  its line) and overruns its SE end by ~3 km, so the plane is extended 10 km
  SE. The 2026 joint inversion on that fixed plane gives joint RMS 2.46 cm
  (trace-supported free-trace fit: 2.42 cm), maximum opening 2.44 m, and
  volume 0.156 km³ (`results/distributed_fit_geometry_b643c0de55aa.md`);
  on the reprocessed 30-m (GSLC) inputs the same comparison is 2.81 cm
  against 2.80 cm, with maximum opening 2.42 m and volume 0.154 km³
  (`results/distributed_fit_geometry_b643c0de55aa_gslc.md`).
  `results/uperp_2026_wright2006.png` is one 2026 U⊥ cross-section on
  the Wright plane, from the GUNW-input fit (a `--gslc` fit writes a
  `_gslc` name). It uses the along-strike/depth view of an individual
  opening panel in Hamling et al. (2010) Figure 2; it is not a reproduction
  of that multi-event opening-and-stress figure. Its fixed 0–2.5 m scale can
  support a later comparison after the historical coordinate definitions
  are checked. The panel labels the 10-km southeast extension and the
  0–10 km model depth. The historical opening distributions are still
  awaited. This is a fixed-historical-plane variant with nominally
  transferred smoothing, not a co-planarity result; see `LIMITATIONS.md`.
  The publication-style versions of the fixed-plane results are
  `results/summary_los_model_residual_wright2006.png` (observed, model,
  residual, and profile, in the layout of the free-trace figure) and
  `results/summary_uperp_wright2006.png` (the cross-section restyled),
  both from the 30-m-input fit.

Sixteen pre-event windows from November 2025 through 15 July 2026 have a
maximum absolute 1-km probe mean of 3.4 cm (`results/probe_history.md`). Nine
pre-event Sentinel-1 pairs bound the post-event projection null at ±0.022
(`results/s1_tracks.md`). The post-event signal appears in both Sentinel-1
geometries and in NISAR, and it decays; no added-volume estimate is made.

## Reproduction

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
# Unpack the data archives as described in data/README.md.

python -m mh26.okada
python -m mh26.los_fields
python -m mh26.invert
python -m mh26.sensitivity
python -m mh26.s1_check
python -m mh26.history
python -m mh26.continuation
python -m mh26.s1_tracks
python -m mh26.distributed            # primary fit; selects the smoothing weight
python -m mh26.decompose3d
python -m mh26.distributed --gslc     # needs the gslc data archive (or a rebuild)
python -m mh26.decompose3d --gslc
python -m mh26.faults                 # phase-gradient ridges (needs gslc archive)
python -m mh26.distributed --mapped   # ridge-anchored fault-geometry variant
python -m mh26.distributed --mapped --gslc
python -m mh26.distributed --geometry FILE  # complete external patch geometry
python -m mh26.wright2006             # 2005 dike-plane contract from data/*.dat
python -m mh26.distributed --geometry results/geometry_wright2006.json
python -m mh26.distributed --geometry results/geometry_wright2006.json --gslc
python -m mh26.wright2006 --panel results/distributed_fit_geometry_b643c0de55aa.json
python -m mh26.stress
python -m mh26.uncertainty --n 1000
python -m mh26.figures
python -m mh26.summary_figures        # needs the terrain tile (data/README.md)
python -m mh26.export_geotiff --all
python -m mh26.validate --profile full
```

The shipped `data/gslc/*_ifg.npz` and `*_unw.npz` reproduce from the raw
GSLC granules with `python -m mh26.gslc --crop --ifg --unwrap --compare`
(Earthdata Login in `~/.netrc`; `snaphu` on the PATH, e.g.
`apt install snaphu`; tens of GB of ranged transfer).
Run `python -m mh26.gslc --unwrap-sensitivity` for the documented
single-tile SMOOTH/MCF sensitivity on the two model-input pairs. Each full
grid can take several minutes.

Run all commands from the repository root. `requirements-lock.txt` records the
environment used for the committed outputs. `requirements.txt` gives broader
minimum versions. The committed outputs provide a reference for comparison
with a new run. The sensitivity, Sentinel-1, continuation, and two-track
steps read the model written by `mh26.invert`.

To rebuild the NISAR crops from ASF, configure an Earthdata Login in
`~/.netrc` and run `python -m mh26.fetch`. Exact product identifiers are in
`granules/pairs.tsv`.

## Repository map

```text
mh26/params.py          shared constants and data groups
mh26/gunw.py            screening, corrections, referencing, and LOS geometry
mh26/okada.py           tensile-dislocation forward model and sign test
mh26/invert.py          scaled multistart NISAR inversion
mh26/sensitivity.py     processing-choice model runs
mh26/s1_check.py        Sentinel-1/LiCSAR comparison
mh26/continuation.py    19 July–12 August pattern projection
mh26/s1_tracks.py       two-track Sentinel-1 post-event and null tests
mh26/distributed.py     distributed opening + graben-fault inversion (4 geometries)
mh26/decompose3d.py     joint NISAR+Sentinel-1 co-event decomposition
mh26/gslc.py            full-resolution GSLC interferograms and unwrapping
mh26/faults.py          graben fault traces from wrapped-phase gradients
mh26/stress.py          dike-induced stress change on the fault surfaces
mh26/uncertainty.py     correlated-noise Monte Carlo for the distributed model
mh26/figures.py         rewrapped-fringe data/model comparison figure
mh26/summary_figures.py publication-style summary figures from saved results
mh26/fetch.py           ASF GUNW crop extraction
mh26/validate.py        data, source-manifest, result, and checksum validation
scripts/build_release.py deterministic code and data archive builder
granules/               NISAR and LiCSAR source manifests
results/                numerical reports and figures
data/                   archive layout and verification guide
```

## Processing conventions

- NISAR LOS displacement is `−phase × 0.2384 / (4π)` metres, positive toward
  the satellite. The supplied ionosphere screen is subtracted.
- Pixels require coherence greater than 0.35, a real connected component, and
  valid reference and secondary samples in the GUNW valid-sample mask. Model
  residuals are not used to reject pixels or blocks.
- A plane is fitted beyond 28 km from 12.18°N, 40.66°E and removed from each
  pair. NISAR grids are UTM 37N (EPSG:32637), normally at 80 m posting.
- Local product LOS vectors are interpolated at 500 m ellipsoidal height. The
  two-geometry decompositions fix north displacement to zero and are reported
  as quasi-East and quasi-Up displacement; the four-geometry decomposition
  instead constrains the horizontal direction to rift-perpendicular, and its
  free three-component variant is also reported in the dike frame
  (`results/decomposition_joint.md`).
- The source is a vertical, uniform-opening rectangle in a homogeneous elastic
  half-space with λ = μ = 30 GPa. `μV` is a geodetic comparison scale, not a
  measured seismic moment.

## Collaborator-supplied Sentinel-1 time-series maps

`Sentinel_1_Cumulative_Disp/` holds four cumulative LOS displacement maps
from the Sentinel-1 time series of Pedro Espin Bedon (Leeds), committed by
Ben Ireland (Bristol) on 28 August 2026 for the pre-dike context panels:
one pre-dike and one dike-spanning map per LiCSAR frame (014A:
2014-10-30→2026-07-06 and 2014-10-06→2026-08-11; 079D: 2014-10-23→2026-07-11
and 2014-10-23→2026-08-16). They are EPSG:4326 GeoTIFFs at ≈0.00054°
posting, float32 millimetres, NaN where empty, without a nodata tag. They
are outside the data archives and checksum manifests, and no `mh26` result
uses them. One empirical check was made: the dike-spanning minus pre-dike
difference of each frame correlates with this repository's referenced
LiCSAR co-event field at r = 0.97 (014A) and 0.99 (079D), with slopes 1.22
and 1.10 in the positive-toward-satellite sense used here, so the two
products share a sign convention, and the excess over 1 is consistent with
the post-event tail inside those longer windows. The sign, the reference area, the
temporal reference (the 014A pair starts 24 days apart: 2014-10-06 and
2014-10-30, so their difference contains a short unmatched interval), and
the corrections still need confirmation from the processors before the maps
are used in a figure.

## Data sources and reuse

NISAR inputs are L2 GUNW Provisional V1 products distributed by ASF DAAC. See
the [NISAR GUNW guide](https://nisar-docs.asf.alaska.edu/gunw/) and the exact
granule manifest. Sentinel-1 inputs are Copernicus-derived COMET-LiCSAR
products. Follow the [LiCSAR portal use conditions, acknowledgement, and
citations](https://comet.nerc.ac.uk/COMET-LiCS-portal/). This repository does
not label LiCSAR products as CC-BY.

Analysis code and derived content are © 2026 Yowon Labs. No public licence is
provided. Obtain permission before public redistribution.

## Release assets

`scripts/build_release.py` builds the archives deterministically: entry
timestamps and the archive comment come from the build commit, so a rebuild
at a new commit changes every archive hash even when the packed data is
unchanged.

- Release `analysis-2026-08-18` (commit `f9799dd`): code, core, and
  pre-event archives with their checksum file.
- Release `analysis-2026-08-20`: revised code archive, the two-track
  `s1tracks` add-on, and the `gslc` add-on with the 30-m reprocessed
  interferograms and unwrapped fields. The 2026-08-18 **core** and
  **pre-event** archives are unchanged, are not re-uploaded, and keep their
  published checksums:

      b245d47e063eb278cef203e4305a93d54ed72f721fda1e834f67eb2a51861e1b  manda-hararo-2026_data-core.zip
      f6dc9271b5aaae87300687ad7d729fde33d13bc1a127eadf7a7e350ba6a45a8f  manda-hararo-2026_data-preevent.zip

From `analysis-2026-08-20` on, each release ships a
`SHA256SUMS-release-<tag>.txt` manifest listing the exact hash of every
asset the release refers to, including retained assets from earlier
releases. A locally generated `dist/SHA256SUMS.txt` describes a full rebuild
at the current commit only; it is not uploaded, because its core/pre-event
lines cannot match retained assets built at an earlier commit.
