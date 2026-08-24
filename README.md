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

| quantity | result | output |
|---|---|---|
| co-event LOS | maximum +76.6 cm ascending and +58.2 cm descending | `results/los_summary.md` |
| onset | after 15 July and present by 18 July | NISAR control and LiCSAR event pairs |
| uniform dike | strike 152.6°, length 17.15 km, depth 1.49–4.99 km, opening 2.81 m; RMS 3.05 cm | `results/okada_fit.md` |
| processing sensitivity | length 17.1–17.5 km, depth 1.3–6.0 km, opening 2.31–2.92 m | `results/model_sensitivity.md` |
| geometric volume | 0.169 km³; processing-choice range 0.167–0.188 km³ | `results/model_sensitivity.md` |
| Sentinel-1 ascending comparison | r = 0.954, amplitude ratio 0.88, RMS 2.3 cm | `results/s1_crosscheck.md` |
| Sentinel-1 descending comparison | r = 0.917, amplitude ratio 0.85, RMS 2.8 cm | `results/s1_tracks.md` |
| NISAR post-event projection | +0.107 for 19 July–12 August; processing-choice range +0.100–+0.123 | `results/continuation.md` |
| Sentinel-1 post-event projections | +0.140 (ascending 18–30 July) and +0.053 (descending 23 July–4 August); later windows at the null level | `results/s1_tracks.md` |
| post-event decay | one fit across six windows gives a total ≈0.15 of the co-event pattern and a 4.3-day decay time from 18 July; this is a descriptive fit with no formal uncertainty | `results/s1_tracks.md` |
| post-event opening geometry | +18 cm flank divergence in the Sentinel-1-only East/Up decomposition (illustrative: its two windows are offset by five days; see `LIMITATIONS.md`) | `results/s1_tracks.md` |

## Model extensions and checks

- The distributed opening and fault-basis fit reduces joint RMS from 2.59 to
  2.30 cm. NISAR→Sentinel-1 is the held-out-sensor transfer test. The reverse
  transfer is conditional because NISAR affected geometry selection. See
  `results/distributed_fit.md`.
- The joint decomposition, free three-component variant, closure diagnostics,
  and 30-m-input repeat are in `results/decomposition_joint*.md`. Counts now
  distinguish exactly three geometries from four geometries.
- The GSLC reprocessing comparison and its near-axis differences are in
  `results/gslc_gunw_compare_*.json`.
- Phase-gradient ridges, threshold sensitivity, and an along-strike alignment
  null are in `results/fault_traces.json`. Only the east ridge has consistent
  LOS-step support. The default trace-column fit excludes the unsupported west
  ridge.
- The trace-column fit is a geometry sensitivity, not a fault map. Its source
  path and SHA-256 hash are stored in `results/distributed_fit_mapped.json`.
- The stress report uses the source and receiver geometry stored in that fit.
  It reports stress by receiver depth and does not use the uniform rectangle's
  fitted top. See `results/stress_change.md`.
- The Monte Carlo applies the far-field datum and ramp projection to the model
  and to every covariance draw. See `results/distributed_uncertainty.md`.
- The rewrapped fringe comparison is in
  `results/fringes_model_comparison.png`.

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
python -m mh26.stress
python -m mh26.uncertainty --n 1000
python -m mh26.figures
python -m mh26.export_geotiff --all
python -m mh26.validate --profile full
```

The shipped `data/gslc/*_ifg.npz` and `*_unw.npz` reproduce from the raw
GSLC granules with `python -m mh26.gslc --crop --ifg --unwrap --compare`
(Earthdata Login in `~/.netrc`; `snaphu` on the PATH, e.g.
`apt install snaphu`; tens of GB of ranged transfer).

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
  instead constrains the horizontal direction to rift-perpendicular
  (`results/decomposition_joint.md`).
- The source is a vertical, uniform-opening rectangle in a homogeneous elastic
  half-space with λ = μ = 30 GPa. `μV` is a geodetic comparison scale, not a
  measured seismic moment.

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
