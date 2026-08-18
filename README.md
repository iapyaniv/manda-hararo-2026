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
| Sentinel-1 comparison | r = 0.954, amplitude ratio 0.88, RMS 2.3 cm | `results/s1_crosscheck.md` |
| post-event pattern projection | +0.107 for 19 July–12 August; processing-choice range +0.100–+0.123 | `results/continuation.md` |

Sixteen pre-event windows from November 2025 through 15 July 2026 have a
maximum absolute 1-km probe mean of 3.4 cm (`results/probe_history.md`). The
post-event projection uses one 24-day descending geometry and is not an
estimate of added volume.

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
python -m mh26.export_geotiff --all
python -m mh26.validate --profile full
```

Run all commands from the repository root. `requirements-lock.txt` records the
environment used for the committed outputs. `requirements.txt` gives broader
minimum versions. The committed outputs provide a reference for comparison
with a new run. The sensitivity, Sentinel-1, and continuation steps read the
model written by `mh26.invert`.

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
  two-geometry decomposition fixes north displacement to zero and is reported
  as quasi-East and quasi-Up displacement.
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
