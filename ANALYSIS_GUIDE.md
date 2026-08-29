# Analysis guide

This repository contains NISAR displacement fields, a uniform-opening source
model, a Sentinel-1 comparison, a pre-event probe history, a post-event
pattern projection, and exact source-product identifiers for the July 2026
Manda Hararo intrusion.

## Package contents

- `results/los_panels.png` and `results/decomposition.png` show the NISAR LOS
  fields and the two-geometry quasi-East/quasi-Up decomposition.
- `results/okada_fit_panels.png` and `results/okada_fit.json` contain the
  uniform-opening dike fit.
- `results/model_sensitivity.md` reports the effect of nine processing choices.
- `results/s1_crosscheck_panels.png` compares the NISAR-derived model with a
  Sentinel-1 LiCSAR interferogram.
- `results/probe_history.md` reports 16 pre-event windows and the co-event and
  post-event windows.
- `results/continuation_panels.png` compares the 19 July–12 August descending
  pair with the co-event dike pattern.
- `results/s1_tracks.md` and its three `s1_tracks_*.png` figures contain the
  two-track Sentinel-1 tests: the descending co-event cross-check, six
  post-event projections with nine pre-event nulls, the decay fit, and the
  Sentinel-1-only East/Up decompositions.
- `results/distributed_fit.md`, `results/distributed_planes.png`, and
  `results/distributed_fit_panels.png` contain the distributed opening and
  graben-fault model fitted jointly to the two NISAR and two Sentinel-1
  co-event geometries. The mapped GSLC run also writes
  `results/distributed_fit_fringe_panels_mapped_gslc.png`, with rewrapped
  observations and model fields plus continuous residuals for all four looks.
- `results/decomposition_joint.md` and `results/decomposition_joint.png`
  contain the four-geometry decomposition (rift-perpendicular horizontal +
  vertical) with the cross-sensor closure map, the free three-component
  solve, and that solve rotated into the dike frame (U⊥, U∥, up) with a
  strike sensitivity.
- `mh26/gslc.py` rebuilds the co-event interferograms from the 10-m GSLC
  products (crops, interferogram, snaphu unwrapping); see
  `PROCESSING_NOTES.md`. The saved SMOOTH-cost sensitivity for the ascending
  model-input pair is in
  `results/gslc_unwrap_sensitivity_172_008_A.{json,png}`.
- `results/fault_traces.json` and `results/fault_traces.png` extract the
  phase-gradient ridges flanking the axis (candidate graben-fault scarps)
  with a LOS-step diagnostic per ridge; `mh26.distributed --mapped` uses only
  trace sides that pass that diagnostic (columns end at the dike plane;
  `--traces FILE` tests all traces in an external file and uses content-hash
  output names). `mh26.distributed --geometry FILE` loads a complete patch
  contract and records the source and geometry hashes. `mh26.wright2006`
  builds that contract from the 2005 dike-plane source model in
  `data/complex_geom2_revised_test.dat` (provenance and format in
  `data/README.md`); the fixed-plane 2026 fit is
  `results/distributed_fit_geometry_b643c0de55aa.{json,md}` and
  `results/uperp_2026_wright2006.png` is one 2026 U⊥ cross-section on
  the Wright et al. (2006) plane. It uses the along-strike/depth view of an
  individual U⊥ opening panel in Hamling et al. (2010) Figure 2; it does not
  reproduce that figure's multi-event opening-and-stress layout. The panel
  uses a fixed 0–2.5 m scale and labels the added southeast section. Its
  0–10 km plane is not the later Hamling computational grid.
  `results/stress_change.md` resolves the dike's stress change as spatial
  maps and a depth summary on the supported fault surface;
  `results/distributed_uncertainty.md` gives
  correlated-noise Monte Carlo intervals for the distributed model
  (pivotal, about the noiseless refit); `results/fringes_model_comparison.png`
  compares rewrapped observed and modelled fringes with a cross-axis
  profile.
- `results/summary_*.png` (with matching PDF files) are six draft
  publication-style figures built by `mh26.summary_figures` from the saved
  results: the four-geometry trace-supported fit with residuals and
  profiles, the free E–N–U decomposition with uncertainties and closure,
  the same solve in the dike frame (U⊥, U∥, up;
  `summary_decomposition_dike_frame`), the physical source geometry, the
  signed phase-gradient faulting evidence, and the dike-induced stress on
  the east-fault receivers.
  `results/summary_figures.json` records the input and output hashes.
- `granules/` records the NISAR (GUNW + GSLC) and LiCSAR source products.
- `data/` contains archive manifests and instructions for the data files.

The core and s1tracks data archives together support all reported results
except the GSLC reprocessing outputs, which the gslc archive provides (or a
rebuild from ASF, see `data/README.md`). The pre-event archive is needed only
to regenerate the 16-window history from the NPZ layers. The core archive
already contains the corresponding GeoTIFFs.

## Current source model

The reference model is a vertical, uniform-opening rectangle in a homogeneous
elastic half-space. Its fitted parameters are:

- strike: 152.6°;
- length: 17.15 km;
- depth: 1.49–4.99 km;
- opening: 2.81 m;
- geometric volume: 0.169 km³;
- NISAR RMS residual: 3.05 cm.

The processing-choice runs span 0.167–0.188 km³. This range is not a
confidence interval. The separate depth and opening estimates are more
sensitive to processing choices than their product.

The model uses the GUNW valid-sample mask and scaled optimization variables.
The effects of these choices are described in `PROCESSING_NOTES.md`.

## Post-event deformation

For the NISAR descending 19 July–12 August pair, the projection onto the
co-event dike pattern is +0.107 with a processing-choice range of +0.100 to
+0.123.

The two-track Sentinel-1 tests confirm this signal and resolve it in time.
The early windows give +0.140 (ascending, 18–30 July) and +0.053 (descending,
23 July–4 August). The late windows (from 30 July and from 4 August) lie at
the level of the nine pre-event nulls (maximum magnitude 0.022). A
window-integral exponential fit gives a total of about 0.15 of the co-event
pattern with a decay time of about 4–5 days from 18 July (point fit 4.3 days;
the one fit uses all six windows, is descriptive, has no formal uncertainty,
and weights overlapping windows equally). The point-fit curve predicts +0.120
for the independent NISAR
window, against +0.107 observed. The Sentinel-1-only East/Up decomposition of
the 24-day post-event pairs shows +18 cm of flank divergence, about one sixth
of the co-event value; the ascending and descending windows are offset by
five days at each end, so the decomposition is illustrative and the
per-geometry projections are the quantitative result (see `LIMITATIONS.md`).

The deformation therefore continued after the main intrusion in a dike-like
pattern and stopped by early August. The coefficients are still not converted
to an added-volume estimate; that requires source modelling of the post-event
signal itself.

## Scope

- Long-term regional InSAR time series are not included.
- Local seismic waveform and catalogue data are not included.
- The reference uniform model does not include fault slip, topography,
  elastic layering, or a coupled reservoir. The distributed model adds graben
  fault slip on a geometry fixed a priori (see `LIMITATIONS.md`); topography,
  layering, and reservoir coupling remain unmodelled. The post-event signal
  is not source-modelled.
- No seismic-to-geodetic moment ratio is calculated.
- No relation to other Afar events is tested.
- Public reuse terms have not been selected.

See `LIMITATIONS.md` before interpreting the model parameters or the
post-event projection.

## Suggested processing order

```bash
python -m unittest discover -s tests -v
python -m mh26.validate --profile core
python -m mh26.los_fields
python -m mh26.invert
python -m mh26.sensitivity
python -m mh26.s1_check
python -m mh26.history
python -m mh26.continuation
python -m mh26.s1_tracks
python -m mh26.distributed          # primary fit; selects the smoothing weight
python -m mh26.decompose3d
python -m mh26.distributed --gslc   # needs the gslc archive or a rebuild
python -m mh26.decompose3d --gslc
python -m mh26.faults
python -m mh26.distributed --mapped
python -m mh26.distributed --mapped --gslc
python -m mh26.distributed --geometry FILE  # after an external contract arrives
python -m mh26.stress
python -m mh26.uncertainty --n 1000
python -m mh26.figures
python -m mh26.summary_figures   # needs the terrain tile (data/README.md)
```

Run `python -m mh26.gslc --unwrap-sensitivity` only when a full-grid
cost-mode sensitivity is needed. It is not required to reproduce the primary
result, and each pair can take several minutes.

The Markdown reports and JSON files in `results/` contain the same numerical
outputs in human-readable and machine-readable forms.
