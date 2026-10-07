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
- `results/along_rift_halves.md` and `results/along_rift_halves.png`
  contain the north/south-half U∥ averages on the two flanks (the flank
  probes are local disks and a full-flank average cancels the
  antisymmetric along-rift pattern, so neither measures it), against the
  fitted model passed through the observations' far-field referencing
  and solve, with a correlated-noise floor from far-field control
  regions, per-differential single-cell influence, and a strike
  sensitivity.
- `mh26/gslc.py` rebuilds the co-event interferograms from the 10-m GSLC
  products (crops, interferogram, snaphu unwrapping); see
  `PROCESSING_NOTES.md`. The saved SMOOTH-cost sensitivity for the ascending
  model-input pair is in
  `results/gslc_unwrap_sensitivity_172_008_A.{json,png}`.
- `mh26/closure.py` runs the three-epoch loop-closure check of the
  reprocessed unwrapping per NISAR track (pre, co-event, and skip legs on
  the shared GSLC grid); `results/loop_closure.md`,
  `results/loop_closure_<tag>.json`, and `results/loop_closure_nisar.png`
  report the integer-cycle inconsistency statistics and maps. The nine
  legacy raw caches have unknown snaphu settings. The separate closure
  archive supplies the 18 inputs in `results/loop_closure_inventory.json`
  for `mh26.closure --check --cached-only`, which verifies the saved
  baseline and replays it without snaphu or file writes. Replacing reports
  after an analysis change requires the explicit `--regenerate` option.
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
  `results/distributed_fit_geometry_4a7b2386c8b5.{json,md}` (GUNW inputs;
  `--gslc` writes `..._gslc.{json,md}`, the chain the summary figures
  draw) and
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
- `mh26/hamling.py` decodes Ian Hamling's `dyke_files.zip` delivery (the
  September 2005 and fifteen later opening files on the ~74 km Hamling
  computational grid; the primary cumulative sums the canonical thirteen
  files of his 17 September 2026 update, which the data-derived version
  pairs corroborate, and the report records the superseded and
  untrimmed variants; the format is documented in the module docstring
  and the tracked tables live in `data/ihamling/opening_files/` with the
  17 September update in `data/ihamling/2026-09-17/`,
  see `data/README.md`), registers that grid on the Wright-plane s axis,
  and writes the three-epoch U⊥ comparison per input chain:
  `results/uperp_three_epoch{,_gslc}.{json,md,png,pdf}` (September 2005 —
  drawn from the updated 10-km model of the 17 September delivery, with
  the raw 20-km field tabled — cumulative 2006–2010, and July 2026 on
  one along-plane axis, with per-epoch volume and maximum-opening
  tables).
- `mh26/opening_profiles.py` renders the same inputs as top-10-km
  depth-averaged opening per intrusion in the style of Fig. 13 of
  Hamling et al. (2009): one labelled panel per event (September 2005,
  the canonical thirteen 2006–2010 intrusions, July 2026) plus the
  cumulative stack, per input chain:
  `results/opening_profiles{,_gslc}.{json,md,png,pdf}` and the
  standalone cumulative panel
  `results/opening_profiles_cumulative{,_gslc}.{png,pdf}` (sized to
  sit under the three-epoch figure). A profile integrates exactly to
  the epoch's above-10-km volume, cross-checked against `mh26.hamling`
  on every run; the delivered inversions carry no uncertainties, so
  the panels have no error bars.
- `mh26/opening_history.py` lays the three-epoch U⊥ sections and the
  cumulative depth-averaged stack out as one figure, per input chain:
  `results/opening_history{,_gslc}.{png,pdf}` (per-section colour
  scales) and `results/opening_history_common_scale{,_gslc}.{png,pdf}`
  (one shared scale topped at the September 2005 panel's 16.5 m, on
  which the July 2026 section reads mainly through its contour
  lines). Every drawn quantity comes from the `mh26.hamling` and
  `mh26.opening_profiles` reports; the module computes nothing new.
- `mh26/export_2026_opening.py` samples the fixed-plane 2026 fit at
  the patch centres of the Hamling computational grid and writes it in
  the delivery's `dyke_only_<epoch>_ex.txt` file convention (one value
  per geometry row, row order of `plot_3d_dyke_only.dat`), per input
  chain: `results/opening_2026_hamling_grid{,_gslc}.{txt,json}`. The
  sampled field integrates to 0.999 of the native fit volume; the
  JSON records the sampling rule, coverage, and that ratio.
- `mh26/fault_depth.py` sweeps the depth extent of the mapped east-fault
  term (truncation caps re-solved at the reused smoothing weight) and
  reports changes in misfit and fault potency/moment, per
  input chain: `results/fault_depth_mapped.{json,md,png}` and
  `results/fault_depth_mapped_gslc.{json,md,png}`. These are conditional
  sensitivities, not a data-required depth or a physical moment bound.
- `results/summary_*.png` (with matching PDF files) are eight draft
  publication-style figures built by `mh26.summary_figures` from the saved
  results: the four-geometry trace-supported fit with residuals and
  profiles, the free E–N–U decomposition with uncertainties and closure,
  the same solve in the dike frame (U⊥, U∥, up;
  `summary_decomposition_dike_frame`), the physical source geometry, the
  signed phase-gradient faulting evidence, the dike-induced stress on
  the east-fault receivers, and the fixed 2005-plane fit — map view
  (`summary_los_model_residual_wright2006`) and U⊥ cross-section
  (`summary_uperp_wright2006`).
  `summary_decomposition_dike_frame_pub` is the same dike-frame figure
  in manuscript form (`--only dike_frame_pub`): no in-file title,
  subtitle, or footnote, and no display-bookkeeping annotations — a
  manuscript caption carries them, and the counts stay in the result
  payloads; the weak-mask note and legends stay.
  `results/summary_figures.json` records the input and output hashes.
- `granules/` records the NISAR (GUNW + GSLC) and LiCSAR source products.
- `data/` contains archive manifests and instructions for the data files.

The core and s1tracks data archives support the standard-product analysis.
The gslc archive adds the 30-m co-event inputs, and the separate closure
archive supplies the raw arrays needed to reproduce loop closure (see
`data/README.md`). The pre-event archive is needed only
to regenerate the 16-window history from the NPZ layers. The core archive
already contains the corresponding GeoTIFFs.

## Current source model

The reference model is a vertical, uniform-opening rectangle in a homogeneous
elastic half-space. Its fitted parameters are:

- strike: 153.0°;
- length: 17.31 km;
- depth: 1.30–5.94 km;
- opening: 2.32 m;
- geometric volume: 0.187 km³;
- NISAR RMS residual: 2.91 cm.

The processing-choice runs span 0.186–0.208 km³. This range is not a
confidence interval. The separate depth and opening estimates are more
sensitive to processing choices than their product.

The model uses the GUNW valid-sample mask and scaled optimization variables.
The effects of these choices are described in `PROCESSING_NOTES.md`.

## Post-event deformation

For the NISAR descending 19 July–12 August pair, the projection onto the
co-event dike pattern is +0.107 with a processing-choice range of +0.100 to
+0.121.

The two-track Sentinel-1 tests confirm this signal and resolve it in time.
The early windows give +0.140 (ascending, 18–30 July) and +0.053 (descending,
23 July–4 August). The late windows (from 30 July and from 4 August) lie at
the level of the nine pre-event nulls (maximum magnitude 0.021). A
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
pattern; no dike-template deformation was detected after early August. The coefficients are still not converted
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
python -m mh26.along_rift
python -m mh26.distributed --gslc   # needs the gslc archive or a rebuild
python -m mh26.decompose3d --gslc
python -m mh26.along_rift --gslc
python -m mh26.faults
python -m mh26.distributed --mapped
python -m mh26.distributed --mapped --gslc
python -m mh26.fault_depth
python -m mh26.fault_depth --gslc
python -m mh26.wright2006
python -m mh26.distributed --geometry results/geometry_wright2006.json
python -m mh26.distributed --geometry results/geometry_wright2006.json --gslc
python -m mh26.wright2006 --panel results/distributed_fit_geometry_4a7b2386c8b5.json
python -m mh26.stress
python -m mh26.uncertainty --n 1000
python -m mh26.figures
python -m mh26.summary_figures   # needs the terrain tile (data/README.md)
python -m mh26.export_geotiff --all
python -m mh26.validate --profile full
```

With the closure data add-on installed, check the saved closure result with:

```bash
python -m mh26.closure --check --cached-only
python -m mh26.validate --profile full --require-closure-inputs
```

Run `python -m mh26.gslc --unwrap-sensitivity` only when a full-grid
cost-mode sensitivity is needed. It is not required to reproduce the primary
result, and each pair can take several minutes.

The Markdown reports and JSON files in `results/` contain the same numerical
outputs in human-readable and machine-readable forms.
