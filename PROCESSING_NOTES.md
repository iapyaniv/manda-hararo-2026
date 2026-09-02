# Processing notes

This file records processing details that affect the numerical results.

## GUNW extraction

The GUNW HDF5 coherence dataset is named `coherenceMagnitude`. The extraction
code maps source dataset names to the shorter NPZ keys and checks the required
datasets, raster shapes, projection, metadata, and output file before use.

## Valid-sample mask

Coherence and connected-component screens do not remove every sample that the
GUNW product marks as invalid. All analysis paths and GeoTIFF exports therefore
apply the reference and secondary valid-sample flags before further screening.

Without this mask, the main ascending fit pair contains 21,183 invalid pixels
(2.8% of the earlier analysis mask). Each track-100 pair contains about 45,000
such pixels (7.3–7.6%). The main descending fit pair contains none. These
pixels produce a diagonal subswath-edge feature in the ascending field.

No block is rejected on the basis of model residuals.

## Inversion scaling

The inversion contains kilometre-scale position and length parameters and
metre-scale opening and offset parameters. Optimization uses dimensionless
scaled variables and explicit finite-difference steps. All six starting points
converge to the same solution basin. The best optimality is 0.0043, and the
maximum length difference among the six solutions is less than 1 m.

## Sentinel-1 coherence

The local LiCSAR coherence raster uses UInt8 values from 0 to 255. The
comparison normalizes this raster to 0–1 before it applies the 0.25 threshold.
It also requires identical grids for all five input rasters and checks the LOS
vector norms.

## Two-track Sentinel-1 processing

The two-track module reuses the pixel statistic of `mh26.s1_check` for the
co-event cross-checks and the block statistic of `mh26.continuation` for the
projections: 12-pixel (≈1.33 km) blocks on the 0.001° LiCSAR grid, a 25-km
fit radius, a 1-cm prediction threshold, and a fitted intercept.

The LiCSAR unwrapped-phase sign is resolved empirically once per frame from
the co-event regression against the model (both frames resolve to −1) and is
then applied unchanged to that frame's pre- and post-event pairs. A weak pair
must not resolve its own sign; that would bias every projection positive.
The validator refuses a sign change.

The pre-event pairs give an empirical null distribution for the projection
coefficient. Atmospheric and ionospheric noise enters the null pairs in the
same way as the post-event pairs.

The decay fit is a least-squares fit of
`A (exp(-t1/tau) - exp(-t2/tau))` to the six Sentinel-1 window-integrated
projections with t in days after 18 July. The East/Up decomposition bins each
geometry's screened, sign-aligned LOS and LOS vectors onto a shared 1.28-km
UTM lattice (at least 40 samples per cell and geometry) and solves the
zero-north two-component system per cell.

The 2026-08-17 LiCSAR regeneration republished the ascending co-event pair
byte-identically, so the core-archive copies match the current portal files.

## Distributed model (mh26.distributed)

The distributed inversion fixes all geometry a priori (dike plane and strike
from the uniform fit; graben faults ±2 km from the trace, dip 65° inward,
0–4 km depth) and solves only dislocation magnitudes: patch opening on the
dike (28 km × 10 km, 1.75 × 1 km patches) and normal slip on the two faults,
non-negative, with one Laplacian smoothing weight (selected as described
below) and a per-dataset offset. Data are the two NISAR geometries (mh26.invert blocks)
plus the two LiCSAR co-event pairs (mh26.s1_tracks blocks), each divided by
its decay-normalization factor (window end vs 18 July; factors 1.000–1.105)
and weighted by its far-field robust scatter.

The module selftests Green's-function conventions before solving: a patch
grid tiling the uniform rectangle reproduces `okada.forward` to < 0.1 mm,
and positive normal slip on the east fault subsides the axis side.

The smoothing weight is selected by held-out prediction, not by the
L-curve: on this problem the discrete L-curve curvature maximum is not
bracketed by any moderate grid (it keeps climbing into the heavily
over-damped branch, where held-out prediction is clearly worse), so a
curvature "knee" is a boundary artifact and is reported for reference
only. The criterion is the mean 10-km checkerboard held-out RMS over
λ = 0.03–1000 (19 logarithmic steps); it has an interior minimum at
λ = 17.8, selection raises an error at a grid edge, and the full selection
table plus a cross-sensor transfer check are stored with the primary
result. Every variant (mapped, GSLC) reuses the primary λ so comparisons
are at matched regularization, and a λ×1/3.2–λ×3.2 envelope is carried in
the sensitivity list so the λ-dependence of each quoted quantity is
visible. The 10-km cells are not fully independent under 15–30 km noise
correlation lengths. NISAR→Sentinel-1 is therefore reported as the stronger
held-out-sensor check. S1→NISAR is only a conditional reverse transfer because
the NISAR fields affected geometry selection.

Joint RMS on the four geometries: uniform rectangle 2.49 cm → distributed
opening 2.37 cm → opening + graben faults 2.09 cm. The graben faults carry
about a fifth of the dike's potency (Σ area·dislocation; ratio 0.196,
not a moment-tensor ratio). The dike volume
of the distributed model (0.170 km³; 0.164–0.187 across the smoothing
envelope, geometry variants, sensor subsets, and the
no-decay-normalization variant) is not the same quantity as the
uniform-rectangle volume (0.187 km³): part of the near-axis signal moves
from opening into fault slip. A volume is meaningful only together with
its model class.

## Joint decomposition (mh26.decompose3d)

The four-geometry decomposition solves, per 1.28-km lattice cell, horizontal
motion along the rift-perpendicular azimuth (strike − 90°) plus vertical.
This replaces the zero-north assumption of the earlier two-geometry
decompositions: near-polar LOS vectors barely observe north, and zero-north
aliases real rift-normal motion. Cells need at least two geometries (exact
solve); with three or more, the weighted residual is a cross-sensor closure
diagnostic (median 0.55 σ). Inputs are decay-normalized as in the
distributed inversion.

The free three-component solve (cells with at least three looks) is also
rotated into the dike frame of `mh26.distributed._strike_vectors`: U∥
along the fitted strike (153.0°, positive toward 333.0° — the near-north
direction, so northward motion draws warm), U⊥ perpendicular-left (ENE
positive; the same axis as the constrained h), and up. A linear
reparameterization leaves the least-squares solution unchanged, so
coverage and closure are identical; the rotated formal variances include
the east–north cross-covariance, which the near-polar looks make large, and
the module checks that the rotated U⊥ flank probes equal the projected
rift-perpendicular probes. U∥ is nearly north–south and inherits the north
uncertainty (median σ 4.1 cm against 2.0–2.1 cm for U⊥ and 0.8 cm for up). Two
saved sensitivities bound it: the input chain (GUNW versus 30-m fields)
moves the flank U∥ probes by several cm, and using the southernmost Wright
et al. (2006) segment strike (148.5°, 4.5° from the fit) moves them by about
5 cm, because a strike error of δ moves U⊥·sin δ into U∥ (≈5 cm for 60 cm
of U⊥). The east-flank U∥ lobe also sits in the cells where the four looks
disagree on the east-lobe amplitude (probe-disk closure ≈3 σ), where a
cross-sensor amplitude difference can only be absorbed by the weakly
constrained component. U∥ is therefore reported with these caveats and not
interpreted.

## GSLC reprocessing (mh26.gslc)

The GUNW chain works at the product's 80-m unwrapped posting and inherits
the product's unwrapping. `mh26.gslc` rebuilds the two co-event
interferograms from the 10-m L2 GSLC products: ranged-read crops of the HH
layer (chunk-aligned row bands through earthaccess/fsspec, with the GSLC
`mask` and `inputDataExceptionMask` rasters kept — the provisional mask
encodes 0 = no subswath, 1 = subswath 1, 255 = fill; invalid epoch samples
zero the interferogram pixel), per-pixel interferogram (GSLC epochs share
one UTM grid and are DEM-flattened; the grid is 10 m in x and 5 m in y, so
6 × 3 looks give the 30-m posting), Goldstein filter (α = 0.6, 64-px
tiles), boxcar coherence, snaphu, the GUNW ionosphere screen of the same
pair interpolated and subtracted, and far-field plane referencing. The GSLC
products carry the same provisional calibration maturity as the GUNW
products.

Unwrapping configuration matters at the first-fringe level. The primary run
uses one tile, DEFO costs, and MCF initialization. For the ascending
model-input pair, its comparison with GUNW gives slope +1.01, r = 0.993,
and 1.5 cm RMS difference on 1.28-km cells. The remaining difference is
concentrated near the axis, where the 30-m and 80-m processings can keep or
unwrap different cells. `--compare` writes region statistics and acceptance
limits for each pair. The limits are fatal only for the two inversion pairs.

`--unwrap-sensitivity` keeps the same interferogram, corrections, one-tile
layout, and MCF initialization, but changes the cost mode from DEFO to
SMOOTH. For the ascending model-input pair, the two cases differ by 4.60 cm
RMS near the axis, and 6.15% of near-axis pixels change by more than half a
fringe. Off axis within 25 km, the RMS difference is 0.14 cm and 0.014% of
pixels change by more than half a fringe. The command records input and
primary-output hashes in
`results/gslc_unwrap_sensitivity_172_008_A.json`. This is a processing
sensitivity, not an independent unwrap validation.

For the third pair (track 100 ascending, 15→27 July) the comparison itself
is the finding: the two unwrapped fields agree at the millimetre level over
three quarters of the common cells, and differ by one L-band fringe
(≈11.9 cm) over a sharp-edged island covering the east-flank lobe
(difference quartiles −0.1/0.0/+0.1 cm, 95th percentile +11.7 cm). The
reprocessed field is smooth and continuous across that lobe, while the GUNW
field contains an internal step. These observations alone do not decide
which fringe count is correct. No reported model quantity uses the GUNW
amplitude of this pair; it enters the analysis as a qualitative timing
bracket. The probe-history row for this window uses the GUNW value.

## Phase-gradient ridges (mh26.faults) and the ridge-anchored fault test

Wrapped-phase gradients of the 30-m reprocessed interferograms resolve two
coherent gradient ridges flanking the axis, detected independently in both
NISAR geometries (median cross-geometry spread 300 m) over ~25 km of
strike: the east ridge ≈1.6 km and the west ridge ≈1.3 km from the dike
trace. The picker selects at most one ridge per side and 500-m bin, and a
smooth steep gradient band can pass its quality threshold, so the ridges
are candidate scarp lines, not mapped faults. A LOS-step diagnostic across
each ridge (unwrapped fields; 0.1–0.4 km bands each side, local trend
removed) tests the discontinuity interpretation: the east ridge shows an
axis-side-down step in both geometries (median −4.8 cm ascending, −1.7 cm
descending; 79%/63% of bins negative) — the signature of a downthrown
scarp — while the west ridge does not (+2.9 cm ascending, −2.9 cm
descending), so the west line remains a candidate only. The strike-frame
sign convention is pinned to `mh26.distributed._strike_vectors`
(ENE-positive); the extraction is label-symmetric, so an inconsistent sign
convention swaps the east/west labels without changing the ridge geometry.
The result also stores sensitivity to the coherence and peak-quality
thresholds, plus a circular-shift null for cross-geometry alignment along
strike. These diagnostics test the picker and do not turn either ridge into
an independently mapped fault.

Anchoring the inversion's fault planes at the ridges is a geometry test,
with one construction rule: each inward-dipping column ends where it meets
the vertical dike plane (depth = |offset| × tan dip), so no fault patch
continues through the dike onto the far side of the rift; ridge offsets
under ~1.9 km truncate their columns. The default mapped fit includes the
33 of 42 nominal east patches that remain and excludes the west trace because
it fails the cross-geometry LOS-step test. This east-trace fit has 2.15 cm
joint RMS, compared with 2.09 cm for the two-plane idealization; its 5-km
checkerboard holdout improves on both splits. The regularized east slip
maximum is 0.67 m, and the fault/dike potency ratio is 0.102 (idealized:
0.196 — the partition between opening and slip shifts with the assumed fault
area). In the saved 55–75° dip sweep, the GUNW chain spans 0.019 cm joint RMS
and the GSLC chain spans 0.058 cm. Their minima do not agree, and the GSLC
minimum is at the grid boundary, so the data do not give a stable dip
estimate. Two scope points apply to the validation. The checkerboard holdout
is conditional — the ridge geometry
was derived from the full fields before the split, and 5-km cells are not
independent under the observed noise correlation lengths. NISAR→Sentinel-1 is
the stronger transfer direction because Sentinel-1 did not set the geometry.
The fault term also improves S1→NISAR prediction, but that reverse direction
is conditional because the NISAR test fields affected geometry selection.
Both parameterizations reduce in-sample and held-out residuals relative to
opening only. These conditional tests do not establish the trace geometry;
the 2005–2010 geometry files remain the discriminating test. The idealized
fit is kept as the reference (`distributed_fit.*`) and the ridge-anchored fit as
`distributed_fit_mapped.*`. That default mapped fit includes only trace sides
that pass the saved LOS-step rule. `--traces FILE` runs the same inversion on
all traces in an externally supplied file with the `fault_traces.json` schema.
It writes content-hash-suffixed outputs and records the source path and SHA-256
hash, so it cannot replace the normal mapped outputs. Stress can use that exact
result with `--fit`.

`--geometry FILE` accepts either the complete `geometry` object stored in a
distributed result or a JSON file that contains that contract directly. It
validates the patch arrays, dislocation rows, surface shapes, masks, axes,
areas, and patch order. It then uses those patches without rebuilding them
from current code constants. The result uses a content-hash suffix and records
the source path, source-file hash, and geometry hash. This closes the generic
input and provenance gap. The 2005 plane arrived exactly this way
(`data/complex_geom2_revised_test.dat` converted to
`results/geometry_wright2006.json`), and its refit is complete on both input
chains. The 2006–2010 geometry and opening files are still missing, and an
unknown legacy format can still require a converter when they arrive.

## Dike-induced stress on the fault surfaces (mh26.stress)

The static stress change of the saved distributed opening, with fault slip
excluded from the source, is resolved for inward-dipping normal slip at
receivers built from the geometry stored in the same fit. The fault patch
centres and a 250-m depth profile follow each saved trace column. Receivers
within 400 m horizontally of the dike plane are excluded as near-singular.
The 65° result is stored as spatial normal-stress and Coulomb-stress grids on
the supported receiver plane and as an along-strike depth summary. The
shallowest receiver row is mainly unclamped, while most deeper rows are
clamped. The distributed opening is parameterized from the surface to 10 km
and has no subsurface hard top. The report therefore makes no conclusion
relative to the uniform rectangle's fitted top. No GSLC stress case or
input-chain-independence claim is reported.

The strain-derivative convention of the dislocation code is checked
against central finite differences by the module self-test
(`python -m mh26.stress --selftest`), which also runs at the start of
every stress computation.

## Correlated-noise Monte Carlo (mh26.uncertainty)

Per-dataset nugget + exponential stationary kernels are fitted to the
far-field semivariograms of the referenced block observations (range searched
to 30 km; the fit parameters, binned semivariograms, and any search-limit hits
are stored in the result). A covariance fit to referenced values does not
impose the plane-removal constraints. The best-model prediction and every
stationary covariance draw are therefore passed through the block-space
far-field datum and ramp projection before each synthetic field is inverted
with fixed saved geometry, bounds, and smoothing. The noiseless refit differs
from the best fit through the smoothing penalty acting again on an
already-smooth input (the far-field projection is part of the estimator
itself, so it adds no refit bias) — so raw
refit percentiles are reported only as estimator dispersion, and the 95%
intervals are pivotal (basic bootstrap): refit errors about the noiseless
refit, re-centred on the best-fit values. The intervals are conditional
data-noise dispersion under this simulator, for the idealized-geometry
GUNW fit at the selected smoothing; they are not a full
physical-parameter confidence interval. Geometry, smoothing, and
processing choices are bracketed separately by the sensitivity envelopes,
and positivity truncates the noise response near zero slip.

## Decomposition frame-coverage rule

Lattice cells that straddle a NISAR swath/subswath boundary or a
product-mask gap average partly filled, edge-degraded samples: a rim of
ascending-frame edge cells (≈−14 cm) combined with a Sentinel-1
tropospheric patch produced a spurious localized feature near UTM
(670, 1338) km in early decomposition runs. The decomposition therefore
requires every NISAR cell to be at least 90% covered by the product
valid-sample mask. The criterion uses the mask alone — it is
value-independent, applies to every NISAR geometry everywhere (coherence
plays no role, so sparse-but-in-frame axial cells survive), and involves
no post-hoc spatial selection. LiCSAR frames have no swath boundary inside
the lattice, so no equivalent rule exists for Sentinel-1; the residual
Sentinel-1 tropospheric patch near (670, 1338) km remains visible in
low-redundancy cells and is disclosed in LIMITATIONS.

## Decomposition probe-disk closure and the east-flank systematic

The radial closure summaries hide local systematics, so closure is also
reported inside the 2-km probe disks that define the reported divergence
and subsidence values. The
EAST-flank redundant cells close at ≈3.1 σ (all above 1 σ) while WEST
closes at ≈0.6 σ and the axis disk has no redundant cells: the four looks
disagree systematically on the east-lobe amplitude, consistent with the
Sentinel-1-vs-NISAR amplitude ratios (0.88–0.92) of the independent S1
cross-checks. Rerunning the decomposition with the 30-m reprocessed NISAR
fields (`--gslc`) moves the divergence by 2–3 cm and leaves the EAST
closure unchanged (≈3.2 σ), which locates the systematic between the
sensors, not between the NISAR processing chains. The divergence values
are therefore cross-sensor compromises at the east flank; the spread
across input chain and normalization choices (+112 to +121 cm) is the
honest width, and the source inversion — which fits each look with its own
noise weight instead of forcing per-cell agreement — is the primary
quantitative product.

## Parameter sensitivity

The formal covariance in `results/okada_fit.json` assumes independent white
noise and is only a numerical diagnostic. Nine processing-choice runs vary
block size, far-field radius, LOS interpolation height, geometry weighting,
and fit radius. Their range describes processing sensitivity, not statistical
uncertainty.

The valid-sample and optimization changes alter the fitted geometry:

| parameter | earlier processing | current processing |
|---|---:|---:|
| length | 18.1 km | 17.31 km |
| depth | 0.8–7.0 km | 1.30–5.94 km |
| uniform opening | 1.58 m | 2.32 m |
| geometric volume | 0.18 km³ | 0.187 km³ |
| NISAR RMS | 3.8 cm | 2.91 cm |

All files in `results/` use the current processing.

## Summary figures (mh26.summary_figures)

The eight `results/summary_*` figures are display products built from the
saved results; they change no analysis quantity. The two Wright-plane
figures draw the saved 30-m-input fixed-plane fit; because the segmented
plane has no single straight frame, their map trace and profile frame come
from the saved patch arrays (opening-weighted centre and local segment
strike), not from the chord fields. Display choices that
matter for reading them:

- Model fields are evaluated on a 200-m display grid. The inversion support
  stays at the ~1.3-km blocks, and the printed RMS values are the saved
  block-level misfits.
- Observations are decay-normalized to the common 18 July reference with
  the per-dataset factors stored in the fit, the same normalization the
  inversion uses.
- Rewrapped maps use one colour cycle per 10 cm of LOS displacement toward
  the satellite. Residual maps are clipped at ±10 cm and keep coherent
  near-axis residuals visible.
- The decomposition panels show the unsmoothed 1.28-km solve. In the
  E–N–U figure the weakly observed North panel dots cells with formal σ
  above 6 cm as a display mark; the values under the dots are unchanged.
  The dike-frame figure is the same solve rotated, with the four
  decay-normalized LOS inputs as its top row; its U∥ panel masks cells
  with formal σ above 5 cm, and its frame arrows show +U⊥ (ENE) and +U∥
  (NNW, near-north).
- The signed cross-rift phase change is the wrapped one-pixel forward
  difference projected ENE (rift-perpendicular). It aliases where the
  cross-rift phase rate exceeds π per 30-m pixel (λ/4 of LOS per pixel);
  the white map veil and the gray profile band mark that zone, estimated
  from the unwrapped LOS on the same stencil. The forward-difference
  stencil shifts the signed map about 15 m WNW of true gradient positions.
- The stress panels resolve the saved distributed dike opening (fault slip
  excluded) on the trace-supported east-fault receiver orientation; the
  receiver parameters are printed on the figure. The GSLC/GUNW comparison
  uses identical receivers from the shared fit geometry.
- The terrain background is a public Copernicus GLO-30 tile
  (`data/README.md`), hash-verified at run time and used for display only.
- Cross-sections in the source and gradient figures sit at the along-strike
  position of the peak-slip east-fault column, derived from the saved fit.

## Archive verification

The core, pre-event, and complete data sets have separate SHA-256 manifests.
The core archive can therefore be verified without the optional pre-event NPZ
files. Release archives use the source commit time, record the source commit
hash, and pass a ZIP CRC check before replacement of an existing archive.

## Additional analysis

The following work is not part of the current results:

- estimation with a spatial covariance model;
- comparison with alternative source geometries;
- source modelling of the post-event signal itself;
- integration with long-term deformation and local seismic data.

The validator checks data schema, coordinates, masks, LOS vectors, raster
alignment, numerical ranges, source manifests, and SHA-256 hashes. The tests
cover the validity mask, coherence encoding, block averaging, and LOS-vector
interpolation.
