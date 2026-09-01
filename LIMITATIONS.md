# Limitations

## Observations

- The inputs are NISAR L2 GUNW Provisional V1 products. The analysis uses the
  supplied valid-sample mask because unreliable subswath-edge samples can
  otherwise survive coherence and unwrapping screens. NISAR also documents
  possible ionosphere, RFI, coregistration, and edge artifacts in early
  products.
- The supplied ionosphere phase is subtracted. Wet/hydrostatic troposphere and
  solid-Earth-tide lookup corrections are not applied; most prepared crops do
  not carry those optional cubes. A per-pair far-field plane removes a datum
  and linear ramp but is not an atmospheric correction. That plane also
  carries away the deformation field's own far-field component, so every
  mapped field and the decomposition are relative to this datum; the
  inversions pass the model through the same projection
  (reference-operator consistency), so the fitted sources are consistent
  with the operator.
- The reprocessed (GSLC) interferograms set invalid source pixels to zero
  before multilooking and store no per-cell valid fraction, so a cell with
  few valid samples can keep high coherence; a stored common-valid
  fraction is planned for the next reprocessing.
- Pixel and block samples are spatially correlated. Counts in reports are
  retained-sample counts, not independent degrees of freedom.
- The ascending/descending decomposition has only two look directions and
  fixes north displacement to zero. It is a quasi-East/quasi-Up diagnostic,
  not a three-component solution.
- The joint four-geometry decomposition constrains horizontal motion to the
  rift-perpendicular azimuth (taken from the fitted strike, so it is not
  independent of the source model) instead of assuming zero north motion.
  That is exact only if horizontal motion is opening-normal; along-strike
  motion maps into the closure residual, not into the solved components. The
  free three-component variant removes the constraint at the price of
  several-cm formal north sigmas. Its dike-frame rotation inherits that:
  U∥ (along strike, nearly north–south) has a median formal σ of 4.1 cm,
  its flank values change by several cm between the GUNW and 30-m input
  chains and with a 4° change of the frame strike, and its east-flank lobe
  sits where the looks disagree on amplitude, so U∥ is reported but not
  resolved; U⊥ is the same axis as the constrained h. The axial
  subsidence probe rests on
  two-look NISAR-only cells: Sentinel-1 does not survive at the axis. All
  joint products also depend on the decay normalization that maps the four
  windows to a common 18 July epoch; it assumes the post-event transient
  shares the co-event dike pattern AND that every deformation component —
  including any purely coseismic fault slip — carries the same transient
  fraction, which is untested; turning it off moves the flank divergence
  from +112 to +117 cm, and both variants are reported as assumptions, not
  as a neutral epoch conversion. The redundant cells of the EAST-flank
  probe disk close at ≈3 σ (a systematic cross-sensor amplitude
  difference, unchanged when the reprocessed 30-m NISAR fields replace the
  GUNW fields), so the divergence values are cross-sensor compromises
  there; see PROCESSING_NOTES. A residual Sentinel-1 tropospheric patch
  near UTM (670, 1338) km survives in low-redundancy cells and is visible
  in the decomposition maps; the flank and axial probes are more than 8 km
  from it.
- The reprocessed 30-m interferograms (mh26.gslc) carry the same provisional
  calibration maturity as the GUNW products, reuse the GUNW ionosphere screen
  (an 80-m product layer interpolated to 30 m), and depend on the documented
  snaphu configuration. In the saved ascending-pair sensitivity, changing
  DEFO costs to SMOOTH costs changes 6.1% of near-axis pixels by more than
  half a fringe, while off-axis agreement remains close. Agreement with the
  GUNW fields is 1.1–1.5 cm RMS on 1.28-km cells. The near-axis band has no
  independent registered unwrap check yet.

## Source model

- The inversion is one vertical, uniform-opening rectangle in a homogeneous,
  isotropic elastic half-space. It omits nonuniform opening, shallow normal
  faulting, topography, elastic layering/heterogeneity, and any coupled
  reservoir source. Near-axis residuals show that the rectangle is not a
  complete deformation model.
- Formal covariance in `okada_fit.json` assumes white residual noise and is
  only a numerical diagnostic. The processing-choice range is also not a
  confidence or credible interval. A statistical uncertainty analysis requires
  spatial covariance and comparison with alternative source geometries.
- Length × height × uniform opening is geometric volume in this idealized
  model. `μV` is a comparison scale; it is not an observed seismic moment and
  should not be used alone to characterize seismic efficiency.
- The distributed model fixes all geometry a priori (dike plane and strike
  from the uniform fit; two graben faults at ±2 km, dip 65°, 0–4 km) and
  inverts only non-negative dislocation magnitudes with Laplacian smoothing.
  The smoothing weight is a held-out-prediction choice, not a resolved
  parameter: the L-curve curvature criterion is boundary-unstable on this
  problem (see PROCESSING_NOTES), the 10-km selection cells are not fully
  independent under the observed noise correlation lengths, and the
  λ×1/3.2–λ×3.2 envelope in the sensitivities shows how strongly the slip
  maxima (and, less so, the volume) depend on it. Positivity plus smoothing
  bias the solution, edge patches are poorly resolved, and no formal
  resolution analysis is computed; the sensitivities bracket the main
  choices only. The summed fault potency is
  repeatable across data subsets, but the east/west partition is not (NISAR
  alone puts the maximum on the east fault, Sentinel-1 alone on the west), so
  per-fault slip values are regularized magnitudes, not resolved slip. The
  fault/dike ratio is a potency ratio (area × dislocation), not a
  moment-tensor ratio. Dike volumes from different model classes are not
  directly comparable: the distributed model moves near-axis signal from
  opening into fault slip, so its 0.17 km³ sits below the uniform-rectangle
  0.187 km³ without contradiction. The 2005 plane itself is now refit
  directly (the fixed 2005-plane entry below); whether the 2026 plane
  coincides with the 2006–2010 dike geometries remains untested until those
  geometry files are available. The complete JSON geometry contract and hash-safe input command
  are ready, but an unknown legacy file format can still require a converter.
- The phase-gradient ridges are candidate scarp lines: the picker selects
  at most one ridge per side and bin, and only the east ridge passes the
  LOS-step discontinuity test in both geometries (the west ridge's step
  changes sign between geometries). The saved threshold sensitivity and
  circular-shift alignment null test the picker, but no comparison against an
  independently mapped scarp database exists yet. The step rule itself
  checks only the sign and the down-side fraction of steps measured at the
  selected maximum-gradient positions after a linear outer detrend; it is
  not calibrated against a smooth-source null, so it is a consistency
  check, not a fault detection. The default trace-column
  fit therefore includes only the east ridge. Its surface keeps one fixed
  strike with per-column offsets (plan-view steps are parameterisation, not
  resolved geometry) and ends at the dike plane. Its RMS is 2.15 cm, compared
  with 2.09 cm for the ±2-km two-plane idealization. The 55–75° dip sweeps do
  not give a stable common minimum between the GUNW and GSLC chains, and the
  slip values remain regularized magnitudes. The checkerboard holdout is
  conditional on this data-derived geometry (the ridges were picked from the
  full fields);
  The held-out-sensor transfer gain in `distributed_fit.md` belongs to the
  idealized two-fault geometry; the mapped east-trace term is constrained
  mainly by NISAR (a Sentinel-1-only fit puts an order of magnitude less
  slip on it). S1→NISAR is conditional because NISAR affected geometry
  selection.
  The slip-depth distribution is also shaped by regularization. The stress
  analysis maps the supported receiver plane and reports a shallow-to-deep
  sign change, but the distributed source has no subsurface hard top. The
  stress-change and Monte Carlo results inherit the geometry
  and noise-model assumptions in their reports. The Monte Carlo intervals
  cover data noise for the idealized GUNW fit only. Two of the four fitted
  noise-covariance ranges reached the former 30-km search limit; the
  search now extends to 60 km.
- The fixed 2005-plane fit (`mh26.wright2006`) is a fixed-historical-plane
  variant, not proof of co-planarity: the 2026 uniform-fit strike differs
  from the southernmost 2005 segment by about 4.5°, and its endpoints lie
  0.45 and 0.92 km off that segment's line. The plane is the Wright et
  al. (2006) source geometry, not the longer 1 km × 1 km computational grid
  of Hamling et al. (2009); same-plane comparisons must say which. Its 10-km
  SE extension is our addition: about 12% of the fitted volume lies SE of
  the historic plane's termination, while the opening tapers to zero well
  before the extended edge (the panel step reports both), so the statement
  "the 2026 dike extends past the 2005 plane's SE end" is a model result on
  an invented extension, not an observation independent of it. The reused
  smoothing weight is a numeric λ on a
  grid with different patch widths and twice the parameters of the primary
  fit, so regularization is matched only nominally; the reported λ-envelope
  brackets the effect.

## Timing and continuation

- The 15–18 July onset bracket combines a NISAR control ending 15 July with a
  Sentinel-1 interferogram ending 18 July. It is an interval bound, not an
  event time or duration estimate.
- The post-event continuation is now observed in three independent looks:
  the NISAR descending pair and six Sentinel-1 pairs in two geometries, with
  nine pre-event Sentinel-1 pairs as an empirical null. This supports the
  statement that dike-patterned deformation continued after 18–19 July and
  decayed within about two weeks. It does not upgrade the coefficients into
  magma volumes: every projection depends on a per-pair intercept, the
  processing envelopes and cluster-wild intervals do not sample atmospheric
  covariance, and the decay fit is a two-parameter descriptive summary of
  window-integrated values with correlated noise.
- The post-event East/Up decomposition combines ascending and descending
  windows that differ by five days at each end. For a decaying signal the
  split blends the two windows, so the divergence and subsidence values are
  indicative, not window-consistent displacements.
- Sixteen pre-event windows constrain only the listed pairs and three 1-km
  probes. “No precursor above 3.4 cm” is a probe/history statement, not proof
  that no spatially distributed or temporally unresolved precursor occurred.

## External data and reuse

- Long-term InSAR time-series results and local seismic products are outside
  this repository. Statements based on unpublished external data require
  permission and attribution.
- LiCSAR files are Copernicus-derived COMET products subject to the portal's
  use conditions, acknowledgement, and requested citations. They are not
  labelled CC-BY here.
- No public licence has been selected for this code or the derived figures.
  Reuse terms and archive metadata remain to be defined before any public
  redistribution.
