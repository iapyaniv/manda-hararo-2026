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
  and linear ramp but is not an atmospheric correction.
- Pixel and block samples are spatially correlated. Counts in reports are
  retained-sample counts, not independent degrees of freedom.
- The ascending/descending decomposition has only two look directions and
  fixes north displacement to zero. It is a quasi-East/quasi-Up diagnostic,
  not a three-component solution.

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

## Timing and continuation

- The 15–18 July onset bracket combines a NISAR control ending 15 July with a
  Sentinel-1 interferogram ending 18 July. It is an interval bound, not an
  event time or duration estimate.
- The 19 July–12 August continuation result is one 24-day descending look.
  Its projection coefficient depends on an intercept because a broad residual
  field remains: +0.107 with the intercept, +0.180 if the intercept is forced
  to zero. The included processing envelope varies block size, reference
  radius, fit radius, and signal threshold but does not sample atmospheric
  covariance. Do not translate it into added magma volume without independent
  geometries.
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
  Define reuse terms and archive metadata before public redistribution.
