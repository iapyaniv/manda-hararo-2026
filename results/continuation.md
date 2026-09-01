# Post-event continuation test

- NISAR descending pair: **2026-07-19→2026-08-12** (24 days; starts at the co-event
  descending pair's secondary acquisition)
- model: results/okada_fit.json (SHA-256 `7f96d2fa7a2076264ca41ab2d7c4c44bc321bde6f59876574eacaf4a56f26f9e`)
- dike-pattern projection coefficient: **+0.107**; processing-choice
  envelope **[+0.100, +0.121]**
- conditional 5-km cluster-wild interval: **[+0.085,
  +0.128]** (87 clusters); this is not a calibrated
  confidence interval for spatially correlated InSAR noise
- correlation with the co-event pattern: +0.629; residual RMS: 1.7 cm
- fitted subset offset: +3.3 cm; forcing it to zero would give a
  coefficient of +0.183
- probe means: EAST +1.3 cm, WEST +12.4 cm,
  SOUTH +7.6 cm

The positive projection coefficient is about one-tenth of the co-event value.
This result uses one provisional descending geometry and does not establish
renewed intrusion.

The test detects only deformation shaped like the fitted dike. The broad residual
field and fitted datum offset show why additional geometries and atmospheric
covariance treatment are required before converting this coefficient to added volume.
