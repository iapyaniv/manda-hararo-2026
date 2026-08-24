# Correlated-noise Monte Carlo (distributed model)

- target: the idealized-geometry GUNW fit (distributed_fit.json); 1000 realizations at fixed geometry and λ = 17.78
- noise: per-dataset nugget + exponential stationary kernel fitted to the far-field semivariogram of the referenced observations; the best-model prediction and every covariance draw receive the block-space far-field datum/ramp projection
- intervals: pivotal (basic bootstrap) about the noiseless refit, re-centred on the best fit. The noiseless refit differs from the best fit through the combined estimator bias — the smoothing penalty acting again plus the referencing removing part of the model's own far-field signal — so raw refit percentiles are reported separately as estimator dispersion. The intervals are conditional data-noise dispersion under this simulator, not a full physical-parameter confidence interval

| quantity | best fit | noiseless refit | 95% interval |
|---|---|---|---|
| dike_volume_km3 | 0.145 | 0.132 | 0.132–0.157 |
| max_opening_m | 2.416 | 2.167 | 2.311–2.537 |
| fault_east_max_slip_m | 0.660 | 0.581 | 0.580–0.752 |
| fault_west_max_slip_m | 0.368 | 0.338 | 0.233–0.468 |
| fault_to_dike_potency_ratio | 0.235 | 0.265 | 0.143–0.304 |
| opening_centroid_depth_km | 3.589 | 3.437 | 3.378–3.774 |

| quantity | raw refit 2.5% | 50% | 97.5% |
|---|---|---|---|
| dike_volume_km3 | 0.120 | 0.132 | 0.144 |
| max_opening_m | 2.046 | 2.163 | 2.271 |
| fault_east_max_slip_m | 0.489 | 0.580 | 0.660 |
| fault_west_max_slip_m | 0.238 | 0.340 | 0.473 |
| fault_to_dike_potency_ratio | 0.195 | 0.268 | 0.356 |
| opening_centroid_depth_km | 3.252 | 3.445 | 3.649 |

The intervals describe data noise only, for fixed geometry and
regularization; positivity truncates the noise response near zero
slip, and the pivotal re-centring assumes the refit error about the
noiseless refit approximates the estimator error about the truth.
Geometry and processing choices are bracketed separately by the
sensitivity envelopes in the distributed-fit reports. Range search hit its 30 km limit for: NISAR desc 07-07→07-19, S1 desc 07-11→07-23.
