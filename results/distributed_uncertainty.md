# Correlated-noise Monte Carlo (distributed model)

- target: the idealized-geometry GUNW fit (distributed_fit.json); 1000 realizations at fixed geometry and λ = 17.78
- noise: per-dataset nugget + exponential stationary kernel fitted to the far-field semivariogram of the referenced observations; the best-model prediction and every covariance draw receive the block-space far-field datum/ramp projection
- intervals: pivotal (basic bootstrap) about the noiseless refit, re-centred on the best fit. The noiseless refit differs from the best fit through the combined estimator bias — the smoothing penalty acting again plus the referencing removing part of the model's own far-field signal — so raw refit percentiles are reported separately as estimator dispersion. The intervals are conditional data-noise dispersion under this simulator, not a full physical-parameter confidence interval

| quantity | best fit | noiseless refit | 95% interval |
|---|---|---|---|
| dike_volume_km3 | 0.170 | 0.176 | 0.154–0.186 |
| max_opening_m | 2.340 | 2.145 | 2.231–2.449 |
| fault_east_max_slip_m | 0.640 | 0.482 | 0.569–0.724 |
| fault_west_max_slip_m | 0.383 | 0.348 | 0.254–0.485 |
| fault_to_dike_potency_ratio | 0.196 | 0.182 | 0.133–0.243 |
| opening_centroid_depth_km | 4.144 | 4.301 | 3.867–4.411 |

| quantity | raw refit 2.5% | 50% | 97.5% |
|---|---|---|---|
| dike_volume_km3 | 0.160 | 0.177 | 0.193 |
| max_opening_m | 2.036 | 2.146 | 2.254 |
| fault_east_max_slip_m | 0.398 | 0.482 | 0.553 |
| fault_west_max_slip_m | 0.246 | 0.351 | 0.477 |
| fault_to_dike_potency_ratio | 0.135 | 0.183 | 0.245 |
| opening_centroid_depth_km | 4.034 | 4.318 | 4.579 |

The intervals describe data noise only, for fixed geometry and
regularization; positivity truncates the noise response near zero
slip, and the pivotal re-centring assumes the refit error about the
noiseless refit approximates the estimator error about the truth.
Geometry and processing choices are bracketed separately by the
sensitivity envelopes in the distributed-fit reports.
