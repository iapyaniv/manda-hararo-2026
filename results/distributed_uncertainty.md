# Correlated-noise Monte Carlo (distributed model)

- target: the idealized-geometry GUNW fit (distributed_fit.json); 1000 realizations at fixed geometry and λ = 17.78
- noise: per-dataset nugget + exponential stationary kernel fitted to the far-field semivariogram of the referenced observations; the best-model prediction and every covariance draw receive the block-space far-field datum/ramp projection
- intervals: pivotal (basic bootstrap) about the noiseless refit, re-centred on the best fit. The noiseless refit differs from the best fit through the estimator bias of the smoothing penalty acting again (the reference projection is applied inside the estimator itself, so it adds no bias here), and raw refit percentiles are reported separately as estimator dispersion. The intervals are conditional data-noise dispersion under this simulator, not a full physical-parameter confidence interval

| quantity | best fit | noiseless refit | 95% interval |
|---|---|---|---|
| dike_volume_km3 | 0.170 | 0.176 | 0.153–0.186 |
| max_opening_m | 2.337 | 2.143 | 2.228–2.446 |
| fault_east_max_slip_m | 0.638 | 0.481 | 0.567–0.722 |
| fault_west_max_slip_m | 0.384 | 0.349 | 0.255–0.486 |
| fault_to_dike_potency_ratio | 0.196 | 0.182 | 0.133–0.243 |
| opening_centroid_depth_km | 4.143 | 4.300 | 3.866–4.411 |

| quantity | raw refit 2.5% | 50% | 97.5% |
|---|---|---|---|
| dike_volume_km3 | 0.160 | 0.177 | 0.193 |
| max_opening_m | 2.034 | 2.144 | 2.252 |
| fault_east_max_slip_m | 0.397 | 0.481 | 0.552 |
| fault_west_max_slip_m | 0.247 | 0.351 | 0.478 |
| fault_to_dike_potency_ratio | 0.135 | 0.183 | 0.245 |
| opening_centroid_depth_km | 4.033 | 4.317 | 4.578 |

The intervals describe data noise only, for fixed geometry and
regularization; positivity truncates the noise response near zero
slip, and the pivotal re-centring assumes the refit error about the
noiseless refit approximates the estimator error about the truth.
Geometry and processing choices are bracketed separately by the
sensitivity envelopes in the distributed-fit reports.
