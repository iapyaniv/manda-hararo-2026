# Okada dike fit — Manda Hararo July 2026

- observations: 7717 (block-avg 1280 m), RMS misfit **3.1 cm**
- screening is entirely a priori: coherence, connected component, and the
  product valid-sample mask; no blocks were rejected using model residuals
- optimizer: 6/6 starts completed; best optimality 0.00429

| param | value | formal 1σ* | unit |
|---|---|---|---|
| east_c | 683225.4 | ±25.2 | m |
| north_c | 1350366.9 | ±39.0 | m |
| strike | 152.6 | ±0.2 | deg |
| length | 17150.2 | ±119.9 | m |
| top | 1494.4 | ±42.9 | m |
| bottom | 4988.7 | ±115.5 | m |
| opening | 2.81 | ±0.11 | m |
| offs_asc | -0.0159 | ±0.0006 | m |
| offs_desc | -0.0073 | ±0.0005 | m |

*The formal errors are a numerical white-noise diagnostic only. InSAR noise is
spatially correlated and the uniform-opening half-space is an approximation, so
these values are **not** parameter confidence intervals. Processing-choice
sensitivity is reported separately in `model_sensitivity.md`.

- dike area 17.2 km × 3.5 km, opening 2.81 m ⇒ **geometric volume ≈ 0.17 km³**
- geodetic moment scale μV ≈ 5.06e+18 N·m (Mw≈6.4 only as a μV comparison; this is not a measured seismic moment)
