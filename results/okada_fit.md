# Okada dike fit — Manda Hararo July 2026

- observations: 7717 (block-avg 1280 m), RMS misfit **2.9 cm**
- screening is entirely a priori: coherence, connected component, and the
  product valid-sample mask; no blocks were rejected using model residuals
- the model LOS passes through the same far-field plane projection as the
  referenced observations (reference-operator consistency)
- optimizer: 6/6 starts completed; best optimality 0.00553

| param | value | formal 1σ* | unit |
|---|---|---|---|
| east_c | 683169.9 | ±23.9 | m |
| north_c | 1350402.3 | ±37.1 | m |
| strike | 153.0 | ±0.2 | deg |
| length | 17314.6 | ±114.4 | m |
| top | 1300.9 | ±31.4 | m |
| bottom | 5941.9 | ±113.0 | m |
| opening | 2.32 | ±0.05 | m |
| offs_asc | -0.0029 | ±0.0005 | m |
| offs_desc | -0.0079 | ±0.0004 | m |

*The formal errors are a numerical white-noise diagnostic only. InSAR noise is
spatially correlated and the uniform-opening half-space is an approximation, so
these values are **not** parameter confidence intervals. Processing-choice
sensitivity is reported separately in `model_sensitivity.md`.

- dike area 17.3 km × 4.6 km, opening 2.32 m ⇒ **geometric volume ≈ 0.19 km³**
- geodetic moment scale μV ≈ 5.60e+18 N·m (Mw≈6.4 only as a μV comparison; this is not a measured seismic moment)
