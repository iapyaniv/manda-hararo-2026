# Distributed opening + graben-fault model (four geometries)

- fault geometry: idealized traces ±2 km from the dike, dip 65° inward
- data: NISAR asc 07-08→07-20, NISAR desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14053 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 chosen by minimum 10-km checkerboard held-out RMS on the grid 0.032–1000 (the discrete L-curve curvature maximum is not bracketed by any moderate grid on this problem and is reported for reference only)

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.71 cm | 3.01 cm | 2.89 cm | 2.25 cm | 2.56 cm |
| distributed opening | 2.59 cm | 2.96 cm | 3.02 cm | 1.90 cm | 2.08 cm |
| distributed opening + graben faults | 2.30 cm | 2.76 cm | 2.45 cm | 1.86 cm | 1.96 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.145 km³** (opening-only variant 0.162); max opening **2.42 m** at 3.6 km centroid depth
- graben faults: regularized maximum normal slip east 0.66 m, west 0.37 m; fault/dike potency ratio 0.235 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.73 cm → with faults 2.46 cm
- train parity 1: held-out RMS opening-only 2.84 cm → with faults 2.56 cm

## Sensitivities (volume km³ / max opening m / east–west max slip m / RMS cm)

- smoothing λ×1/3.2: 0.137 / 2.69 / 0.94–0.52 / 2.19
- smoothing λ×3.2: 0.161 / 1.81 / 0.39–0.23 / 2.57
- NISAR only: 0.166 / 2.31 / 0.77–0.34 / 2.55
- Sentinel-1 only: 0.137 / 2.04 / 0.10–0.51 / 1.82
- faults ±1.5 km: 0.149 / 2.48 / 0.55–0.44 / 2.38
- faults ±2.5 km: 0.141 / 2.31 / 0.72–0.33 / 2.25
- fault dip 55°: 0.142 / 2.38 / 0.66–0.51 / 2.28
- fault dip 75°: 0.149 / 2.45 / 0.64–0.31 / 2.35
- no decay normalization: 0.151 / 2.49 / 0.72–0.36 / 2.41

## Smoothing selection (10-km checkerboard held-out RMS)

| λ | fit RMS | held-out RMS (parity mean) |
|---|---|---|
| 0.0316 | 2.12 cm | 4.179 cm |
| 0.0562 | 2.12 cm | 4.189 cm |
| 0.1 | 2.12 cm | 4.232 cm |
| 0.178 | 2.12 cm | 4.143 cm |
| 0.316 | 2.12 cm | 3.981 cm |
| 0.562 | 2.13 cm | 3.700 cm |
| 1 | 2.13 cm | 3.470 cm |
| 1.78 | 2.14 cm | 3.220 cm |
| 3.16 | 2.16 cm | 3.008 cm |
| 5.62 | 2.19 cm | 2.828 cm |
| 10 | 2.23 cm | 2.708 cm |
| 17.8 | 2.30 cm | 2.674 cm ← |
| 31.6 | 2.42 cm | 2.707 cm |
| 56.2 | 2.57 cm | 2.810 cm |
| 100 | 2.77 cm | 2.970 cm |
| 178 | 2.99 cm | 3.249 cm |
| 316 | 3.43 cm | 4.050 cm |
| 562 | 4.54 cm | 5.433 cm |
| 1e+03 | 5.92 cm | 6.530 cm |

- curvature knee for reference: λ = 178
- cross-sensor transfer at the chosen λ (offsets refit): NISAR->S1: opening-only 2.29 → with faults 2.18 cm; S1->NISAR: opening-only 3.57 → with faults 3.33 cm
- validation role: NISAR→S1 is the held-out-sensor direction. S1→NISAR is conditional because the NISAR fields affected geometry selection.

## Scope

- The dike plane, strike, and fault geometry are fixed a priori; the
  inversion solves only dislocation magnitudes (convex, positive).
- The graben-fault geometry is a two-plane idealization chosen from the
  residual trough; a direct test against the 2005-2010 dike and fault
  geometries requires those geometry files.
- The total fault-like contribution is repeatable across data subsets,
  but the east/west partition is not (see the sensor-subset rows);
  treat per-fault slip values as regularized magnitudes, not resolved
  slip. The west-fault maximum sits at the model-grid edge.
- Formal resolution/uncertainty of the slip distribution is not
  computed here; the sensitivity list brackets the main choices.
