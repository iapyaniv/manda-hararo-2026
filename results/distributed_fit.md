# Distributed opening + graben-fault model (four geometries)

- fault geometry: idealized traces ±2 km from the dike, dip 65° inward
- data: NISAR asc 07-08→07-20, NISAR desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14053 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 chosen by minimum 10-km checkerboard held-out RMS on the grid 0.032–1000 (the discrete L-curve curvature maximum is not bracketed by any moderate grid on this problem and is reported for reference only)

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.49 cm | 2.96 cm | 2.70 cm | 1.70 cm | 2.32 cm |
| distributed opening | 2.37 cm | 2.83 cm | 2.84 cm | 1.32 cm | 1.88 cm |
| distributed opening + graben faults | 2.08 cm | 2.73 cm | 2.27 cm | 1.24 cm | 1.73 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.170 km³** (opening-only variant 0.186); max opening **2.34 m** at 4.1 km centroid depth
- graben faults: regularized maximum normal slip east 0.64 m, west 0.38 m; fault/dike potency ratio 0.196 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.51 cm → with faults 2.21 cm
- train parity 1: held-out RMS opening-only 2.63 cm → with faults 2.38 cm

## Sensitivities (volume km³ / max opening m / east–west max slip m / RMS cm)

- smoothing λ×1/3.2: 0.164 / 2.32 / 1.12–0.53 / 1.94
- smoothing λ×3.2: 0.187 / 1.85 / 0.34–0.22 / 2.33
- NISAR only: 0.184 / 2.30 / 0.70–0.36 / 2.43
- Sentinel-1 only: 0.169 / 1.98 / 0.10–0.48 / 1.41
- faults ±1.5 km: 0.174 / 2.42 / 0.45–0.46 / 2.17
- faults ±2.5 km: 0.166 / 2.23 / 0.70–0.34 / 1.99
- fault dip 55°: 0.168 / 2.29 / 0.64–0.53 / 2.04
- fault dip 75°: 0.175 / 2.38 / 0.61–0.30 / 2.15
- no decay normalization: 0.178 / 2.41 / 0.70–0.36 / 2.18

## Smoothing selection (10-km checkerboard held-out RMS)

| λ | fit RMS | held-out RMS (parity mean) |
|---|---|---|
| 0.0316 | 1.86 cm | 4.371 cm |
| 0.0562 | 1.86 cm | 4.386 cm |
| 0.1 | 1.86 cm | 4.402 cm |
| 0.178 | 1.87 cm | 4.308 cm |
| 0.316 | 1.87 cm | 4.073 cm |
| 0.562 | 1.87 cm | 3.808 cm |
| 1 | 1.87 cm | 3.580 cm |
| 1.78 | 1.88 cm | 3.256 cm |
| 3.16 | 1.90 cm | 2.958 cm |
| 5.62 | 1.94 cm | 2.699 cm |
| 10 | 2.00 cm | 2.537 cm |
| 17.8 | 2.08 cm | 2.462 cm ← |
| 31.6 | 2.19 cm | 2.477 cm |
| 56.2 | 2.33 cm | 2.556 cm |
| 100 | 2.50 cm | 2.675 cm |
| 178 | 2.70 cm | 2.981 cm |
| 316 | 3.21 cm | 3.943 cm |
| 562 | 4.50 cm | 5.453 cm |
| 1e+03 | 5.95 cm | 6.560 cm |

- curvature knee for reference: λ = 178
- cross-sensor transfer at the chosen λ (offsets refit): NISAR->S1: opening-only 1.81 → with faults 1.70 cm; S1->NISAR: opening-only 3.37 → with faults 3.15 cm
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
