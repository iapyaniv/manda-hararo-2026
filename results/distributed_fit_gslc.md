# Distributed opening + graben-fault model (four geometries; NISAR from GSLC reprocessing)

- fault geometry: idealized traces ±2 km from the dike, dip 65° inward
- data: NISAR/GSLC asc 07-08→07-20, NISAR/GSLC desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14388 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 3.07 cm | 3.48 cm | 3.52 cm | 2.25 cm | 2.56 cm |
| distributed opening | 2.98 cm | 3.49 cm | 3.61 cm | 1.91 cm | 2.09 cm |
| distributed opening + graben faults | 2.74 cm | 3.52 cm | 3.06 cm | 1.85 cm | 1.99 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.147 km³** (opening-only variant 0.161); max opening **2.51 m** at 3.6 km centroid depth
- graben faults: regularized maximum normal slip east 0.60 m, west 0.37 m; fault/dike potency ratio 0.197 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 3.13 cm → with faults 2.87 cm
- train parity 1: held-out RMS opening-only 3.41 cm → with faults 3.19 cm

## Sensitivities (volume km³ / max opening m / east–west max slip m / RMS cm)

- smoothing λ×1/3.2: 0.138 / 2.75 / 1.08–0.49 / 2.63
- smoothing λ×3.2: 0.164 / 1.85 / 0.36–0.23 / 3.01
- NISAR only: 0.168 / 2.41 / 0.66–0.31 / 3.21
- Sentinel-1 only: 0.137 / 2.04 / 0.10–0.51 / 1.82
- faults ±1.5 km: 0.151 / 2.60 / 0.43–0.43 / 2.83
- faults ±2.5 km: 0.142 / 2.44 / 0.67–0.33 / 2.66
- fault dip 55°: 0.145 / 2.51 / 0.50–0.51 / 2.73
- fault dip 75°: 0.150 / 2.53 / 0.66–0.32 / 2.74
- no decay normalization: 0.153 / 2.58 / 0.66–0.34 / 2.86

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
