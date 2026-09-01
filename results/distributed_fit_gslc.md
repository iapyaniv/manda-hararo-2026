# Distributed opening + graben-fault model (four geometries; NISAR from GSLC reprocessing)

- fault geometry: idealized traces ±2 km from the dike, dip 65° inward
- data: NISAR/GSLC asc 07-08→07-20, NISAR/GSLC desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14388 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.87 cm | 3.40 cm | 3.37 cm | 1.70 cm | 2.32 cm |
| distributed opening | 2.76 cm | 3.30 cm | 3.44 cm | 1.36 cm | 1.90 cm |
| distributed opening + graben faults | 2.51 cm | 3.34 cm | 2.88 cm | 1.25 cm | 1.77 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.171 km³** (opening-only variant 0.185); max opening **2.40 m** at 4.1 km centroid depth
- graben faults: regularized maximum normal slip east 0.57 m, west 0.37 m; fault/dike potency ratio 0.164 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.92 cm → with faults 2.64 cm
- train parity 1: held-out RMS opening-only 3.15 cm → with faults 2.91 cm

## Sensitivities (volume km³ / max opening m / east–west max slip m / RMS cm)

- smoothing λ×1/3.2: 0.166 / 2.49 / 1.04–0.47 / 2.40
- smoothing λ×3.2: 0.189 / 1.88 / 0.34–0.23 / 2.75
- NISAR only: 0.183 / 2.36 / 0.62–0.33 / 3.03
- Sentinel-1 only: 0.169 / 1.98 / 0.10–0.48 / 1.41
- faults ±1.5 km: 0.175 / 2.48 / 0.39–0.45 / 2.60
- faults ±2.5 km: 0.166 / 2.33 / 0.64–0.33 / 2.42
- fault dip 55°: 0.169 / 2.39 / 0.46–0.53 / 2.49
- fault dip 75°: 0.175 / 2.42 / 0.62–0.31 / 2.53
- no decay normalization: 0.178 / 2.47 / 0.63–0.35 / 2.62

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
