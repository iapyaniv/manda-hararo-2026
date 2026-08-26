# Distributed opening + graben-fault model (four geometries; complete external geometry)

- fault geometry: complete patch geometry from results/geometry_wright2006.json (file SHA-256 86ffda58012256c63fe8df4c762604927a8dbf62ba3997eedad55644268f13c0; geometry SHA-256 0e587ed44c3ee32399410621d706e814f877fb6b18b2eabf129fbd7fae53ba7a)
- data: NISAR asc 07-08→07-20, NISAR desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14053 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.71 cm | 3.01 cm | 2.89 cm | 2.25 cm | 2.56 cm |
| distributed opening | 2.64 cm | 3.38 cm | 2.83 cm | 1.83 cm | 2.20 cm |
| distributed opening + graben faults | 2.46 cm | 3.19 cm | 2.47 cm | 1.86 cm | 2.17 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.156 km³** (opening-only variant 0.165); max opening **2.44 m** at 3.6 km centroid depth
- retained fault term: regularized maximum normal slip east 0.71 m; fault/dike potency ratio 0.120 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.77 cm → with faults 2.64 cm
- train parity 1: held-out RMS opening-only 2.92 cm → with faults 2.69 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.149 / 2.59 / 1.04 / 2.36
- smoothing λ×3.2: 0.171 / 1.86 / 0.42 / 2.71
- NISAR only: 0.180 / 2.30 / 0.81 / 2.73
- Sentinel-1 only: 0.150 / 2.08 / 0.11 / 1.89
- no decay normalization: 0.161 / 2.50 / 0.77 / 2.57

## Scope

- The dike plane, strike, and fault geometry are fixed a priori; the
  inversion solves only dislocation magnitudes (convex, positive).
- The patch rows, surface masks, strike, dip, and depth axes come
  from the recorded geometry file. The inversion does not rebuild
  them from current source-code constants.
- Source-specific geometry sensitivity cases need separate geometry
  files. No hard-coded offsets or dips are applied to this run.
- The retained east-fault term improves both sensor-subset fits, but
  its magnitude changes between them; treat its slip values as
  regularized magnitudes, not resolved slip. The west trace is absent
  because it fails the cross-geometry LOS-step support rule.
- Formal resolution/uncertainty of the slip distribution is not
  computed here; the sensitivity list brackets the main choices.
