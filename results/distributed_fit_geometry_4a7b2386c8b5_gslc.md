# Distributed opening + graben-fault model (four geometries; NISAR from GSLC reprocessing; complete external geometry)

- fault geometry: complete patch geometry from results/geometry_wright2006.json (file SHA-256 4a7b2386c8b5eaf9d94b477492de7f5d3d5b312b14b77654d41a086873b86550; geometry SHA-256 23e42105d82ea9cdcadb770cccdbf9d9585b670e132b6c18c6e99d74daf8a01a)
- data: NISAR/GSLC asc 07-08→07-20, NISAR/GSLC desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14388 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.87 cm | 3.40 cm | 3.37 cm | 1.70 cm | 2.32 cm |
| distributed opening | 2.85 cm | 3.61 cm | 3.42 cm | 1.30 cm | 2.03 cm |
| distributed opening + graben faults | 2.61 cm | 3.56 cm | 2.88 cm | 1.29 cm | 2.01 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.180 km³** (opening-only variant 0.190); max opening **2.31 m** at 4.1 km centroid depth
- retained fault term: regularized maximum normal slip east 0.67 m; fault/dike potency ratio 0.107 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 3.04 cm → with faults 2.86 cm
- train parity 1: held-out RMS opening-only 3.12 cm → with faults 2.90 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.179 / 2.47 / 1.05 / 2.52
- smoothing λ×3.2: 0.196 / 1.90 / 0.40 / 2.82
- NISAR only: 0.209 / 2.27 / 0.75 / 3.10
- Sentinel-1 only: 0.181 / 2.00 / 0.15 / 1.52
- no decay normalization: 0.187 / 2.38 / 0.74 / 2.73

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
