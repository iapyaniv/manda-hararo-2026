# Distributed opening + graben-fault model (four geometries; NISAR from GSLC reprocessing; complete external geometry)

- fault geometry: complete patch geometry from results/geometry_wright2006.json (file SHA-256 b643c0de55aa324e5683817d587cada1e9b6daebc4383849c45ef813ca5cbbee; geometry SHA-256 9e2429ca3944c943b193c461279ecab180e7facb47e9d9d1be4a046734e969c2)
- data: NISAR/GSLC asc 07-08→07-20, NISAR/GSLC desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14388 blocks)
- decay normalization to 2026-07-18 (factors 1.055, 1.031, 1.000, 1.102)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.87 cm | 3.41 cm | 3.37 cm | 1.70 cm | 2.32 cm |
| distributed opening | 2.85 cm | 3.61 cm | 3.42 cm | 1.30 cm | 2.03 cm |
| distributed opening + graben faults | 2.63 cm | 3.59 cm | 2.91 cm | 1.28 cm | 2.01 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.181 km³** (opening-only variant 0.190); max opening **2.33 m** at 4.1 km centroid depth
- retained fault term: regularized maximum normal slip east 0.66 m; fault/dike potency ratio 0.104 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 3.05 cm → with faults 2.90 cm
- train parity 1: held-out RMS opening-only 3.13 cm → with faults 2.91 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.180 / 2.50 / 1.13 / 2.54
- smoothing λ×3.2: 0.197 / 1.90 / 0.39 / 2.83
- NISAR only: 0.211 / 2.29 / 0.73 / 3.13
- Sentinel-1 only: 0.181 / 2.00 / 0.15 / 1.52
- no decay normalization: 0.188 / 2.40 / 0.72 / 2.74

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
