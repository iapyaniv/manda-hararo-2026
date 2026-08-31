# Distributed opening + graben-fault model (four geometries; NISAR from GSLC reprocessing; complete external geometry)

- fault geometry: complete patch geometry from results/geometry_wright2006.json (file SHA-256 b643c0de55aa324e5683817d587cada1e9b6daebc4383849c45ef813ca5cbbee; geometry SHA-256 9e2429ca3944c943b193c461279ecab180e7facb47e9d9d1be4a046734e969c2)
- data: NISAR/GSLC asc 07-08→07-20, NISAR/GSLC desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14388 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 3.07 cm | 3.48 cm | 3.52 cm | 2.25 cm | 2.56 cm |
| distributed opening | 3.03 cm | 3.76 cm | 3.52 cm | 1.85 cm | 2.21 cm |
| distributed opening + graben faults | 2.81 cm | 3.71 cm | 3.03 cm | 1.86 cm | 2.19 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.154 km³** (opening-only variant 0.164); max opening **2.42 m** at 3.6 km centroid depth
- retained fault term: regularized maximum normal slip east 0.68 m; fault/dike potency ratio 0.126 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 3.21 cm → with faults 3.07 cm
- train parity 1: held-out RMS opening-only 3.29 cm → with faults 3.08 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.147 / 2.58 / 1.16 / 2.72
- smoothing λ×3.2: 0.170 / 1.87 / 0.40 / 3.04
- NISAR only: 0.177 / 2.32 / 0.75 / 3.27
- Sentinel-1 only: 0.150 / 2.08 / 0.11 / 1.89
- no decay normalization: 0.160 / 2.49 / 0.75 / 2.94

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
