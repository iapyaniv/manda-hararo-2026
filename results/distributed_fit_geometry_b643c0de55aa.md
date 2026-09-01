# Distributed opening + graben-fault model (four geometries; complete external geometry)

- fault geometry: complete patch geometry from results/geometry_wright2006.json (file SHA-256 b643c0de55aa324e5683817d587cada1e9b6daebc4383849c45ef813ca5cbbee; geometry SHA-256 9e2429ca3944c943b193c461279ecab180e7facb47e9d9d1be4a046734e969c2)
- data: NISAR asc 07-08→07-20, NISAR desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14053 blocks)
- decay normalization to 2026-07-18 (factors 1.055, 1.031, 1.000, 1.102)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.49 cm | 2.96 cm | 2.70 cm | 1.70 cm | 2.32 cm |
| distributed opening | 2.43 cm | 3.22 cm | 2.70 cm | 1.27 cm | 2.01 cm |
| distributed opening + graben faults | 2.24 cm | 3.05 cm | 2.30 cm | 1.28 cm | 1.99 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.183 km³** (opening-only variant 0.192); max opening **2.34 m** at 4.1 km centroid depth
- retained fault term: regularized maximum normal slip east 0.69 m; fault/dike potency ratio 0.099 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.57 cm → with faults 2.43 cm
- train parity 1: held-out RMS opening-only 2.73 cm → with faults 2.48 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.182 / 2.42 / 1.02 / 2.13
- smoothing λ×3.2: 0.198 / 1.89 / 0.41 / 2.47
- NISAR only: 0.218 / 2.28 / 0.79 / 2.56
- Sentinel-1 only: 0.181 / 2.00 / 0.15 / 1.52
- no decay normalization: 0.190 / 2.41 / 0.75 / 2.34

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
