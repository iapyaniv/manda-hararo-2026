# Distributed opening + graben-fault model (four geometries; complete external geometry)

- fault geometry: complete patch geometry from results/geometry_wright2006.json (file SHA-256 f9c9032796101047f2c5e0ca0292d766bcd5ac300d4ed49ce3677c29f675f3d6; geometry SHA-256 48ee382a4f6e0431229308cf5e896b64b1fe0a31a4eb6ce57c95bad0447c5051)
- data: NISAR asc 07-08→07-20, NISAR desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14053 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.49 cm | 2.96 cm | 2.70 cm | 1.70 cm | 2.32 cm |
| distributed opening | 2.43 cm | 3.21 cm | 2.69 cm | 1.27 cm | 2.01 cm |
| distributed opening + graben faults | 2.20 cm | 2.96 cm | 2.26 cm | 1.28 cm | 1.99 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.182 km³** (opening-only variant 0.191); max opening **2.30 m** at 4.1 km centroid depth
- retained fault term: regularized maximum normal slip east 0.71 m; fault/dike potency ratio 0.107 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.57 cm → with faults 2.35 cm
- train parity 1: held-out RMS opening-only 2.73 cm → with faults 2.46 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.181 / 2.39 / 0.92 / 2.10
- smoothing λ×3.2: 0.197 / 1.88 / 0.44 / 2.44
- NISAR only: 0.215 / 2.23 / 0.83 / 2.49
- Sentinel-1 only: 0.181 / 2.00 / 0.15 / 1.52
- no decay normalization: 0.189 / 2.37 / 0.78 / 2.30

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
