# Distributed opening + graben-fault model (four geometries; fault-trace geometry)

- fault geometry: anchored at trace columns from results/fault_traces.json (SHA-256 2bab6737597afb4633bee4b716ef8ed1cdf002b9b3c8edeadd3fef81a16cc09e); dip 65° inward; each column ends where it meets the dike plane
- trace evidence rule: west lacks a consistent LOS step across geometries and is excluded from this fit; the retained trace geometry remains a sensitivity, not a fault map
- data: NISAR asc 07-08→07-20, NISAR desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14053 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.49 cm | 2.96 cm | 2.70 cm | 1.70 cm | 2.32 cm |
| distributed opening | 2.37 cm | 2.83 cm | 2.84 cm | 1.32 cm | 1.88 cm |
| distributed opening + graben faults | 2.15 cm | 2.54 cm | 2.46 cm | 1.36 cm | 1.89 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.177 km³** (opening-only variant 0.186); max opening **2.33 m** at 4.1 km centroid depth
- retained fault term: regularized maximum normal slip east 0.67 m; fault/dike potency ratio 0.102 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.51 cm → with faults 2.32 cm
- train parity 1: held-out RMS opening-only 2.63 cm → with faults 2.38 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.175 / 2.37 / 0.80 / 2.05
- smoothing λ×3.2: 0.194 / 1.90 / 0.41 / 2.38
- NISAR only: 0.188 / 2.25 / 0.82 / 2.43
- Sentinel-1 only: 0.179 / 2.04 / 0.07 / 1.53
- trace columns 250 m inward: 0.179 / 2.39 / 0.60 / 2.22
- trace columns 250 m outward: 0.176 / 2.32 / 0.69 / 2.14
- fault dip 55°: 0.176 / 2.29 / 0.71 / 2.16
- fault dip 75°: 0.180 / 2.36 / 0.62 / 2.17
- no decay normalization: 0.185 / 2.39 / 0.74 / 2.24

## Dip sweep (trace-column anchoring; truncation depth scales with tan dip)

| dip | RMS all | potency ratio | east max slip |
|---|---|---|---|
| 55° | 2.1634 cm | 0.104 | 0.71 m |
| 60° | 2.1513 cm | 0.107 | 0.68 m |
| 65° | 2.1487 cm | 0.102 | 0.67 m |
| 70° | 2.1533 cm | 0.096 | 0.65 m |
| 75° | 2.1675 cm | 0.087 | 0.62 m |

## Scope

- The dike plane, strike, and fault geometry are fixed a priori; the
  inversion solves only dislocation magnitudes (convex, positive).
- The fault surfaces keep one fixed strike while their column
  offsets follow the supplied traces, so plan-view steps between
  columns are parameterisation steps, not resolved geometry; the
  smoothing couples stepped columns as neighbours. Columns end at
  the dike plane. A direct test against the 2005-2010 dike and
  fault geometries requires those geometry files.
- The retained east-fault term improves both sensor-subset fits, but
  its magnitude changes between them; treat its slip values as
  regularized magnitudes, not resolved slip. The west trace is absent
  because it fails the cross-geometry LOS-step support rule.
- Formal resolution/uncertainty of the slip distribution is not
  computed here; the sensitivity list brackets the main choices.
