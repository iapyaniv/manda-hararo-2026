# Distributed opening + graben-fault model (four geometries; NISAR from GSLC reprocessing; fault-trace geometry)

- fault geometry: anchored at trace columns from results/fault_traces.json (SHA-256 76328d6ce925d064359d0fd2befe02a7aa02084e55fe258f152f891fc18dc5cf); dip 65° inward; each column ends where it meets the dike plane
- trace evidence rule: west lacks a consistent LOS step across geometries and is excluded from this fit; the retained trace geometry remains a sensitivity, not a fault map
- data: NISAR/GSLC asc 07-08→07-20, NISAR/GSLC desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14388 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 3.07 cm | 3.48 cm | 3.52 cm | 2.25 cm | 2.56 cm |
| distributed opening | 2.98 cm | 3.49 cm | 3.61 cm | 1.91 cm | 2.09 cm |
| distributed opening + graben faults | 2.80 cm | 3.45 cm | 3.18 cm | 1.93 cm | 2.08 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.152 km³** (opening-only variant 0.161); max opening **2.52 m** at 3.6 km centroid depth
- retained fault term: regularized maximum normal slip east 0.61 m; fault/dike potency ratio 0.110 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 3.13 cm → with faults 2.99 cm
- train parity 1: held-out RMS opening-only 3.41 cm → with faults 3.21 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.145 / 2.76 / 1.03 / 2.70
- smoothing λ×3.2: 0.169 / 1.88 / 0.38 / 3.06
- NISAR only: 0.172 / 2.42 / 0.71 / 3.24
- Sentinel-1 only: 0.149 / 2.12 / 0.07 / 1.92
- trace columns 250 m inward: 0.154 / 2.55 / 0.56 / 2.83
- trace columns 250 m outward: 0.151 / 2.49 / 0.65 / 2.77
- fault dip 55°: 0.153 / 2.54 / 0.53 / 2.86
- fault dip 75°: 0.154 / 2.51 / 0.67 / 2.77
- no decay normalization: 0.158 / 2.59 / 0.68 / 2.92

## Dip sweep (trace-column anchoring; truncation depth scales with tan dip)

| dip | RMS all | potency ratio | east max slip |
|---|---|---|---|
| 55° | 2.8595 cm | 0.087 | 0.53 m |
| 60° | 2.8257 cm | 0.102 | 0.57 m |
| 65° | 2.7982 cm | 0.110 | 0.61 m |
| 70° | 2.7773 cm | 0.115 | 0.65 m |
| 75° | 2.7672 cm | 0.113 | 0.67 m |

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
