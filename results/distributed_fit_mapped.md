# Distributed opening + graben-fault model (four geometries; fault-trace geometry)

- fault geometry: anchored at trace columns from results/fault_traces.json (SHA-256 76328d6ce925d064359d0fd2befe02a7aa02084e55fe258f152f891fc18dc5cf); dip 65° inward; each column ends where it meets the dike plane
- trace evidence rule: west lacks a consistent LOS step across geometries and is excluded from this fit; the retained trace geometry remains a sensitivity, not a fault map
- data: NISAR asc 07-08→07-20, NISAR desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14053 blocks)
- decay normalization to 2026-07-18 (factors 1.057, 1.032, 1.000, 1.105)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.71 cm | 3.01 cm | 2.89 cm | 2.25 cm | 2.56 cm |
| distributed opening | 2.59 cm | 2.96 cm | 3.02 cm | 1.90 cm | 2.08 cm |
| distributed opening + graben faults | 2.42 cm | 2.74 cm | 2.69 cm | 1.95 cm | 2.07 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.153 km³** (opening-only variant 0.162); max opening **2.46 m** at 3.6 km centroid depth
- retained fault term: regularized maximum normal slip east 0.67 m; fault/dike potency ratio 0.112 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.73 cm → with faults 2.62 cm
- train parity 1: held-out RMS opening-only 2.84 cm → with faults 2.63 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.146 / 2.72 / 0.98 / 2.32
- smoothing λ×3.2: 0.169 / 1.87 / 0.39 / 2.67
- NISAR only: 0.174 / 2.32 / 0.80 / 2.66
- Sentinel-1 only: 0.149 / 2.12 / 0.07 / 1.92
- trace columns 250 m inward: 0.155 / 2.49 / 0.61 / 2.46
- trace columns 250 m outward: 0.151 / 2.43 / 0.69 / 2.38
- fault dip 55°: 0.153 / 2.43 / 0.70 / 2.44
- fault dip 75°: 0.156 / 2.50 / 0.64 / 2.43
- no decay normalization: 0.159 / 2.53 / 0.73 / 2.52

## Dip sweep (trace-column anchoring; truncation depth scales with tan dip)

| dip | RMS all | potency ratio | east max slip |
|---|---|---|---|
| 55° | 2.4397 cm | 0.104 | 0.70 m |
| 60° | 2.4249 cm | 0.112 | 0.69 m |
| 65° | 2.4192 cm | 0.112 | 0.67 m |
| 70° | 2.4190 cm | 0.107 | 0.66 m |
| 75° | 2.4269 cm | 0.100 | 0.64 m |

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
