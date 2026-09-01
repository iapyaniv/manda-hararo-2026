# Distributed opening + graben-fault model (four geometries; NISAR from GSLC reprocessing; fault-trace geometry)

- fault geometry: anchored at trace columns from results/fault_traces.json (SHA-256 2bab6737597afb4633bee4b716ef8ed1cdf002b9b3c8edeadd3fef81a16cc09e); dip 65° inward; each column ends where it meets the dike plane
- trace evidence rule: west lacks a consistent LOS step across geometries and is excluded from this fit; the retained trace geometry remains a sensitivity, not a fault map
- data: NISAR/GSLC asc 07-08→07-20, NISAR/GSLC desc 07-07→07-19, S1 asc 07-12→07-18, S1 desc 07-11→07-23 (14388 blocks)
- decay normalization to 2026-07-18 (factors 1.055, 1.031, 1.000, 1.102)
- smoothing: λ = 17.78 reused from the primary fit so variants compare at matched regularization

| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |
|---|---|---|---|---|---|
| uniform rectangle (fixed, offsets refit) | 2.87 cm | 3.41 cm | 3.37 cm | 1.70 cm | 2.32 cm |
| distributed opening | 2.76 cm | 3.31 cm | 3.44 cm | 1.35 cm | 1.90 cm |
| distributed opening + graben faults | 2.56 cm | 3.26 cm | 2.98 cm | 1.36 cm | 1.89 cm |

- the fault term's in-sample test is the distributed-opening row against the with-faults row; the uniform row is a fixed reference with many fewer parameters, not a fault test

- dike volume **0.176 km³** (opening-only variant 0.185); max opening **2.38 m** at 4.1 km centroid depth
- retained fault term: regularized maximum normal slip east 0.61 m; fault/dike potency ratio 0.097 (Σ area·slip / Σ area·opening; not a moment-tensor ratio)

## Spatial holdout (5-km checkerboard)

- train parity 0: held-out RMS opening-only 2.92 cm → with faults 2.76 cm
- train parity 1: held-out RMS opening-only 3.15 cm → with faults 2.93 cm

## Sensitivities (volume km³ / max opening m / east max slip m / RMS cm)

- smoothing λ×1/3.2: 0.173 / 2.54 / 0.89 / 2.46
- smoothing λ×3.2: 0.194 / 1.91 / 0.37 / 2.80
- NISAR only: 0.187 / 2.35 / 0.71 / 3.03
- Sentinel-1 only: 0.179 / 2.04 / 0.08 / 1.54
- trace columns 250 m inward: 0.178 / 2.42 / 0.53 / 2.60
- trace columns 250 m outward: 0.175 / 2.37 / 0.64 / 2.54
- fault dip 55°: 0.176 / 2.41 / 0.52 / 2.62
- fault dip 75°: 0.178 / 2.38 / 0.65 / 2.53
- no decay normalization: 0.183 / 2.45 / 0.68 / 2.66

## Dip sweep (trace-column anchoring; truncation depth scales with tan dip)

| dip | RMS all | potency ratio | east max slip |
|---|---|---|---|
| 55° | 2.6222 cm | 0.082 | 0.52 m |
| 60° | 2.5872 cm | 0.092 | 0.56 m |
| 65° | 2.5578 cm | 0.097 | 0.61 m |
| 70° | 2.5372 cm | 0.099 | 0.64 m |
| 75° | 2.5288 cm | 0.095 | 0.65 m |

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
