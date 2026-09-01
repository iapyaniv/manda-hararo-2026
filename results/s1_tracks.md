# Two-track Sentinel-1 tests

## Descending co-event cross-check

- LiCSAR 079D_07694_131313 pair 20260711_20260723 (2026-07-11→2026-07-23): **r = 0.915**, amplitude ratio 0.81, RMS 2.9 cm against the NISAR-only model (44529 pixels)
- ascending frame, same statistic (consistency with `mh26.s1_check`): r = 0.953, ratio 0.85, RMS 2.4 cm
- empirical LiCSAR sign factors: ascending -1, descending -1, resolved from the co-event pairs only

## Dike-pattern projections

| geometry | pair | days | role | blocks | projection | cluster-wild 95% | zero-intercept | r |
|---|---|---:|---|---:|---:|---|---:|---:|
| ascending | 2026-07-12→2026-07-18 | 6 | co-event | 1024 | **+0.900** | [+0.872, +0.929] | +0.889 | +0.98 |
| descending | 2026-07-11→2026-07-23 | 12 | co-event | 916 | **+0.891** | [+0.820, +0.965] | +0.855 | +0.95 |
| ascending | 2026-07-18→2026-07-30 | 12 | post-event | 1020 | **+0.136** | [+0.104, +0.168] | +0.117 | +0.81 |
| ascending | 2026-07-30→2026-08-11 | 12 | post-event | 1056 | **+0.016** | [-0.004, +0.034] | +0.004 | +0.20 |
| ascending | 2026-07-18→2026-08-11 | 24 | post-event | 958 | **+0.151** | [+0.122, +0.179] | +0.122 | +0.81 |
| descending | 2026-07-23→2026-08-04 | 12 | post-event | 1008 | **+0.052** | [+0.043, +0.061] | +0.048 | +0.73 |
| descending | 2026-08-04→2026-08-16 | 12 | post-event | 994 | **-0.015** | [-0.023, -0.007] | -0.030 | -0.30 |
| descending | 2026-07-23→2026-08-16 | 24 | post-event | 893 | **+0.039** | [+0.028, +0.051] | +0.016 | +0.58 |
| ascending | 2026-06-05→2026-06-12 | 7 | pre-event | 1068 | **-0.021** | [-0.032, -0.010] | -0.013 | -0.27 |
| ascending | 2026-06-12→2026-06-17 | 5 | pre-event | 1070 | **+0.020** | [+0.010, +0.030] | +0.049 | +0.29 |
| ascending | 2026-06-17→2026-06-24 | 7 | pre-event | 1069 | **-0.016** | [-0.029, -0.004] | -0.024 | -0.20 |
| ascending | 2026-06-24→2026-06-29 | 5 | pre-event | 1072 | **+0.014** | [+0.006, +0.022] | +0.007 | +0.23 |
| ascending | 2026-06-29→2026-07-06 | 7 | pre-event | 1072 | **+0.007** | [+0.003, +0.010] | +0.028 | +0.20 |
| ascending | 2026-07-06→2026-07-12 | 6 | pre-event | 1074 | **-0.021** | [-0.031, -0.012] | -0.052 | -0.37 |
| descending | 2026-06-10→2026-06-22 | 12 | pre-event | 1019 | **+0.008** | [-0.002, +0.019] | +0.012 | +0.16 |
| descending | 2026-06-22→2026-06-29 | 7 | pre-event | 1027 | **+0.011** | [+0.006, +0.017] | +0.020 | +0.27 |
| descending | 2026-06-29→2026-07-11 | 12 | pre-event | 1026 | **+0.004** | [-0.003, +0.012] | -0.006 | +0.10 |

- 9 pre-event null projections: mean +0.001, standard deviation 0.016, maximum magnitude 0.021
- the two early post-event windows exceed every null in both
  geometries; the two late windows (from 30 July and 4 August) are at
  the null level: the signal decays
- sub-window additivity: ascending +0.136 + +0.016 = +0.151 vs +0.151 direct; descending +0.052 + -0.015 = +0.037 vs +0.039 direct

## Post-event decay

- exponential window-integral fit to the six Sentinel-1 projections: total **+0.149** of the co-event pattern, decay constant **4.3 days** from 18 July (fit RMS 0.009); a descriptive summary, not a physical model
- independent NISAR descending pair 2026-07-19→2026-08-12: observed +0.107; the Sentinel-1-only decay curve predicts +0.118 for that window

## Post-event East/Up decomposition

### co-event (ascending 2026-07-12→2026-07-18, descending 2026-07-11→2026-07-23)

- flank quasi-east probes (2-km disks): EAST +54.4 cm, WEST -57.7 cm ⇒ **opening divergence +112.1 cm**
- axis quasi-up probe: decorrelated
- east p01–p99 -56.1 to +52.5 cm; up p01–p99 -2.0 to +26.2 cm (502 cells within 18 km)

### post-event (ascending 2026-07-18→2026-08-11, descending 2026-07-23→2026-08-16)

- flank quasi-east probes (2-km disks): EAST +6.9 cm, WEST -11.1 cm ⇒ **opening divergence +18.0 cm**
- axis quasi-up probe: -1.0 cm
- east p01–p99 -9.9 to +7.7 cm; up p01–p99 -4.3 to +5.7 cm (519 cells within 18 km)

## Scope

- The projection detects only deformation shaped like the fitted dike;
  the intercept absorbs a broad datum shift per pair.
- The ascending and descending post-event windows differ by five days at
  each end; the decomposition assumes the signal is common to both. For
  a decaying signal that assumption fails in part, and the East/Up split
  blends the two windows; the projection table is the quantitative result.
- The cluster-wild intervals and the pre-event null distribution do not
  form a calibrated atmospheric covariance model.
- Near-axis C-band decoherence removes the highest-gradient pixels from
  the co-event comparisons.
