# Loop-closure check of the reprocessed NISAR unwrapping

Per track, the three epochs give three interferograms on the shared GSLC grid
(pre, co-event, and the skip pair that spans both). The unwrapped closure
`unw_pre + unw_co − unw_skip` minus the wrapped closure of the same filtered
phases is an exact integer number of cycles per pixel (snaphu is congruent
modulo 2π). Deformation, atmosphere, and ionosphere are epoch differences and
cancel in the loop, so a nonzero integer isolates unwrap inconsistency among
the three legs. The arbitrary per-unwrap reference makes a constant integer
offset unobservable; the far-field mode is removed. One cycle is λ/2 = 11.9 cm of LOS. This is an internal-consistency test
of the three unwraps against each other, not an accuracy test of any one of
them. The co-event phase inputs are the saved analysis interferograms.
The raw caches are checked for phase congruence; their snaphu settings are reported only when cache metadata records them.
The near-axis band has a 3-km half-width (6 km wide) and a 20-km along-axis half-length.

## 172_008_A

- triangle: pre 2026-06-14 -> 2026-07-08, co 2026-07-08 -> 2026-07-20, skip 2026-06-14 -> 2026-07-20
- cached snaphu settings: pre unknown (no config metadata); co unknown (no config metadata); skip unknown (no config metadata)
- cached snaphu executable identity and version: pre unknown; co unknown; skip unknown
- maximum raw/phase congruence error: pre 9.83e-06 rad, co 5.61e-05 rad, skip 6.66e-05 rad
- whole coherent scene: 0.78% (40,067 of 5,166,670 pixels) off by ≥ 1 cycle (36.06 km²; max 4 cycles); wrapped nonclosure RMS 0.007 cycles
- near-axis band (|⊥| < 3 km, |∥| < 20 km): 16.24% (38,381 of 236,307 pixels) off by ≥ 1 cycle (34.54 km²; max 4 cycles); wrapped nonclosure RMS 0.016 cycles
- off-axis within 25 km: 0.0004% (5 of 1,297,072 pixels) off by ≥ 1 cycle (0.0045 km²; max 1 cycles); wrapped nonclosure RMS 0.007 cycles
- error-region attribution: median |⊥ distance| 0.9 km; median coherence there pre 0.67, co 0.69, skip 0.58 (elsewhere 0.86, 0.87, 0.85) — the loop flags the triple, and cannot assign the error to one leg by itself

## 165_082_D

- triangle: pre 2026-06-13 -> 2026-07-07, co 2026-07-07 -> 2026-07-19, skip 2026-06-13 -> 2026-07-19
- cached snaphu settings: pre unknown (no config metadata); co unknown (no config metadata); skip unknown (no config metadata)
- cached snaphu executable identity and version: pre unknown; co unknown; skip unknown
- maximum raw/phase congruence error: pre 8.89e-06 rad, co 6.36e-05 rad, skip 7.08e-05 rad
- whole coherent scene: 0.05% (3,928 of 8,240,337 pixels) off by ≥ 1 cycle (3.54 km²; max 3 cycles); wrapped nonclosure RMS 0.008 cycles
- near-axis band (|⊥| < 3 km, |∥| < 20 km): 1.13% (2,634 of 233,457 pixels) off by ≥ 1 cycle (2.37 km²; max 3 cycles); wrapped nonclosure RMS 0.017 cycles
- off-axis within 25 km: 0.0000% (0 of 1,741,239 pixels) off by ≥ 1 cycle (0.0000 km²; max 0 cycles); wrapped nonclosure RMS 0.008 cycles
- error-region attribution: median |⊥ distance| 1.6 km; median coherence there pre 0.52, co 0.49, skip 0.43 (elsewhere 0.80, 0.74, 0.85) — the loop flags the triple, and cannot assign the error to one leg by itself

## 100_A_008_STRADDLE

- triangle: pre 2026-07-03 -> 2026-07-15, co 2026-07-15 -> 2026-07-27, skip 2026-07-03 -> 2026-07-27
- cached snaphu settings: pre unknown (no config metadata); co unknown (no config metadata); skip unknown (no config metadata)
- cached snaphu executable identity and version: pre unknown; co unknown; skip unknown
- maximum raw/phase congruence error: pre 1.84e-06 rad, co 4.45e-05 rad, skip 5.73e-05 rad
- whole coherent scene: 0.30% (12,044 of 3,975,077 pixels) off by ≥ 1 cycle (10.84 km²; max 3 cycles); wrapped nonclosure RMS 0.008 cycles
- near-axis band (|⊥| < 3 km, |∥| < 20 km): 4.70% (11,082 of 235,580 pixels) off by ≥ 1 cycle (9.97 km²; max 3 cycles); wrapped nonclosure RMS 0.021 cycles
- off-axis within 25 km: 0.0039% (52 of 1,317,300 pixels) off by ≥ 1 cycle (0.0468 km²; max 1 cycles); wrapped nonclosure RMS 0.007 cycles
- error-region attribution: median |⊥ distance| 1.2 km; median coherence there pre 0.72, co 0.51, skip 0.45 (elsewhere 0.87, 0.86, 0.83) — the loop flags the triple, and cannot assign the error to one leg by itself

## Scope

- `loop_closure_inventory.json` records SHA-256 hashes, file sizes, analysis settings, source files, and report identities. Its 18 interferogram/cache inputs reproduce closure with `python -m mh26.closure --check --cached-only`, without snaphu or file writes. Replay first verifies the saved source, settings, reports, and input hashes, then checks each interferogram's ordered granules, EPSG, looks, array shape, 30-m grid, and raw-array phase congruence. Explicit `--regenerate` replaces reports and the inventory after an authorized analysis change; input hashes remain pinned. The original settings of legacy caches are unknown; the current default is not evidence of those settings. Missing optional raw inputs prevent this replay even when the archive checks pass.
- The mask requires coherence above the analysis threshold (0.25) in all three legs, so the checked area is the intersection. Flagged pixels therefore pass that threshold; coherence screening does not remove every unwrap inconsistency. This check does not show that all flagged axial support was excluded from the model fits.
- The wrapped-closure RMS is the multilook/filter nonclosure floor; it is not an unwrap error.
- Closure inconsistency is nearly zero outside the near-axis band: less than 0.005% on each track within 25 km. Inside the band, inconsistencies reach 1–4 cycles (about 12–48 cm LOS). The loop cannot identify the faulty leg or certify the flank or axial displacement accuracy.
- A triple can contain unwrapping errors that cancel in the loop. The GUNW comparison (results/gslc_gunw_compare_*.json) and the Sentinel-1 cross-checks add evidence, but do not exclude every consistent error.
- The pre legs span 24 days on tracks 172/165 (the 12-day cycle-023 GSLC epochs are absent from the provisional archive) and 12 days on track 100.
