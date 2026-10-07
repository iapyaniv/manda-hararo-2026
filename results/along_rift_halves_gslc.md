# Along-rift (U∥) north/south-half averages on the flanks (30-m reprocessed NISAR inputs)

Along-rift stretching moves the northern half of the deforming zone NNW and
the southern half SSE. The flank probes are local 2-km disks at the opening
peaks (3.8 and 1.9 km south of the fitted centre) and a full-flank average
integrates across the sign change, so neither measures that contrast.
These are half-flank averages of the dike-frame U∥ field (σ∥ ≤ 5 cm, the
summary-figure mask), against the fitted model passed through the same
processing (per-geometry LOS, per-dataset far-field reference plane, same
solve and rotation) over the same cells. The elastic dike itself predicts
the antisymmetric pattern, so the test is observed vs model, not observed
vs zero.

Regions: |s| < L/2 = 8.7 km from the fitted centre, 1.3–8 km off axis.

| region | n | observed (cm) | model (cm) |
|---|---|---|---|
| WEST north | 31 | +9.5 ± 2.6 (median +10.3) | +3.5 |
| WEST south | 29 | -4.2 ± 2.6 (median -5.1) | -5.1 |
| EAST north | 32 | +9.2 ± 2.6 (median +9.5) | +3.1 |
| EAST south | 17 | +5.5 ± 2.6 (median +10.6) | -2.6 |
| BOTH north | 63 | +9.4 (median +9.8) | +3.3 |
| BOTH south | 46 | -0.6 (median -0.6) | -4.2 |

- **WEST north − south: +13.7 ± 1.9 cm** (model +8.6); relative to the control floor: 7.1 σ from zero, +2.6 σ from the model; largest single-cell influence 0.8 cm
- **EAST north − south: +3.7 ± 1.9 cm** (model +5.7); relative to the control floor: 1.9 σ from zero, -1.1 σ from the model; largest single-cell influence 3.5 cm
- BOTH flanks north − south: +9.9 cm (model +7.5) — no control floor is assigned: the statistic spans the two flanks with unequal weights

- control floor: 34 far-field boxes of the regions' dimensions (8.7 × 6.7 km) on the observed − model residual give a region-mean scatter of 2.6 cm and an adjacent-pair differential scatter of 1.9 cm (28 pairs). 31 of 34 boxes lie east of the dike; median cross-sensor closure is 0.39 σ in the boxes vs 0.39/0.41 σ in the west regions and 1.50/2.40 σ in the east regions, so the floor represents the west-flank statistics and understates the east-flank errors; the formal independent-cell errors are several times smaller still and are not quoted
- region sensitivity (|s| < 12 km, 1.3–6 km off axis): north − south WEST +15.7 cm (model +9.0); EAST +6.3 cm (model +5.7)
- strike sensitivity (148.5° from the southernmost 2005 segment): north − south WEST +13.5 cm; EAST +5.1 cm — the differential survives the strike change that moves individual probes by several cm

## Notes

- The free three-component solve observes the near-north direction weakly (median σ∥ ≈ 4 cm per cell), and the cells are not independent; the control-box scatter is a correlated-noise floor, not a complete error model, so the σ statements above are relative to that floor. Model-parameter uncertainty and cross-sensor systematics are not propagated.
- The WEST differential is the robust result: the control floor's closure level matches the west regions, and no single cell changes the differential by more than 0.8 cm.
- The EAST differential is fragile and does not establish the reversal: both east halves average positive, the south region is the raggedest (fewest cells, mean and median farthest apart), one cell moves the differential by 3.5 cm, and the east-flank closure is well above the control level (the east-flank cross-sensor disagreement of results/decomposition_joint*.md).
- A strike error leaks U⊥ into U∥ with one sign per flank, shifting both halves together, so the north-minus-south differential is the strike-robust statistic.
- The model column shares each region's exact cells and the observations' reference treatment; the unreferenced elastic means (model_raw_mean_cm in the JSON) differ by up to about 1 cm per region.
