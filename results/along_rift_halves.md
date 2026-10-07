# Along-rift (U∥) north/south-half averages on the flanks

Along-rift stretching moves the northern half of the deforming zone NNW and
the southern half SSE. The flank probes are local 2-km disks at the opening
peaks (3.8 and 1.9 km south of the fitted centre) and a full-flank average
integrates across the sign change, so neither measures that contrast.
These are half-flank averages of the dike-frame U∥ field (σ∥ ≤ 5 cm, the
summary-figure mask), against the fitted model referenced like the data
(per-geometry LOS, per-dataset far-field reference plane, same solve and
rotation; the plane is fitted at cell resolution, an approximation worth
about 0.2 cm in the contrasts) over the same cells. The elastic dike itself
predicts the antisymmetric pattern, so the test is observed vs model, not
observed vs zero.

Regions: |s| < L/2 = 8.7 km from the fitted centre, 1.3–8 km off axis.

| region | n | observed (cm) | model (cm) |
|---|---|---|---|
| WEST north | 31 | +9.3 ± 2.5 (median +10.2) | +3.5 |
| WEST south | 28 | -5.5 ± 2.5 (median -5.8) | -5.3 |
| EAST north | 32 | +9.4 ± 2.5 (median +10.1) | +3.1 |
| EAST south | 17 | +4.1 ± 2.5 (median +10.8) | -2.6 |
| BOTH north | 63 | +9.3 (median +10.2) | +3.3 |
| BOTH south | 45 | -1.8 (median -2.3) | -4.3 |

- **WEST north − south: +14.7 ± 1.9 cm** (model +8.8); as ratios to the control floor (not established significances): 7.7× from zero, +3.1× from the model; largest single-cell influence 0.4 cm
- **EAST north − south: +5.3 ± 1.9 cm** (model +5.7); as ratios to the control floor (not established significances): 2.8× from zero, -0.2× from the model; largest single-cell influence 5.1 cm
- BOTH flanks north − south: +11.2 cm (model +7.6) — no control floor is assigned: the statistic spans the two flanks with unequal weights

- control floor: 34 far-field boxes of the regions' dimensions (8.7 × 6.7 km) on the observed − model residual give a region-mean scatter of 2.5 cm and an adjacent-pair differential scatter of 1.9 cm (28 pairs). 31 of 34 boxes lie east of the dike; median cross-sensor closure is 0.37 σ in the boxes vs 0.33/0.40 σ in the west regions and 1.44/2.48 σ in the east regions. Closure measures cross-sensor consistency only (an error consistent across the four looks leaves it unchanged), so the match makes the floor plausible for the west flank without certifying it, and east-flank errors certainly exceed the floor; the formal independent-cell errors are several times smaller still and are not quoted
- region sensitivity (|s| < 12 km, 1.3–6 km off axis): north − south WEST +17.1 cm (model +9.2); EAST +8.2 cm (model +5.7)
- strike sensitivity (148.5° from the southernmost 2005 segment): north − south WEST +14.4 cm; EAST +6.6 cm — the differential survives the strike change that moves individual probes by several cm

## Notes

- The free three-component solve observes the near-north direction weakly (median σ∥ ≈ 4 cm per cell), and the cells are not independent; the control-box scatter is a correlated-noise floor, not a demonstrated error model, so the ratios above are to that floor, not established significances. Model-parameter uncertainty and cross-sensor systematics are not propagated.
- The WEST north − south contrast is the stable statistic: no single cell changes it by more than 0.4 cm and the alternative strike moves it by 0.2 cm. The individual half means are NOT strike-stable (the south mean is -5.5 cm at the fitted strike and -1.5 cm at 148.5°), so the contrast — not the southward motion of the south half — is the resolved quantity.
- The EAST differential is fragile and does not establish the reversal: both east halves average positive, the south region is the raggedest (fewest cells, mean and median farthest apart), one cell moves the differential by 5.1 cm, and the east-flank closure is well above the control level (the east-flank cross-sensor disagreement of results/decomposition_joint*.md).
- A strike error leaks U⊥ into U∥ with one sign per flank, shifting both halves together, so the north-minus-south differential is the strike-robust statistic.
- The model column shares each region's exact cells and the observations' reference treatment at cell resolution (an approximation to the native-pixel referencing worth about 0.2 cm in the contrasts); the unreferenced elastic means (model_raw_mean_cm in the JSON) differ by up to about 1 cm per region.
