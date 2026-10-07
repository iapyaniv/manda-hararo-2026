# Along-rift (U∥) north/south-half averages on the flanks (30-m reprocessed NISAR inputs)

The flank probes sit near the dike's along-strike midpoint, the node of any
along-rift pattern: rift-parallel stretching moves the northern half NNW and
the southern half SSE, which a mid-dike probe or full-flank average cancels.
These are half-flank averages of the dike-frame U∥ field (σ∥ ≤ 5 cm, the
summary-figure mask), against the uniform-opening model over the same cells.
The elastic dike itself predicts the antisymmetric pattern, so the test is
observed vs model, not observed vs zero.

Regions: |s| < L/2 = 8.7 km from the fitted centre, 1.3–8 km off axis.

| region | n | observed (cm) | model (cm) |
|---|---|---|---|
| WEST north | 31 | +9.5 ± 2.8 (median +10.3) | +4.5 |
| WEST south | 29 | -4.2 ± 2.8 (median -5.1) | -4.3 |
| EAST north | 32 | +9.2 ± 2.8 (median +9.5) | +3.9 |
| EAST south | 17 | +5.5 ± 2.8 (median +10.6) | -2.0 |
| BOTH north | 63 | +9.4 ± 2.8 (median +9.8) | +4.2 |
| BOTH south | 46 | -0.6 ± 2.8 (median -0.6) | -3.5 |

- **WEST north − south: +13.7 ± 1.9 cm** (model +8.8; detection 7.1 σ, departure from the model +2.5 σ)
- **EAST north − south: +3.7 ± 1.9 cm** (model +5.9; detection 1.9 σ, departure from the model -1.2 σ)
- **BOTH north − south: +9.9 ± 1.9 cm** (model +7.7; detection 5.2 σ, departure from the model +1.2 σ)

- empirical errors: 34 far-field placebo boxes of the regions' dimensions (8.7 × 6.7 km) on the observed − model residual give a region-mean scatter of 2.8 cm and an adjacent-pair differential scatter of 1.9 cm (28 pairs); the formal independent-cell errors are several times smaller and are not quoted
- region sensitivity (|s| < 12 km, 1.3–6 km off axis): north − south WEST +15.7 cm (model +9.3); EAST +6.3 cm (model +6.0)
- strike sensitivity (148.5° from the southernmost 2005 segment): north − south WEST +13.5 cm; EAST +5.1 cm — the differential survives the strike change that moves individual probes by several cm

## Notes

- The free three-component solve observes the near-north direction weakly (median σ∥ ≈ 4 cm per cell), and the cells are not independent; the placebo boxes carry the correlated noise and are the errors to quote.
- A strike error leaks U⊥ into U∥ with one sign per flank, shifting both halves together, so the north-minus-south differential is the strike-robust statistic.
- The EAST south region is the raggedest (fewest cells; its mean and median disagree most) — it sits where the east-flank cross-sensor closure is poorest (see results/decomposition_joint*.md).
- The model column shares each region's exact cells, so coverage imbalance between halves affects observed and model alike.
