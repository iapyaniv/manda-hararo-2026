# Along-rift (U∥) north/south-half averages on the flanks

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
| WEST north | 31 | +9.3 ± 2.7 (median +10.2) | +4.5 |
| WEST south | 28 | -5.5 ± 2.7 (median -5.8) | -4.5 |
| EAST north | 32 | +9.4 ± 2.7 (median +10.1) | +3.9 |
| EAST south | 17 | +4.1 ± 2.7 (median +10.8) | -2.0 |
| BOTH north | 63 | +9.3 ± 2.7 (median +10.2) | +4.2 |
| BOTH south | 45 | -1.8 ± 2.7 (median -2.3) | -3.6 |

- **WEST north − south: +14.7 ± 1.9 cm** (model +9.0; detection 7.8 σ, departure from the model +3.0 σ)
- **EAST north − south: +5.3 ± 1.9 cm** (model +5.9; detection 2.8 σ, departure from the model -0.3 σ)
- **BOTH north − south: +11.2 ± 1.9 cm** (model +7.8; detection 5.9 σ, departure from the model +1.8 σ)

- empirical errors: 34 far-field placebo boxes of the regions' dimensions (8.7 × 6.7 km) on the observed − model residual give a region-mean scatter of 2.7 cm and an adjacent-pair differential scatter of 1.9 cm (28 pairs); the formal independent-cell errors are several times smaller and are not quoted
- region sensitivity (|s| < 12 km, 1.3–6 km off axis): north − south WEST +17.1 cm (model +9.5); EAST +8.2 cm (model +6.0)
- strike sensitivity (148.5° from the southernmost 2005 segment): north − south WEST +14.4 cm; EAST +6.6 cm — the differential survives the strike change that moves individual probes by several cm

## Notes

- The free three-component solve observes the near-north direction weakly (median σ∥ ≈ 4 cm per cell), and the cells are not independent; the placebo boxes carry the correlated noise and are the errors to quote.
- A strike error leaks U⊥ into U∥ with one sign per flank, shifting both halves together, so the north-minus-south differential is the strike-robust statistic.
- The EAST south region is the raggedest (fewest cells; its mean and median disagree most) — it sits where the east-flank cross-sensor closure is poorest (see results/decomposition_joint*.md).
- The model column shares each region's exact cells, so coverage imbalance between halves affects observed and model alike.
