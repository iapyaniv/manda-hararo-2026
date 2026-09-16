# Graben-fault depth-extent sweep (mapped trace; GUNW product chain)

- reference fit: distributed_fit_mapped.json (λ = 17.78 reused; SHA-256 b67920d372db…); the uncapped baseline must reproduce it exactly before the sweep runs
- graben band: 6-km half-width from the dike axis (12 km wide), 20-km along-axis half-length (991 blocks); the fault term acts here, and the all-block RMS dilutes it
- fault moment = μ · Σ area·slip at μ = 30 GPa. Okada surface displacements depend on (λ, μ) only through α = (λ+μ)/(λ+2μ), so a uniform rescaling of both changes no prediction and rescales the moment linearly: the data do not constrain μ. A depth-varying soft layer has not been modelled.

| fault depth extent | patches | max slip | potency | μ·potency | dike volume | RMS all | RMS band | band improvement retained |
|---|---|---|---|---|---|---|---|---|
| none (opening only) | 0 | — | — | — | 0.186 km³ | 2.37 cm | 7.64 cm | 0.00 |
| ≤ 0.5 km | 14 | 0.34 m | 1.7×10⁶ m³ | 0.5×10¹⁷ N·m | 0.186 km³ | 2.33 cm | 7.48 cm | 0.18 |
| ≤ 1.0 km | 14 | 0.47 m | 4.7×10⁶ m³ | 1.4×10¹⁷ N·m | 0.184 km³ | 2.27 cm | 7.23 cm | 0.46 |
| ≤ 1.5 km | 25 | 0.62 m | 8.8×10⁶ m³ | 2.6×10¹⁷ N·m | 0.182 km³ | 2.21 cm | 7.01 cm | 0.72 |
| ≤ 2.0 km | 25 | 0.60 m | 11.2×10⁶ m³ | 3.4×10¹⁷ N·m | 0.181 km³ | 2.19 cm | 6.89 cm | 0.85 |
| ≤ 2.5 km | 25 | 0.60 m | 13.3×10⁶ m³ | 4.0×10¹⁷ N·m | 0.180 km³ | 2.17 cm | 6.82 cm | 0.93 |
| ≤ 3.0 km | 33 | 0.74 m | 16.5×10⁶ m³ | 4.9×10¹⁷ N·m | 0.178 km³ | 2.15 cm | 6.75 cm | 1.01 |
| ≤ 3.5 km | 33 | 0.71 m | 17.4×10⁶ m³ | 5.2×10¹⁷ N·m | 0.177 km³ | 2.15 cm | 6.75 cm | 1.01 |
| trace-limited (≤ 4 km) | 33 | 0.67 m | 18.1×10⁶ m³ | 5.4×10¹⁷ N·m | 0.177 km³ | 2.15 cm | 6.76 cm | 1.00 |

## Depth reach of the uncapped solution

- cumulative potency fraction above 1 km: 0.27; above 2 km: 0.56; above 3 km: 0.82 (piecewise-linear in depth within patches)
- shallowest tested cap retaining ≥ 90% of the fault term's graben-band improvement: 2.5 km, with potency 13.3×10⁶ m³ (μ·potency 4.0×10¹⁷ N·m, vs 5.4×10¹⁷ N·m uncapped)

## Scope

- The sweep truncates the depth extent of the fixed mapped-trace, fixed-dip fault surface and re-solves at the reused smoothing weight. The retained-improvement threshold is a comparison rule, not a confidence level or a data-required depth. The sweep does not test other fault geometries (listric surfaces, distributed fissuring, the unsupported west trace), and the slip values remain regularized magnitudes. Neither the capped nor the uncapped moment is a physical upper confidence bound.
- Truncation changes patch dimensions, Green's functions, and the retained smoothing operator. These are not nested fits on one fixed basis, so a capped fit can have lower RMS; some saved rows improve slightly on the uncapped fit. The dike-volume column records the change in the opening/slip partition.
- Lowering the uniform μ rescales the moment only if both Lamé constants are scaled together. A lower stiffness value is an assumed sensitivity, not a value inferred from these data. A depth-varying soft layer would also change the displacement calculation and has not been modelled.
- The graben band includes axial residuals and support where loop closure flags 1–4-cycle inconsistencies (about 12–48 cm LOS). The check cannot identify the faulty leg or show that all flagged axial support was excluded from the fits. Both the RMS levels and their differences are conditional on the retained observations and model assumptions.
- Both input chains run the same sweep (fault_depth_mapped.md / fault_depth_mapped_gslc.md).
