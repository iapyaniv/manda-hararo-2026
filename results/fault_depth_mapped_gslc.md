# Graben-fault depth-extent sweep (mapped trace; 30-m reprocessed (GSLC) chain)

- reference fit: distributed_fit_mapped_gslc.json (λ = 17.78 reused; SHA-256 eebae6c090fa…); the uncapped baseline must reproduce it exactly before the sweep runs
- graben band: 6-km half-width from the dike axis (12 km wide), 20-km along-axis half-length (1011 blocks); the fault term acts here, and the all-block RMS dilutes it
- fault moment = μ · Σ area·slip at μ = 30 GPa. Okada surface displacements depend on (λ, μ) only through α = (λ+μ)/(λ+2μ), so a uniform rescaling of both changes no prediction and rescales the moment linearly: the data do not constrain μ. A depth-varying soft layer has not been modelled.

| fault depth extent | patches | max slip | potency | μ·potency | dike volume | RMS all | RMS band | band improvement retained |
|---|---|---|---|---|---|---|---|---|
| none (opening only) | 0 | — | — | — | 0.185 km³ | 2.76 cm | 9.25 cm | 0.00 |
| ≤ 0.5 km | 14 | 0.24 m | 1.4×10⁶ m³ | 0.4×10¹⁷ N·m | 0.184 km³ | 2.73 cm | 9.15 cm | 0.12 |
| ≤ 1.0 km | 14 | 0.39 m | 4.1×10⁶ m³ | 1.2×10¹⁷ N·m | 0.183 km³ | 2.69 cm | 8.97 cm | 0.35 |
| ≤ 1.5 km | 25 | 0.54 m | 7.8×10⁶ m³ | 2.3×10¹⁷ N·m | 0.181 km³ | 2.64 cm | 8.77 cm | 0.60 |
| ≤ 2.0 km | 25 | 0.54 m | 10.4×10⁶ m³ | 3.1×10¹⁷ N·m | 0.180 km³ | 2.60 cm | 8.63 cm | 0.78 |
| ≤ 2.5 km | 25 | 0.56 m | 12.6×10⁶ m³ | 3.8×10¹⁷ N·m | 0.178 km³ | 2.58 cm | 8.54 cm | 0.89 |
| ≤ 3.0 km | 33 | 0.69 m | 15.9×10⁶ m³ | 4.8×10¹⁷ N·m | 0.177 km³ | 2.55 cm | 8.44 cm | 1.01 |
| ≤ 3.5 km | 33 | 0.65 m | 16.7×10⁶ m³ | 5.0×10¹⁷ N·m | 0.176 km³ | 2.55 cm | 8.44 cm | 1.01 |
| trace-limited (≤ 4 km) | 33 | 0.61 m | 17.1×10⁶ m³ | 5.1×10¹⁷ N·m | 0.176 km³ | 2.55 cm | 8.45 cm | 1.00 |

## Depth reach of the uncapped solution

- cumulative potency fraction above 1 km: 0.27; above 2 km: 0.57; above 3 km: 0.83 (piecewise-linear in depth within patches)
- shallowest tested cap retaining ≥ 90% of the fault term's graben-band improvement: 3.0 km, with potency 15.9×10⁶ m³ (μ·potency 4.8×10¹⁷ N·m, vs 5.1×10¹⁷ N·m uncapped)

## Scope

- The sweep truncates the depth extent of the fixed mapped-trace, fixed-dip fault surface and re-solves at the reused smoothing weight. The retained-improvement threshold is a comparison rule, not a confidence level or a data-required depth. The sweep does not test other fault geometries (listric surfaces, distributed fissuring, the unsupported west trace), and the slip values remain regularized magnitudes. Neither the capped nor the uncapped moment is a physical upper confidence bound.
- Truncation changes patch dimensions, Green's functions, and the retained smoothing operator. These are not nested fits on one fixed basis, so a capped fit can have lower RMS; some saved rows improve slightly on the uncapped fit. The dike-volume column records the change in the opening/slip partition.
- Lowering the uniform μ rescales the moment only if both Lamé constants are scaled together. A lower stiffness value is an assumed sensitivity, not a value inferred from these data. A depth-varying soft layer would also change the displacement calculation and has not been modelled.
- The graben band includes axial residuals and support where loop closure flags 1–4-cycle inconsistencies (about 12–48 cm LOS). The check cannot identify the faulty leg or show that all flagged axial support was excluded from the fits. Both the RMS levels and their differences are conditional on the retained observations and model assumptions.
- Both input chains run the same sweep (fault_depth_mapped.md / fault_depth_mapped_gslc.md).
