# Graben-fault depth-extent sweep (mapped trace; GUNW product chain)

- reference fit: distributed_fit_mapped.json (λ = 17.78 reused; SHA-256 b67920d372db…); the uncapped baseline must reproduce it exactly before the sweep runs
- graben band: within 6 km of the trace, 20 km along strike (991 blocks); the fault term acts here, and the all-block RMS dilutes it
- fault moment = μ · Σ area·slip at μ = 30 GPa. Okada surface displacements depend on (λ, μ) only through α = (λ+μ)/(λ+2μ), so a uniform rescaling of both changes no prediction and rescales the moment linearly: the data do not constrain μ. A depth-varying μ is outside the half-space model.

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
- shallowest cap retaining ≥ 90% of the fault term's graben-band improvement: 2.5 km, with potency 13.3×10⁶ m³ (μ·potency 4.0×10¹⁷ N·m, vs 5.4×10¹⁷ N·m uncapped)

## Scope

- The sweep truncates the depth extent of the fixed mapped-trace, fixed-dip fault surface and re-solves at the reused smoothing weight. It does not test other fault geometries (listric surfaces, distributed fissuring, the unsupported west trace), and the slip values remain regularized magnitudes, not resolved slip.
- Truncation removes parameters, so a capped fit can only misfit more; the question is how little the graben band degrades while the potency (and so the moment) drops. The dike-volume column shows how much removed fault signal the opening absorbs.
- The band's RMS level is dominated by axial signal no model term fits (graben-floor cells, where the loop-closure check also flags a one-fringe unwrap risk), so differences between rows are meaningful; the level is not.
- Both input chains run the same sweep (fault_depth_mapped.md / fault_depth_mapped_gslc.md).
