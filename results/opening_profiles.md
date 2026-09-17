# Depth-averaged opening per intrusion (after Hamling et al. 2009, Fig. 13)

NISAR chain of the 2026 fit: gunw (results/distributed_fit_geometry_4a7b2386c8b5.json).

Top-10-km depth-averaged opening of every intrusion of the 2005–2010 Dabbahu–Manda Hararo sequence and of the July 2026 dike, versus distance along the Wright et al. (2006) plane (NW→SE s axis of the 2026 contract), in the style of Fig. 13 of Hamling et al. (2009). The 2005 profile averages Ian Hamling's 10-km September 2005 model, the 2006–2010 profiles the upper 10 km of his canonical thirteen files, and the 2026 profile this repository's fixed-plane fit (its plane spans exactly 0–10 km). On the native grid a profile integrates exactly to the epoch's above-10-km volume: average × 10 km × column width = column volume (cross-checked against `uperp_three_epoch*.json` on every run). The drawn curves interpolate linearly between column centres on the projected axis and are not volume-preserving: the largest drawn-area deviation is February 2009 at 5.7% and the cumulative stack draws as 3.20 km3 against the native 3.19 km3 (every drawn area is in the JSON's display block), so quoted volumes always come from the native-grid tables. The delivered inversions carry no uncertainties, so unlike the published figure the panels have no error bars. Magnitude and extent comparisons also inherit a threshold asymmetry, stated on the cumulative panels: every canonical post-2005 file floors small opening (0.31–0.72 m, consistent with the papers' discard of opening below the local 2σ uncertainty), while the September 2005 model and the 2026 fit keep all fitted opening. As one uniform illustrative floor — a sensitivity, not a reconstruction of the historical thresholds — 0.5 m reduces the July 2026 volume from 0.182 to 0.156 km3 (about 14%); the matching cumulative sensitivity is in `uperp_three_epoch*.md`, whose remaining comparability caveats (different regularization and greens functions) apply unchanged.

| intrusion | file | volume above 10 km (km3) | max average opening (m) | at s (km) |
|---|---|---|---|---|
| September 2005 | model_slip_sept05.txt | 2.230 | 10.29 | 3.6 |
| June 2006 | dyke_only_june06_ex.txt | 0.098 | 1.47 | 38.0 |
| July 2006 | dyke_only_july06_ex.txt | 0.025 | 0.54 | 44.5 |
| September 2006 | dyke_only_sept06_ex.txt | 0.044 | 0.69 | 49.2 |
| December 2006 | dyke_only_dec06_ex.txt | 0.044 | 0.79 | 43.6 |
| January 2007 | dyke_only_jan07_ex.txt | 0.013 | 0.31 | 47.3 |
| August 2007 | dyke_only_aug07_ex.txt | 0.040 | 0.74 | 44.5 |
| November 2007 | dyke_only_nov07_ex.txt | 0.111 | 1.01 | 53.9 |
| April 2008 | dyke_only_apr08_ex.txt | 0.067 | 0.70 | 57.7 |
| July 2008 | dyke_only_july08_ex.txt | 0.021 | 0.39 | 51.1 |
| October 2008 | dyke_only_oct08_ex.txt | 0.182 | 2.08 | 33.0 |
| February 2009 | dyke_only_feb09_ex.txt | 0.069 | 1.15 | 39.6 |
| June 2009 | dyke_only_july09_ex.txt | 0.035 | 0.62 | 44.5 |
| May 2010 | dyke_only_may10_ex.txt | 0.035 | 0.70 | 40.6 |
| July 2026 (gunw) | distributed_fit_geometry_4a7b2386c8b5.json | 0.182 | 1.30 | 56.0 |

The cumulative stack (bottom panel of `opening_profiles.png`, also standalone as `opening_profiles_cumulative.png`) sums the fifteen profiles in event order: 3.19 km3 above 10 km in total. Panel labels are the intrusion months of the Hamling et al. (2009, 2010) sequence — the july09 file is the re-estimate of the June 2009 intrusion — and exact onset dates for captions belong to the source catalogue. Beyond the SE end of the 2026 model plane (71.0 km) there is no 2026 model, and that region is shaded grey in the July 2026 panel and in both cumulative panels (the stack there holds historical models only). The dashed line is the Wright et al. (2006) model-plane SE edge: no 2005 model column southeast of it holds opening, and the drawn 2005 taper crosses it only by the interpolation span between column centres (0.10 m at the line, positive to 61.4 km).
