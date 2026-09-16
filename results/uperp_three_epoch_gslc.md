# Three-epoch opening comparison (Hamling dike files)

NISAR chain of the 2026 fit: gslc (results/distributed_fit_geometry_4a7b2386c8b5_gslc.json).

Openings of the September 2005 megadike and the fifteen post-2005 opening files are the inversions provided by Ian Hamling (GNS) on the ~74 km Hamling computational grid (75 columns x 20 one-km depth rows, UTM 37N). The 2026 distribution is this repository's fixed-plane fit on the Wright et al. (2006) plane. The Hamling columns project monotonically into that plane's s axis (median off-plane offset 93 m, maximum 1412 m mid-plane), so `uperp_three_epoch_gslc.png` compares the three epochs by along-axis position on one axis.

| epoch | volume (km3) | above 10 km (km3) | max opening (m) |
|---|---|---|---|
| sept05 | 3.554 | 1.815 | 12.92 |
| june06 | 0.102 | 0.098 | 2.16 |
| july06 | 0.025 | 0.025 | 1.15 |
| sept06 | 0.050 | 0.044 | 1.03 |
| dec06 | 0.049 | 0.044 | 1.81 |
| jan07 | 0.013 | 0.013 | 0.70 |
| aug07 | 0.040 | 0.040 | 2.46 |
| nov07 | 0.116 | 0.111 | 2.03 |
| apr08 | 0.118 | 0.067 | 1.27 |
| july08 | 0.021 | 0.021 | 1.26 |
| oct08 | 0.192 | 0.182 | 2.47 |
| feb09 | 0.076 | 0.069 | 1.57 |
| june09 | 0.040 | 0.040 | 1.67 |
| july09 | 0.035 | 0.035 | 1.64 |
| apr10 | 0.056 | 0.056 | 1.67 |
| may10 | 0.035 | 0.035 | 1.67 |
| 2006-2010 cumulative (13 files, primary) | 0.898 | 0.809 | 6.67 |
| — alternative: all 15 delivered files | 0.968 | 0.879 | 8.44 |
| — alternative: 14 files (only may10 dropped) | 0.933 | 0.845 | 7.72 |
| — alternative: keep july09, drop june09 | 0.893 | 0.804 | 6.67 |
| — alternative: keep july09 + may10, drop june09 + apr10 | 0.871 | 0.783 | 6.67 |
| — alternative: catalogue 13 events (filtered may10) | 0.877 | 0.788 | 6.67 |
| July 2026 (gslc) | 0.180 | 0.180 | 2.31 |

The delivery has 15 post-2005 files against the 13 published post-megadike intrusions, and the data flag two as duplicate content: may10 equals apr10 with every cell below 0.5 m set to zero (38% of that field's volume), and july09's positive cells all lie inside june09's at correlation 0.994 while no other pair exceeds |r| = 0.44 — two estimates of one event. The primary cumulative therefore sums 13 files (drops may10 and july09); the alternatives above bracket the mapping ambiguity (0.871-0.968 km3), which Ian Hamling has been asked to resolve.

Comparability caveats: whether each file is the published model version is unconfirmed (asked of Ian Hamling). Every post-2005 file carries a positive floor of 0.29-0.72 m, consistent with the papers' discard of opening below the local 2-sigma uncertainty (typically 0.3-0.4 m), while sept05 is unfiltered (cells down to 4e-4 m), integrates to 3.55 km3 against the published ~2.5 km3, and carries 1.74 km3 below 10 km depth (deep lobes peaking at 13-20 km, up to 12.9 m) that the published opening figure does not show and that the papers attribute largely to the inversion fitting long-wavelength signals unrelated to dike opening. The 2026 totals keep all fitted opening. As a uniform-floor sensitivity — not a matched reconstruction — a 0.5-m floor gives 0.821 km3 for the primary cumulative and 0.154 km3 for 2026. The 2026 model plane ends at 10 km, so like-for-like volume comparisons should quote the above-10-km column. See also the notes block in the JSON and LIMITATIONS.md. Exact event dates are not part of the delivery; take them from the Hamling et al. (2009, 2010) catalogue when captioning.
