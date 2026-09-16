# Three-epoch opening comparison (Hamling dike files)

NISAR chain of the 2026 fit: gslc (results/distributed_fit_geometry_4a7b2386c8b5_gslc.json).

Openings of the September 2005 megadike and the fifteen post-2005 opening files are the published inversions provided by Ian Hamling (GNS) on the ~74 km Hamling computational grid (75 columns x 20 one-km depth rows, UTM 37N). The 2026 distribution is this repository's fixed-plane fit on the Wright et al. (2006) plane. The Hamling columns project monotonically into that plane's s axis (median off-plane offset 93 m, maximum 1412 m mid-plane), so `uperp_three_epoch_gslc.png` compares the three epochs by along-axis position on one axis.

| epoch | volume (km3) | above 10 km (km3) | max opening (m) |
|---|---|---|---|
| sept05 | 3.554 | 1.815 | 12.91 |
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
| 2006-2010 cumulative (14 files, duplicate-free) | 0.933 | 0.845 | 7.72 |
| — alternative: all 15 delivered files | 0.968 | 0.879 | 8.44 |
| — alternative: catalogue 13 events | 0.877 | 0.788 | 6.67 |
| July 2026 (gslc) | 0.180 | 0.180 | 2.31 |

The delivery has 15 post-2005 files against the 13 published post-megadike intrusions. Exact duplication: may10 is an exact sub-array of apr10 (36 nonzero cells). The cumulative drops the contained file(s) so no delivered value is counted twice; the file-to-event mapping (apr10, july09) awaits Ian Hamling's confirmation, and the table quotes the all-files and catalogue-count sums as alternatives.

The 2005 model carries 1.74 km3 below 10 km depth (a deep lobe peaking at 13-20 km in its northern half); the 2026 model plane ends at 10 km, so like-for-like volume comparisons should quote the above-10-km column. Caveats: see the notes block in the JSON and LIMITATIONS.md. Exact event dates are not part of the delivery; take them from the Hamling et al. (2009, 2010) catalogue when captioning.
