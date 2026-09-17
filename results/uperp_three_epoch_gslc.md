# Three-epoch opening comparison (Hamling dike files)

NISAR chain of the 2026 fit: gslc (results/distributed_fit_geometry_4a7b2386c8b5_gslc.json).

Openings of the September 2005 megadike and the fifteen post-2005 opening files are the inversions provided by Ian Hamling (GNS) on the ~74 km Hamling computational grid (75 columns x 20 one-km depth rows, UTM 37N). The 2026 distribution is this repository's fixed-plane fit on the Wright et al. (2006) plane. The Hamling columns project monotonically into that plane's s axis (median off-plane offset 93 m, maximum 1412 m mid-plane), so `uperp_three_epoch_gslc.png` compares the three epochs by along-axis position on one axis.

| epoch | volume (km3) | above 10 km (km3) | max opening (m) |
|---|---|---|---|
| sept05 | 3.554 | 1.815 | 12.92 |
| sept05 — updated 10-km model (drawn in the figure) | 2.230 | 2.230 | 16.09 |
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
| 2006-2010 cumulative (13 files, primary) | 0.871 | 0.783 | 6.67 |
| — alternative: all 15 delivered files | 0.968 | 0.879 | 8.44 |
| — alternative: superseded primary (kept june09 + apr10) | 0.898 | 0.809 | 6.67 |
| — alternative: untrimmed apr10 in place of may10 | 0.893 | 0.804 | 6.67 |
| July 2026 (gslc) | 0.180 | 0.180 | 2.31 |

The delivery has 15 post-2005 files against the 13 published post-megadike intrusions. Ian Hamling's update of 17 September 2026 lists the canonical 13 files and calls june09 and apr10 defunct earlier versions, so the primary cumulative keeps july09 (the re-estimate of the June 2009 event) and may10 (apr10 with every cell below 0.5 m set to zero, 38% of that field's volume — kept for consistency with the other floored files). The data-derived pair diagnostics match that identification: july09's positive cells all lie inside june09's at correlation 0.994 while no other pair exceeds |r| = 0.44. The alternatives above record the superseded selection and the untrimmed-apr10 variant.

The row order of the opening files follows plot_3d_dyke_only.dat; Ian Hamling's own plotting script (doplot3d.m, provided 17 September 2026, tracked in data/ihamling/2026-09-17/) indexes an opening file by the geometry file's row number, confirming that reading.

Comparability caveats: every post-2005 file carries a positive floor of 0.29-0.72 m, consistent with the papers' discard of opening below the local 2-sigma uncertainty (typically 0.3-0.4 m), while the delivered sept05 file is the raw inversion regridded to 20 km (Ian Hamling, 17 September 2026): it is unfiltered (cells down to 4e-4 m), integrates to 3.55 km3 against the published ~2.5 km3, and carries 1.74 km3 of ill-constrained opening below 10 km (deep lobes peaking at 13-20 km, up to 12.9 m). The figure therefore shows the upper 10 km of every panel, and its September 2005 panel draws his updated 10-km model (`model_slip_sept05.txt`, provided 17 September 2026), decoded per his layout note of the same day: the dyke section is lines 171-920 of the 931 values — 10 one-km depth rows per column, depth-fastest, at the same 75 horizontal positions — and the trimmed ends are fault sections. The decoded field correlates at r = 0.69 with the raw field's upper 10 km (two inversions of the same event) and integrates to 2.230 km3 with the geometry file's column widths (maximum 16.09 m, all above 10 km), near the published ~2.5 km3 — the earlier whole-file total of 2.478 km3 at 1-km2 patches included the fault values now trimmed. The updated model holds 0.000 km3 southeast of the Wright 2006 plane edge (the raw field holds 0.22 km3 there, 0.04 km3 of it above 10 km), so the dashed edge line also bounds the drawn 2005 opening. Like-for-like volume comparisons should quote the above-10-km column. The 2026 totals keep all fitted opening. As a uniform-floor sensitivity — not a matched reconstruction — a 0.5-m floor gives 0.817 km3 for the primary cumulative and 0.154 km3 for 2026. The 2026 model plane ends at 10 km. See also the notes block in the JSON and LIMITATIONS.md. Exact event dates are not part of the delivery; take them from the Hamling et al. (2009, 2010) catalogue when captioning.
