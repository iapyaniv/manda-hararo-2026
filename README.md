# Manda Hararo 2026

Code, results and inputs behind the analysis of the dyke that intruded the
Dabbahu–Manda Hararo rift segment, Afar, Ethiopia, on 16–17 July 2026: NISAR
L-band and Sentinel-1 C-band interferometry, a uniform and a distributed dyke
and fault model, the displacement decomposition, timing and post-intrusion
checks, and the figures of the paper "Dabbahu–Manda Hararo reawakening"
(Pai et al., in preparation for *Geophysical Research Letters*).

Archived at Zenodo: software DOI to be filled in; data DOI to be filled in.

Maintainer: Vinay Pai, Yowon Labs, <vinay@yowon.org>.

## Key results

| Quantity | Value | Where |
|---|---|---|
| Co-event line-of-sight maxima | +76.6 cm ascending, +58.2 cm descending | `results/los_summary.md` |
| Onset | after 15 July, present by 18 July 2026 | `results/los_summary.md`, `results/probe_history.md` |
| Preferred model (distributed opening with a candidate east graben fault, 30-m NISAR inputs) | 0.176 km³, peak opening 2.38 m, joint RMS 2.55 cm; fault slip up to 0.61 m, conditional on the assumed trace | `results/distributed_fit_mapped_gslc.md` |
| Uniform dyke | strike 153°, length 17.3 km, depth 1.3–5.9 km, opening 2.32 m, 0.187 km³ | `results/okada_fit.md`, `results/model_sensitivity.md` |
| Sentinel-1 cross-check of the NISAR model | r = 0.974 ascending, 0.933 descending | `results/s1_crosscheck.md`, `results/s1_tracks.md` |
| Post-event motion with the dyke pattern | 0.107 of the co-event pattern in the NISAR pair 19 July–12 August; 0.140 and 0.053 in the first Sentinel-1 windows; a descriptive 4.3-day decay | `results/continuation.md`, `results/s1_tracks.md` |
| Pre-event NISAR windows | sixteen windows, November 2025 to 15 July 2026, no signal above 3.4 cm | `results/probe_history.md` |
| 2005–2010 openings against 2026 | the 2026 dyke sits where the earlier sequence added the least opening | `results/opening_history.md`, `results/uperp_three_epoch.md` |

Each report in `results/` states its inputs, settings and limits. Read
[`LIMITATIONS.md`](LIMITATIONS.md) before using any of these numbers.

## Figures

![Figure 1 draft: the segment, the 2026 interferogram, the 2005 and 2008 dykes, and the 2014–2026 Sentinel-1 series](results/figure1_draft.png)

*Figure 1 (draft).* The segment and the July 2026 intrusion in NISAR and
Envisat interferograms, with the Sentinel-1 cumulative displacement of
2014–2026. `python -m mh26.figure1`.

![The four-geometry displacement decomposition in the dyke frame](results/summary_decomposition_dike_frame_pub.png)

*Displacement decomposition.* Opening-perpendicular, along-dyke and vertical
displacement from the two NISAR and two Sentinel-1 geometries.
`python -m mh26.summary_figures`.

![The joint dyke and fault model](results/summary_source_geometry.png)

*The model.* Distributed opening on the dyke plane and slip on the candidate
east graben fault, with the observed, modelled and residual fields.
`python -m mh26.summary_figures`.

![Opening of the September 2005 dyke, the 2006–2010 sequence and the July 2026 dyke on one section](results/opening_history_common_scale_gslc.png)

*Three epochs of opening.* September 2005, the cumulative 2006–2010 sequence
and July 2026 on the same section of the rift, with the depth-averaged
profiles. `python -m mh26.opening_history`.

More figures, mostly for the supplement: the rewrapped data–model–residual
panels (`results/summary_los_model_residual.png`), the loop-closure check
(`results/loop_closure_nisar.png`), the fault depth sweep
(`results/fault_depth_mapped_gslc.png`), the post-event decay fit
(`results/decay_fit_windows.png`), the phase-gradient fault evidence
(`results/summary_phase_gradient.png`), the stress change on the fault
(`results/summary_fault_stress.png`) and the along-rift flank averages
(`results/along_rift_halves_gslc.png`).

## Reproduce

```bash
python3 -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements-lock.txt
# unpack the data archives into data/ as described in data/README.md
python -m mh26.validate --profile full
```

The validation profile recomputes the committed results from the data
archives and compares them with the stored checksums. To rerun the analysis
step by step, follow the processing order in
[`ANALYSIS_GUIDE.md`](ANALYSIS_GUIDE.md); `requirements-lock.txt` is the
environment the committed outputs were made with.

## Data

- The data record holds the NISAR GUNW crops, the reprocessed 30-m NISAR
  interferograms and unwrapped fields, the loop-closure inputs, the analysis
  GeoTIFFs and the LiCSAR Sentinel-1 pairs, each archive with a SHA-256
  manifest. Layout, verification and rebuild instructions:
  [`data/README.md`](data/README.md).
- Contributed inputs, kept in the repository with their authors' credit: the
  2005–2010 dyke opening models and the Envisat track-300 interferograms of
  Hamling et al. (2009, 2010) (`data/ihamling/`); the 2005 dyke plane of
  Wright et al. (2006) (`data/complex_geom2_revised_test.dat`); the
  Sentinel-1 cumulative displacement maps of Pedro Espín Bedón and Ben
  Ireland (`data/sentinel1-cumulative/`); and the site series digitized from
  Ben Ireland's figure (`data/ireland-series/`).
- Sources: NISAR L2 GSLC and GUNW provisional products, NASA/JPL via ASF DAAC
  (10.5067/NIL2GSLC-P1, 10.5067/NIL1GUNW-P1); Sentinel-1 interferograms from
  COMET LiCSAR under the [portal's conditions](https://comet.nerc.ac.uk/COMET-LiCS-portal/);
  Copernicus GLO-30 terrain for the figures; Envisat ASAR data © ESA.

## Documents

- [`ANALYSIS_GUIDE.md`](ANALYSIS_GUIDE.md): what each step does, in processing order, and the code map.
- [`PROCESSING_NOTES.md`](PROCESSING_NOTES.md): the conventions and choices that affect the numbers.
- [`LIMITATIONS.md`](LIMITATIONS.md): what the results do and do not establish.
- [`data/README.md`](data/README.md): archives, checksums, contributed inputs.
- [`CHANGELOG.md`](CHANGELOG.md): releases, the commit history by day, and the former README's analysis notes.

## Licence and citation

The code is released under the MIT licence (`LICENSE`); the data record under
CC BY 4.0. The contributed inputs remain their authors' and are cited to
their papers. Cite the Zenodo records (DOIs above) and the paper.
