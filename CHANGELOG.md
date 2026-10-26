# Changelog

The dated record of this repository: archives and releases first, then the
changes by day (those that altered data, processing, results, figure content
or the documents; the git log has the rest), then the analysis notes that the
README carried until 22 October 2026, kept word for word.

## Archives and releases

- 2026-10-26: Zenodo records, version 1.0.0, software (DOI to be filled in)
  and data (DOI to be filled in), built with `scripts/build_release.py` from
  the commit tagged `v1.0.0`; checksums in `dist/SHA256SUMS.txt`. The data
  record has the four analysis archives, the closure add-on and the new
  figure-input add-on `manda-hararo-2026_data-figures.zip` (the Hamling
  models and track-300 interferograms, the 2005 dike plane, the Sentinel-1
  cumulative maps and digitized series, the Envisat track-49 interferogram
  of 2005 and the NISAR display strip north of the analysis crop).
- 2026-09-14: the closure-data add-on `manda-hararo-2026_data-closure.zip`
  prepared locally (SHA-256
  `a13584043d975a616fd4a22593c4aee8c3d0b7ac9319d5f3ad1592cc25923b0a`).

`scripts/build_release.py` builds the archives deterministically: entry
timestamps and the archive comment come from the build commit, so a rebuild
at a new commit changes every archive hash even when the packed data is
unchanged.

- Release `analysis-2026-08-18` (commit `f9799dd`): code, core, and
  pre-event archives with their checksum file.
- Release `analysis-2026-08-20`: revised code archive, the two-track
  `s1tracks` add-on, and the `gslc` add-on with the 30-m reprocessed
  interferograms and unwrapped fields. The 2026-08-18 **core** and
  **pre-event** archives are unchanged, are not re-uploaded, and keep their
  published checksums:

      b245d47e063eb278cef203e4305a93d54ed72f721fda1e834f67eb2a51861e1b  manda-hararo-2026_data-core.zip
      f6dc9271b5aaae87300687ad7d729fde33d13bc1a127eadf7a7e350ba6a45a8f  manda-hararo-2026_data-preevent.zip

From `analysis-2026-08-20` on, each release ships a
`SHA256SUMS-release-<tag>.txt` manifest listing the exact hash of every
asset the release refers to, including retained assets from earlier
releases. A locally generated `dist/SHA256SUMS.txt` records locally built
archives; its core/pre-event hashes can differ from retained release assets
built at an earlier commit. It is not the published release manifest.

A separate closure-data add-on was prepared locally on 14 September 2026.
It contains the 18 exact interferogram/cache files, an input checksum
manifest, and provenance metadata. It has not been uploaded. The archive
uses the source-base timestamp and records the modified working tree and
exact inventory hash in its ZIP comment. Build this add-on alone with
`python scripts/build_release.py --closure-only`; see `data/README.md` for
installation and read-only reproduction. The builder validates the saved
source, settings, reports, and all 18 inputs before writing. The prepared
archive SHA-256 is:

```text
a13584043d975a616fd4a22593c4aee8c3d0b7ac9319d5f3ad1592cc25923b0a  manda-hararo-2026_data-closure.zip
```

## Changes by day

### 2026-08-18
- Manda Hararo 2026 analysis

### 2026-08-24
- Add two-track Sentinel-1 analysis
- Add GSLC interferogram reprocessing
- Add distributed opening and fault inversion
- Add joint displacement decomposition
- Add fault-trace geometry tests
- Add distributed-source stress analysis
- Add projected-noise uncertainty analysis
- Add rewrapped fringe model comparison
- Strengthen scientific result validation
- Add data archives for new analysis products
- Add external geometry and four-look fit products
- Map stress on saved receiver planes
- Record GSLC unwrap cost sensitivity
- Document method readiness and open inputs
- Add Xu et al. (2026)-style summary figures

### 2026-08-26
- Fit the 2026 intrusion on the Wright et al. (2006) dike plane
- Harden the Wright-2006 contract against frame and input misuse

### 2026-08-27
- Refine the 2026 Wright-plane opening cross-section

### 2026-08-28
- Sentinel-1 cumulative displacement files, with their README

### 2026-08-29
- Rotate the free three-component decomposition into the dike frame

### 2026-08-31
- Render the fixed 2005-plane fit at publication grade

### 2026-09-01
- Apply the observation reference operator to the model in both inversions and in the diagnostics; every fit, figure and report regenerated from the corrected estimator
- Rebuild the Wright-2006 contract from the current mapped fit

### 2026-09-02
- Flip the dike-frame U-parallel sign to near-north and add the LOS inputs to its figure
- Add a three-epoch loop-closure check of the reprocessed NISAR unwrapping

### 2026-09-03
- Sentinel-1 cumulative displacement maps (LiCSBAS, frames 014A_07688_131313 and 079D_07694_131313, October 2014 to July and to August 2026, with and without the ERA5 correction) and their README, from Pedro Espín Bedón

### 2026-09-09
- Accept the loop-closure pre epochs in the GSLC granule manifest check
- Add a depth-extent sweep of the mapped graben-fault term

### 2026-09-16
- Add the three-epoch opening comparison from the Hamling dike files
- Resolve the duplicate files in the Hamling delivery by evidence, de-duplicate the cumulative, and qualify the volume comparison
- Recompute the summed statistics and the bounds in validation by full recomputation
- Question the 2005 model version and derive fit stats from the vector

### 2026-09-17
- Adopt Hamling's canonical event mapping and cut the display at 10 km
- Decode the updated September 2005 model and draw it on the section
- Add the depth-averaged opening profiles in the Hamling 2009 Fig. 13 style
- Record the drawn areas of the display resampling and mark missing 2026 coverage
- Put the floored cumulative on the report's 0-10 km interval; 2005 model as supplied
- Export the 2026 opening model on the Hamling computational grid

### 2026-09-30
- Combine the epoch sections and cumulative stack into one opening-history figure

### 2026-10-01
- Add a manuscript variant of the dike-frame figure

### 2026-10-07
- Add an AGU word-count estimator
- Add north/south-half along-rift averages on the flanks
- Reference the model like the data and scope the half-flank errors
- Correct overstated claims in the along-rift reports
- Draw the decay fit as a one-panel supplement figure

### 2026-10-08
- Decay figure: fitted values as dashed bars over their windows, with intervals, track identity and role-based selection

### 2026-10-22
- Move the Sentinel-1 cumulative displacement maps under data/
- Add the draft of Figure 1 for the GRL paper: the segment and the July 2026 NISAR interferogram, the September 2005 dyke from the Envisat track-49 pair processed here, the October 2008 dyke in Envisat track 300, and the 2014–2026 Sentinel-1 cumulative displacement wrapped at one fringe, with a display strip north of the NISAR analysis crop
- Figure 4: the manuscript layout of the opening history, its caption written from the report numbers
- Add the MIT licence for the code
- Reorganize the documentation: a short README with the figures, the history in a changelog
- Prefer an Earthdata user token for the download paths

### 2026-10-23
- Processing notes: the baseline and geolocation checks of the 2005 Envisat pair

### 2026-10-25
- Figure 1: the manuscript version without the draft notes (results/figure1.png and .pdf; `--notes` writes the annotated draft); the IGSSA location only, the ISC location dropped
- Figure 4: the opening values labelled on the contours of the manuscript sections

### 2026-10-26
- Release build: the figure-input add-on and its manifest `data/SHA256SUMS-figures.txt`

## Analysis notes kept from the README (as of 22 October 2026)

### Results

- Co-event LOS maxima are +76.6 cm ascending and +58.2 cm descending.
  See `results/los_summary.md`.
- The onset is after 15 July and present by 18 July, based on the NISAR
  control and LiCSAR event pairs.
- The uniform dike has strike 153.0°, length 17.31 km, depth 1.30–5.94 km,
  opening 2.32 m, and RMS 2.91 cm. See `results/okada_fit.md`.
- Processing cases span 17.3–17.7 km length, 1.1–6.9 km depth, and
  2.05–2.40 m opening. Geometric volume is 0.187 km³, with a
  processing-choice range of 0.186–0.208 km³. See
  `results/model_sensitivity.md`.
- The Sentinel-1 ascending comparison gives r = 0.974, amplitude ratio 0.92,
  and RMS 1.8 cm. The descending comparison gives r = 0.933, amplitude ratio
  0.88, and RMS 2.7 cm. The model passes through the same far-field plane
  projection as the observations in these comparisons. See `results/s1_crosscheck.md` and
  `results/s1_tracks.md`.
- The NISAR post-event projection is +0.107 for 19 July–12 August, with a
  processing-choice range of +0.100–+0.121. See `results/continuation.md`.
- Sentinel-1 post-event projections are +0.140 ascending for 18–30 July and
  +0.053 descending for 23 July–4 August. Later windows are at the null
  level. One descriptive fit across all six windows gives a total near 0.15
  of the co-event pattern and a 4.3-day decay time from 18 July. It has no
  formal uncertainty. See `results/s1_tracks.md`.
- The illustrative Sentinel-1-only East/Up decomposition gives +18 cm flank
  divergence. Its two windows are offset by five days. See
  `LIMITATIONS.md` and `results/s1_tracks.md`.

### Model extensions and checks

- The preferred fit uses the reprocessed 30-m GSLC inputs and a candidate
  eastern fault trace. It gives 0.176 km³ dike volume, 2.38 m peak opening,
  and 2.55 cm joint RMS. The 0.61 m maximum fault slip is conditional on
  the assumed trace, dip, depth extent, and smoothing. See
  [`results/distributed_fit_mapped_gslc.md`](results/distributed_fit_mapped_gslc.md).
- With the standard 80-m GUNW inputs, the trace-supported fit includes only
  the east fault. It reduces joint RMS from 2.37 to 2.15 cm and improves both
  5-km holdout splits. The 2.08 cm two-fault idealization remains a geometry sensitivity.
  The held-out-sensor transfer gain belongs to the idealized two-fault
  geometry; the mapped east-fault term is constrained mainly by NISAR. See
  `results/distributed_fit_mapped.md` and `LIMITATIONS.md`.
- The joint decomposition, free three-component variant, closure diagnostics,
  and 30-m-input repeat are in `results/decomposition_joint*.md`. Counts now
  distinguish exactly three geometries from four geometries.
- The free three-component solve is also rotated into the dike frame: U⊥
  perpendicular to the fitted strike (ENE positive, the same axis as the
  constrained solve), U∥ along the dike axis (333° NNW positive), and up. The rotation
  is exact and the rotated formal sigmas carry the east–north
  cross-covariance. U∥ is the weakly observed component: with the 30-m
  inputs its flank probes are −5.3 cm (east) and +4.2 cm (west) against a
  median formal σ of 4.1 cm; the GUNW inputs give −10.0 and +3.4 cm, and
  replacing the fitted strike by the southernmost 2005 segment strike (4.5°
  different) moves them by about 5 cm. U∥ is therefore reported but not
  resolved. See the dike-frame section of `results/decomposition_joint*.md`.
- The GSLC reprocessing comparison and its near-axis differences are in
  `results/gslc_gunw_compare_*.json`.
- A single-tile SMOOTH/MCF sensitivity changes 6.1% of near-axis pixels by
  more than half a fringe in the ascending model-input pair. It agrees
  closely off axis. This is not an independent unwrap validation. See
  `results/gslc_unwrap_sensitivity_172_008_A.json`.
- Phase-gradient ridges, threshold sensitivity, and an along-strike alignment
  null are in `results/fault_traces.json`. Only the east ridge has consistent
  LOS-step support. The default trace-column fit excludes the unsupported west
  ridge.
- The trace-column fit is a geometry sensitivity, not a fault map. Its source
  path and SHA-256 hash are stored in `results/distributed_fit_mapped.json`.
- The stress report uses the source and receiver geometry stored in that fit.
  It includes spatial stress maps and a depth summary. It does not use the
  uniform rectangle's fitted top. See `results/stress_change.md`.
- The Monte Carlo applies the far-field datum and ramp projection to the model
  and to every covariance draw. See `results/distributed_uncertainty.md`.
- Same-chain rewrapped observations, model fields, and residuals for all four
  geometries are in `results/distributed_fit_fringe_panels_mapped_gslc.png`.
  The separate GUNW-fit-to-GSLC check is in
  `results/fringes_model_comparison.png`.
- Eight draft publication-style figures summarize the trace-supported fit,
  the free E–N–U decomposition, the same solve in the dike frame (U⊥, U∥,
  up), the physical source, the phase-gradient faulting evidence, the
  dike-induced stress on the east-fault receivers, and the fit on the
  fixed 2005 plane (map view and U⊥ cross-section): `results/summary_*.png`
  with matching vector PDF files and the
  identity manifest `results/summary_figures.json`
  (`mh26.summary_figures`; needs the terrain tile described in
  `data/README.md`).
- `mh26.distributed --geometry FILE` accepts a complete saved geometry,
  validates its patch contract, records file and geometry hashes, and uses
  hash-suffixed outputs.
- `mh26.wright2006` builds such a geometry from the 2005 dike-plane source
  model of Wright et al. (2006) in `data/complex_geom2_revised_test.dat`
  (five chained vertical segments, 61 km, 0–10 km depth; see
  `data/README.md`). The 2026 uniform-fit dike is aligned with the
  southernmost segment (strike differs by ~4.5°; endpoints 0.45–0.92 km off
  its line) and overruns its SE end by ~2.8 km, so the plane is extended 10 km
  SE. The 2026 joint inversion on that fixed plane gives joint RMS 2.20 cm
  (trace-supported free-trace fit: 2.15 cm), maximum opening 2.30 m, and
  volume 0.182 km³ (`results/distributed_fit_geometry_f9c903279610.md`);
  on the reprocessed 30-m (GSLC) inputs the same comparison is 2.61 cm
  against 2.55 cm, with maximum opening 2.31 m and volume 0.180 km³
  (`results/distributed_fit_geometry_f9c903279610_gslc.md`).
  `results/uperp_2026_wright2006.png` is one 2026 U⊥ cross-section on
  the Wright plane, from the GUNW-input fit (a `--gslc` fit writes a
  `_gslc` name). It uses the along-strike/depth view of an individual
  opening panel in Hamling et al. (2010) Figure 2; it is not a reproduction
  of that multi-event opening-and-stress figure. Its fixed 0–2.5 m scale can
  support a later comparison after the historical coordinate definitions
  are checked. The panel labels the 10-km southeast extension and the
  0–10 km model depth. The historical opening distributions are still
  awaited. This is a fixed-historical-plane variant with nominally
  transferred smoothing, not a co-planarity result; see `LIMITATIONS.md`.
  The publication-style versions of the fixed-plane results are
  `results/summary_los_model_residual_wright2006.png` (observed, model,
  residual, and profile, in the layout of the free-trace figure) and
  `results/summary_uperp_wright2006.png` (the cross-section restyled),
  both from the 30-m-input fit.

Sixteen pre-event windows from November 2025 through 15 July 2026 have a
maximum absolute 1-km probe mean of 3.4 cm (`results/probe_history.md`). Nine
pre-event Sentinel-1 pairs bound the post-event projection null at ±0.022
(`results/s1_tracks.md`). The post-event signal appears in both Sentinel-1
geometries and in NISAR, and it decays; no added-volume estimate is made.

The separate visual review of 24 additional January–June NISAR previews is
in `results/preevent_preview_review.md`; its input list is
`granules/preevent_preview_pairs.tsv`. It does not extend the quantitative
16-window result.
