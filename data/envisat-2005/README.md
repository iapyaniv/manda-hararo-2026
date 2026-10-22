# The September 2005 dyke in Envisat track 49 (processed here)

`t049_20050506_20051028.npz` is the wrapped interferogram of Figure 1
panel (b): the descending Envisat track-49 pair of 6 May and 28 October
2005 (frame 3357, swath IS2, V/V), the pair shown in Fig. 2a of Wright et
al. (2006), processed from the ESA archive products for this study.

## Inputs

- Scenes: the ESA online dissemination delivers slices generated on demand,
  named differently from the catalogue entries. Downloaded on 22 October
  2026 with an ESA EO Sign In (Envisat ASAR data © ESA; the products are
  not redistributed here):

  | Delivered product | Bytes | SHA-256 | Catalogue identifier |
  |---|---|---|---|
  | `ASA_IMS_1PNESA20050506_071515_000000182037_00049_16636_0000.N1` | 626,464,409 | `d5ae7200a44323e2736adbed39f544d0ac6a02a18ea127d247b130dbd9a98d53` | `ASA_IMS_1PNESA20050506_071516_000000152037_00049_16636_0000` |
  | `ASA_IMS_1PNESA20051028_071513_000000182042_00049_19141_0000.N1` | 625,989,850 | `18c49b5fa6c650aa6608c2fbcd23e4fb08b5259c38aab6b0543aa1ca7f634330` | `ASA_IMS_1PNESA20051028_071514_000000152042_00049_19141_0000` |

  Absolute orbits 16636 and 19141; sensing 07:15:14–07:15:34 UTC and
  07:15:12–07:15:32 UTC; both carry `PRODUCT_ERR=1` with the processor
  warning flags common to the archive's IMS products.
  The baseline between the two orbits is about 250 m (GMTSAR
  `baseline_start` 251 m).
- Orbits: DORIS precise orbits `DOR_VOR_AXVF-P20120424_180800_20050505_215528_20050507_002328`
  and `DOR_VOR_AXVF-P20120424_195100_20051027_215528_20051029_002328` from
  the public ESA archive (`https://step.esa.int/auxdata/orbits/Doris/vor/2005.zip`).
- Terrain: SRTMGL3 at 3 arcseconds (GMT remote `@earth_relief_03s`) with the
  EGM2008 geoid (`@earth_geoid_01m`) added for ellipsoidal heights,
  39.5–41.7 E, 11.2–13.4 N.

## Processing

GMTSAR 6.2.0 (source commit `527dce7d`), `p2p_processing.csh ENVI_SLC`
with the default configuration (`pop_config.csh ENVI_SLC`) except
`mask_water = 0`: the products decoded with the orbits above, aligned by
cross-correlation, the topographic phase removed with the terrain grid,
a 200-m Gaussian filter, decimation 2, no unwrapping, geocoded with
GMTSAR's coherence mask at 0.10; then the cosine and sine of the
filtered phase (`phasefilt.grd`, Gaussian plus Goldstein filtering)
geocoded separately with `proj_ra2ll.csh`, because GMTSAR projects grids
through a block median, which on wrapped radians mixes values across the
±π boundary (its `phasefilt_ll.grd` is wrong by up to half a cycle at
such nodes). `python -m mh26.envisat_2005 --gmtsar <work dir>` then takes
the phase from the geocoded cosine and sine, keeps nodes with coherence of 0.15 or
more, fixes the sign by the step-agreement
test against the forward model of I. Hamling's revised September 2005
opening in the pair's own look geometry, and averages the unit phasors
onto the product grid. The report is `results/envisat_2005.md`.

## Product

- `range_change_cm`: the range change wrapped at one fringe (2.8118 cm,
  lambda/2), increasing with motion away from the satellite, the
  convention of the Hamling track-300 grids (`data/ihamling/`); NaN where
  blank.
- `coherence`: the processor's coherence averaged onto the same grid.
- `lon`, `lat`: node centres of the 0.0005-degree grid (40.20–41.14 E,
  11.84–12.92 N).
- `heading_deg`, `los_azimuth_deg`, `incidence_deg`, `unit_east_north_up`:
  the look geometry at the scene centre from the reference orbit.
- `scenes`, `catalogue_scenes`, `dates`, `track`, `sense`, `fringe_cm`,
  `coherence_min`: provenance.
