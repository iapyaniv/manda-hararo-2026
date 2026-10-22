# The September 2005 dyke in Envisat track 49 (this processing)

Pair 20050506_20051028 (descending track 49, frame 3357, IS2): `ASA_IMS_1PNESA20050506_071515_000000182037_00049_16636_0000` and `ASA_IMS_1PNESA20051028_071513_000000182042_00049_19141_0000` (ESA's on-demand slices of the catalogue scenes `ASA_IMS_1PNESA20050506_071516_000000152037_00049_16636_0000` and `ASA_IMS_1PNESA20051028_071514_000000152042_00049_19141_0000`), the pair of Wright et al. (2006) Fig. 2a. Processed with GMTSAR 6.2 (`p2p_processing.csh ENVI_SLC`, DORIS precise orbits, SRTM 3-arcsecond terrain on the ellipsoid, 200-m filter, Goldstein filter, decimation 2, no unwrapping), the cosine and sine of the filtered phase geocoded separately and the phase taken after geocoding; product `data/envisat-2005/t049_20050506_20051028.npz`.

Processor grid: 2960 x 2880 nodes at 0.00042 x 0.00042 degrees, 39.809-41.008 E, 11.742-12.975 N; 41.9% of nodes kept (coherence >= 0.15, median kept coherence 0.41). Product grid: 2160 x 1880 nodes at 0.0005 degrees over 40.2-41.14 E, 11.84-12.92 N, unit phasors and coherence averaged from the processor grid, 46.1% valid.

Look geometry at the scene centre from the reference orbit: heading 192.0 deg, ground-to-satellite azimuth 101.3 deg, incidence 22.6 deg.

Sense check against the forward model of I. Hamling's revised September 2005 opening (455 patches, maximum 16.1 m) in this geometry (model LOS -57 to +220 cm): mean cosine of the wrapped phase-step difference between data and model, for model steps of 0.4-2.6 rad, at lags of 1, 2 and 4 nodes:

| processor phase increases with motion | lag 1 | lag 2 | lag 4 |
|---|---|---|---|
| toward the satellite | +0.363 (463780) | +0.147 (896348) | +0.086 (1331822) |
| away from the satellite | +0.744 (463780) | +0.655 (896348) | +0.587 (1331822) |

The processor's phase increases with motion away from the satellite; the product stores the range change wrapped at one fringe (2.8118 cm), increasing with motion away, the sense of the Hamling track-300 grids (established by the same test on the October 2008 and September 2006 grids).

The agreement is a check of the interferogram's identity and sign, not a fit: the model was derived from other data and the near-axis fringes are unresolved at this spacing.
