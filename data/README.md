# Data package layout

The release is split so the core archive verifies independently and the larger
pre-event NPZ history is optional.

```text
data/
├── nisar/      20 NISAR GUNW crops when all archives are installed
├── s1/         40 LiCSAR GeoTIFFs for frames 014A_07688_131313 and
│               079D_07694_131313
├── geotiff/    40 analysis-ready NISAR LOS/coherence GeoTIFFs
├── gslc/       4 reprocessed 30-m products; the closure add-on adds raw loop inputs
├── SHA256SUMS-core.txt
├── SHA256SUMS-preevent.txt
├── SHA256SUMS-s1tracks.txt
├── SHA256SUMS-gslc.txt
├── SHA256SUMS-closure.txt      18 exact interferogram/cache files
└── SHA256SUMS.txt              files in the four earlier data archives
```

## Data archives

From the repository root, install and verify the core archive:

```bash
unzip manda-hararo-2026_data-core.zip -d data
(cd data && sha256sum -c SHA256SUMS-core.txt)
```

The validator runs after the two-track add-on below is also installed.

The core archive contains five NISAR event/control/post-event NPZ crops, the
five original ascending Sentinel-1 files, and all 40 compact NISAR GeoTIFF
exports. It is unchanged from the 2026-08-18 release and verifies with the
same checksums.

The two-track Sentinel-1 add-on contains the other 35 LiCSAR files: the
descending frame vectors and the co-event, post-event, and pre-event null
pairs used by `mh26.s1_tracks`. Install it next to the core archive:

```bash
unzip manda-hararo-2026_data-s1tracks.zip -d data
(cd data && sha256sum -c SHA256SUMS-s1tracks.txt)
python -m mh26.validate --profile core
```

The 35 files can also be rebuilt directly from the public LiCSAR paths in
`../granules/licsar.txt`. Core plus s1tracks supports the standard GUNW
analysis. The 30-m GSLC analysis, loop closure, and the full pre-event
history need the separate inputs described below.

The GSLC add-on contains the reprocessed 30-m co-event interferograms
(filtered phase + coherence) and unwrapped LOS fields for both NISAR
geometries, produced by `mh26.gslc` from the 10-m L2 GSLC granules in
`../granules/gslc.tsv`:

```bash
unzip manda-hararo-2026_data-gslc.zip -d data
(cd data && sha256sum -c SHA256SUMS-gslc.txt)
```

The 1.6-GB GSLC crops themselves are not shipped; rebuild everything from
ASF with `python -m mh26.gslc --crop --ifg --unwrap --compare` (Earthdata
Login in `~/.netrc`, `snaphu` on the PATH). Rebuilt unwrapped fields can
differ at the numerical-noise level from the shipped ones if snaphu or
library versions differ; `--compare` reports agreement against the GUNW
fields either way.

The separate `manda-hararo-2026_data-closure.zip` contains the exact nine
interferograms and nine raw unwrap caches needed to reproduce loop
closure. The 18 NPZ files total 608,484,654 bytes. It also contains
`SHA256SUMS-closure.txt`, this README, and `closure-provenance.json`, a copy
of the tracked `results/loop_closure_inventory.json`. The two inversion
pair interferograms also occur in the GSLC archive and have the same
hashes. Install and replay the closure add-on with:

```bash
unzip -n manda-hararo-2026_data-closure.zip -d data
(cd data && sha256sum -c SHA256SUMS-closure.txt)
python -m mh26.closure --check --cached-only
python -m mh26.validate --profile full --require-closure-inputs
```

The raw caches contain no snaphu configuration metadata or original
interferogram hashes. Their original settings, executable identities,
and versions remain unknown. Replay first verifies the saved source,
settings, reports, and input hashes. It then checks each interferogram's
ordered granules, EPSG, looks, array shape, and 30-m grid, and each raw
array's phase congruence. The command and the Python `closure.check()`
API are read-only by default. They reproduce the saved statistics
without running snaphu or replacing files. They do not reconstruct the
original snaphu run. The existing archive checks alone cannot reproduce
closure when these optional inputs are absent. Validation then reports
which tracks could be replayed; `--require-closure-inputs` makes missing
inputs an error.
Partial `--tags` checks also keep the combined reports and inventory
unchanged. After an authorized analysis change, explicitly regenerate
the complete saved set with
`python -m mh26.closure --check --cached-only --regenerate`. This replaces
the reports and their source/settings inventory; it still rejects
changed input hashes. Regeneration requires all three tracks.

The closure archive is prepared locally and has not been uploaded.
Its hash is listed in the repository README and in
`dist/SHA256SUMS-closure.txt`. To build this archive alone after a full
closure replay, run `python scripts/build_release.py --closure-only`.
The builder requires the complete saved baseline to pass validation,
including all input hashes and raw-array replay, before it writes a ZIP.
Stale source, settings, or reports and malformed inventory paths fail
before archive replacement.
The ZIP uses the normal sorted entries, stored compression, source-base
timestamp, and CRC check. Its comment records the base commit, working
tree status, and inventory hash; the inventory hashes the actual code,
model input, reports, and raw inputs and records the analysis settings.

The earlier complex GSLC epoch crops are needed only to form new pre
and skip interferograms. Rebuild those with
`python -m mh26.closure --crop --ifg`. `--check --regenerate` can run snaphu
when a raw cache is missing. New caches record their settings, input hash,
and the executed binary's resolved path, SHA-256, reported version, and
version-command exit status. The binary is resolved once and that exact
path is used for execution. A nonzero `--version` exit status does not
discard a reported version. No such metadata is added to legacy caches.
That produces a new unwrapping run rather than recovering unknown
legacy settings.

Install the 15-crop pre-event add-on when that regeneration is needed:

```bash
unzip manda-hararo-2026_data-preevent.zip -d data
(cd data && sha256sum -c SHA256SUMS-preevent.txt)
(cd data && sha256sum -c SHA256SUMS.txt)
python -m mh26.validate --profile full
```

Each archive contains a checksum manifest for its own files. The core archive
can therefore be verified without the separate pre-event archive. Top-level
archive checksums are in `dist/SHA256SUMS.txt`.

## Hamling dike files (2005 and 2006–2010 opening distributions)

`data/ihamling/` holds Ian Hamling's `dyke_files.zip` delivery (provided
via Tim Wright on 16 September 2026):

```bash
curl -L -o data/ihamling/dyke_files.zip \
  https://sftp.gns.cri.nz/pub/ihamling/dyke_files.zip
# sha256 28a866cb1712c325e0f9b22326d9389c51a58bde5a4b97704dedac4ebd4c7d59
python -c "import zipfile; zipfile.ZipFile('data/ihamling/dyke_files.zip')\
.extractall('data/ihamling')"
```

The 17 small `opening_files/` tables (the `dyke_only_*_ex.txt` openings and
the `plot_3d_dyke_only.dat` patch geometry, ~0.5 MB) are tracked in git, so
`mh26.hamling` runs from a plain clone; `mh26.validate` checks them against
the hashes recorded in `results/uperp_three_epoch*.json`. The zip itself and
the 15 `t300_*_wrp_dyke.grd` grids (GMT classic netCDF, 1024×1024 at 90 m,
UTM 37N, wrapped Envisat track-300 LOS in centimetres — the z range spans
one 2.81-cm fringe) stay out of git; re-download the zip for them. The
decoded table format is documented in the `mh26/hamling.py` docstring.

`data/ihamling/2026-09-17/` holds the two files of Ian Hamling's
follow-up delivery of 17 September 2026 (both tracked): `doplot3d.m`,
his MATLAB plotting script, which indexes an opening file by the row
number of `plot_3d_dyke_only.dat` and thereby confirms the decoded row
order; and `model_slip_sept05.txt`, his updated September 2005 model with
opening limited to 10 km (931 values on an as-yet-undelivered patch
geometry; order-invariant total 2.478 km³ at 1-km² patches). The same
delivery names the canonical thirteen post-2005 files (june09 and apr10
are defunct earlier versions of july09 and may10).

## Contents

The core NISAR set is:

- ascending 8–20 July and descending 7–19 July co-event pairs used for the fit;
- ascending 3–15 July control and 15–27 July straddling pair used for timing;
- descending 19 July–12 August 24-day post-event pair used for continuation.

The 15-pair add-on supplies the remaining November–December 2025 and June–July
2026 history. Exact product identifiers and dates are in
`../granules/pairs.tsv`.

Each NPZ contains `unwrappedPhase`, `coherence`, `ionospherePhaseScreen`,
`connectedComponents`, `mask`, `x`, `y`, `epsg`, `meta`, and local radar-grid
geometry cubes (`cube_*`). Most crops are 1125×1125; track 100 crops are
1125×1046 because the requested 90-km window reaches the product edge. Posting
is 80 m. Rebuilt products can include extra troposphere or solid-Earth-tide
lookup cubes without changing the required analysis schema.

The Sentinel-1 set contains unwrapped phase and coherence for 17 pairs on the
ascending and descending frames plus frame-level East/North/Up LOS vectors
for both frames. Roles (co-event, post-event, pre-event null), source paths,
and required attribution are recorded in `../granules/licsar.txt`.

The GeoTIFF exports contain screened, ionosphere-corrected,
far-field-referenced LOS displacement in metres and coherence masked to the
same analysis-valid pixels. Regenerate them with:

```bash
python -m mh26.export_geotiff --all
```

## Rebuild NISAR crops from ASF

With an Earthdata Login configured in `~/.netrc`:

```bash
python -m mh26.fetch --list
python -m mh26.fetch
```

The fetcher searches the pinned Provisional V1 granule IDs, downloads one
multi-GB HDF5 at a time, extracts the required crop and geometry, writes the
NPZ atomically, then deletes the HDF5. Expect roughly 45 GB of transfer for all
20 products. Locally rebuilt NPZ hashes can differ if ASF republishes a
product or adds optional lookup cubes; the release manifests identify the
exact distributed files.

## 2005 dike-plane geometry (Wright et al. 2006)

`complex_geom2_revised_test.dat` is the source-model geometry of the
September 2005 Dabbahu mega-dike, provided by Juliet Biggs (University of
Bristol) on 26 August 2026. It is the geometry used in Wright et al. (2006,
Nature) and in the Hamling et al. (2009, 2010) studies of the 2006–2010 dike
sequence. Unlike the large data products, it is tiny and is committed to the
repository directly; it is not part of any data archive.

The file is an ASCII table, one source per row, twelve columns, with
coordinates in UTM zone 37N kilometres. It carries no header; the column
meaning below was inferred from the values and is verified numerically every
time `mh26.wright2006` runs (the five dike-segment traces must chain
end-to-end; the observed junction residuals are at most 0.1 m):

| column | meaning |
|---|---|
| 1–3 | strike, dip, rake (degrees) |
| 4–5 | easting, northing of the segment centre (km, UTM 37N) |
| 6 | segment length (km) |
| 7–8 | top and bottom depth (km); for type-3 point rows, column 7 is the source depth and column 8 is unused |
| 9–10 | number of patches along strike and down dip |
| 11 | group tag |
| 12 | source type: 1 graben fault, 2 dike, 3 point source |

The five type-2 rows form the continuous 61-km dike plane (vertical,
0–10 km depth). The type-1 rows are the graben-bounding faults and the
type-3 rows the two point sources of the 2005 model; `mh26.wright2006`
records both as metadata but does not model them. See the module docstring
for how the plane is discretized, extended, and combined with the 2026
fault term.

## Terrain tile for the summary figures

`mh26.summary_figures` draws its map backgrounds and profile topography from
one public Copernicus GLO-30 tile. It is used for display only; no analysis
quantity depends on it. Download it into `data/dem/`:

```bash
mkdir -p data/dem
curl -o data/dem/Copernicus_DSM_COG_10_N12_00_E040_00_DEM.tif \
  https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N12_00_E040_00_DEM/Copernicus_DSM_COG_10_N12_00_E040_00_DEM.tif
```

The module verifies the tile's SHA-256 hash at run time:

    90ad0f4b67afff380496e9e35a0dd74ea4dcb915ab0e64655b28b30cf2da6d24

Copernicus DEM data: © DLR e.V. 2010–2014 and © Airbus Defence and Space
GmbH 2014–2018, provided under COPERNICUS by the European Union and ESA.
