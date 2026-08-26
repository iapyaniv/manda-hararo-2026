# Data package layout

The release is split so the core archive verifies independently and the larger
pre-event NPZ history is optional.

```text
data/
├── nisar/      20 NISAR GUNW crops when all archives are installed
├── s1/         40 LiCSAR GeoTIFFs for frames 014A_07688_131313 and
│               079D_07694_131313
├── geotiff/    40 analysis-ready NISAR LOS/coherence GeoTIFFs
├── gslc/       4 reprocessed 30-m products (ifg + unw per NISAR geometry)
├── SHA256SUMS-core.txt
├── SHA256SUMS-preevent.txt
├── SHA256SUMS-s1tracks.txt
├── SHA256SUMS-gslc.txt
└── SHA256SUMS.txt              all data files (all four archives)
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
`../granules/licsar.txt`. Core plus s1tracks supports every reported
calculation except regeneration of the full pre-event history from raw NPZ
layers.

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
| 7–8 | top and bottom depth (km) |
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
