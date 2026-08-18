# Data package layout

The release is split so the core archive verifies independently and the larger
pre-event NPZ history is optional.

```text
data/
├── nisar/      20 NISAR GUNW crops when both archives are installed
├── s1/         five LiCSAR GeoTIFFs for frame 014A_07688_131313
├── geotiff/    40 analysis-ready NISAR LOS/coherence GeoTIFFs
├── SHA256SUMS-core.txt
├── SHA256SUMS-preevent.txt
└── SHA256SUMS.txt              comprehensive, after both archives
```

## Data archives

From the repository root, install and verify the core archive:

```bash
unzip manda-hararo-2026_data-core.zip -d data
(cd data && sha256sum -c SHA256SUMS-core.txt)
python -m mh26.validate --profile core
```

The core archive contains five NISAR event/control/post-event NPZ crops, all
five Sentinel-1 files, and all 40 compact NISAR GeoTIFF exports. It therefore
supports every reported calculation except regeneration of the full pre-event
history from raw NPZ layers.

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

The Sentinel-1 set contains pair `20260712_20260718` unwrapped phase and
coherence plus frame-level East/North/Up LOS vectors. Source paths and required
attribution are recorded in `../granules/licsar.txt`.

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
