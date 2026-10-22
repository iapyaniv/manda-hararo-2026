"""rebuild data/nisar from the open NISAR archive at ASF DAAC.

For each row of `granules/pairs.tsv` (tag <TAB> granule ID), this command
downloads the L2 GUNW HDF5 through earthaccess, extracts a 90 × 90 km crop,
writes `data/nisar/<tag>.npz`, and deletes the downloaded HDF5 unless `--keep`
is specified. Earthdata Login as for `mh26.gslc`: a user token in
EARTHDATA_TOKEN or ~/.config/earthdata/token, else ~/.netrc.

Run from the repo root:
    python -m mh26.fetch --list            # check archive availability only
    python -m mh26.fetch                   # fetch everything missing (~50 GB traffic)
    python -m mh26.fetch --tags 172_008_A 165_082_D
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import earthaccess
import h5py
import numpy as np

from .gslc import earthdata_login
from .params import CENTER_E, CENTER_N, DATA_NISAR, REPO

COLLECTION = "NISAR_L2_GUNW_PROVISIONAL_V1"
PAIRS_TSV = REPO / "granules" / "pairs.tsv"
TMP = REPO / "data" / "h5tmp"
CROP_HALF = 45_000.0  # m

GRIDS = "/science/LSAR/GUNW/grids/frequencyA"
RADARGRID = "/science/LSAR/GUNW/metadata/radarGrid"
IDENT = "/science/LSAR/identification"
# Output key -> source HDF5 dataset name. The explicit coherence mapping is
# required because the product dataset is named ``coherenceMagnitude``.
LAYERS = {
    "unwrappedPhase": "unwrappedPhase",
    "coherence": "coherenceMagnitude",
    "ionospherePhaseScreen": "ionospherePhaseScreen",
    "connectedComponents": "connectedComponents",
    "mask": "mask",
}
CUBES = ("incidenceAngle", "losUnitVectorX", "losUnitVectorY",
         "xCoordinates", "yCoordinates", "heightAboveEllipsoid",
         "wetTroposphericPhaseScreen", "hydrostaticTroposphericPhaseScreen",
         "slantRangeSolidEarthTidesPhase")


def read_pairs() -> dict[str, str]:
    pairs = {}
    for line_number, line in enumerate(PAIRS_TSV.read_text().splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 2 or not all(fields):
            raise ValueError(f"{PAIRS_TSV}:{line_number}: expected tag<TAB>granule")
        tag, granule = fields
        if tag in pairs:
            raise ValueError(f"{PAIRS_TSV}:{line_number}: duplicate tag {tag!r}")
        pairs[tag] = granule
    return pairs


def detect_pol(f: h5py.File) -> str:
    base = f[f"{GRIDS}/unwrappedInterferogram"]
    for pol in ("HH", "VV"):
        if pol in base:
            return pol
    raise KeyError(f"no co-pol group found; groups: {list(base.keys())}")


def extract(h5path: Path, tag: str) -> Path:
    """Crop the GUNW layers + geometry cubes around the anomaly into an .npz."""
    with h5py.File(h5path, "r") as f:
        pol = detect_pol(f)
        base = f"{GRIDS}/unwrappedInterferogram"
        x = f[f"{base}/{pol}/xCoordinates"][:]
        y = f[f"{base}/{pol}/yCoordinates"][:]
        ix = np.where((x >= CENTER_E - CROP_HALF) & (x <= CENTER_E + CROP_HALF))[0]
        iy = np.where((y >= CENTER_N - CROP_HALF) & (y <= CENTER_N + CROP_HALF))[0]
        if len(ix) == 0 or len(iy) == 0:
            raise SystemExit(f"{tag}: anomaly window outside grid")
        sl = np.s_[iy[0]: iy[-1] + 1, ix[0]: ix[-1] + 1]

        data = {}
        for key, dataset in LAYERS.items():
            path = f"{base}/mask" if key == "mask" else f"{base}/{pol}/{dataset}"
            if path not in f:
                raise KeyError(f"{tag}: required HDF5 dataset missing: {path}")
            data[key] = f[path][sl]
        analysis_shape = data["unwrappedPhase"].shape
        bad_shapes = {key: value.shape for key, value in data.items()
                      if key in LAYERS and value.shape != analysis_shape}
        if bad_shapes:
            raise ValueError(f"{tag}: analysis-layer shape mismatch: {bad_shapes}")
        data["x"] = x[ix[0]: ix[-1] + 1]
        data["y"] = y[iy[0]: iy[-1] + 1]
        for key in ("projection", "epsg"):
            p = f"{base}/{pol}/{key}"
            if p in f:
                data["epsg"] = np.array(int(np.array(f[p]).item()))
                break
        if "epsg" not in data:
            raise KeyError(f"{tag}: projection/EPSG dataset missing")

        ident = f[IDENT]
        meta = {}
        for k in ("orbitPassDirection", "referenceZeroDopplerStartTime",
                  "secondaryZeroDopplerStartTime", "trackNumber", "frameNumber"):
            if k in ident:
                v = ident[k][()]
                if isinstance(v, bytes):
                    v = v.decode()
                elif isinstance(v, np.generic):
                    v = v.item()
                meta[k] = v
        data["meta"] = np.array(json.dumps(meta, sort_keys=True))

        rg = f[RADARGRID]
        for cube in CUBES:
            if cube in rg:
                data[f"cube_{cube}"] = rg[cube][:]

    DATA_NISAR.mkdir(parents=True, exist_ok=True)
    out = DATA_NISAR / f"{tag}.npz"
    tmp = out.with_name(out.name + ".tmp")
    with tmp.open("wb") as stream:
        np.savez_compressed(stream, **data)
    tmp.replace(out)
    print(f"  {tag}: saved {out.name} ({out.stat().st_size/1e6:.1f} MB), "
          f"crop {data['unwrappedPhase'].shape}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list", action="store_true",
                    help="search the archive and report availability, download nothing")
    ap.add_argument("--keep", action="store_true", help="keep downloaded .h5 files")
    ap.add_argument("--tags", nargs="*", help="only these tags (default: all missing)")
    args = ap.parse_args()

    pairs = read_pairs()
    if args.tags:
        missing_tags = [t for t in args.tags if t not in pairs]
        if missing_tags:
            raise SystemExit(f"unknown tags (not in pairs.tsv): {missing_tags}")
        pairs = {t: pairs[t] for t in args.tags}

    earthdata_login()
    for tag, granule in pairs.items():
        if not args.list and (DATA_NISAR / f"{tag}.npz").exists():
            print(f"{tag}: already extracted")
            continue
        hits = earthaccess.search_data(short_name=COLLECTION, granule_name=granule)
        if not hits:
            print(f"{tag}: {granule} NOT FOUND in archive")
            continue
        if args.list:
            print(f"{tag}: available ({granule})")
            continue
        TMP.mkdir(parents=True, exist_ok=True)
        files = earthaccess.download(hits[:1], str(TMP))
        h5 = Path(files[0])
        extract(h5, tag)
        if not args.keep:
            h5.unlink()


if __name__ == "__main__":
    main()
