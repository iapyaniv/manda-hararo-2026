"""Utility — export analysis-ready GeoTIFFs from the .npz crops.

Per pair, writes two single-band float32 GeoTIFFs (UTM 37N, DEFLATE):
    <tag>.los_disp_m.tif   ionosphere-corrected, far-field-referenced LOS
                           displacement, metres, + = toward the satellite
    <tag>.coherence.tif    interferometric coherence magnitude

Raw unwrapped phase and the geometry cubes stay in the .npz.

Run from the repo root:
    python -m mh26.export_geotiff                # the four labelled pairs
    python -m mh26.export_geotiff --all          # every .npz in data/nisar
"""

from __future__ import annotations

import argparse

import numpy as np
import rasterio
from rasterio.transform import Affine

from .gunw import corrected_los, load_pair, reference_farfield
from .params import DATA_NISAR, PAIRS, REPO

OUT = REPO / "data" / "geotiff"


def write_band(path, arr, x, y, epsg, description):
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])  # negative: y descends
    transform = Affine(dx, 0.0, float(x[0]) - dx / 2,
                       0.0, dy, float(y[0]) - dy / 2)
    with rasterio.open(
        path, "w", driver="GTiff",
        height=arr.shape[0], width=arr.shape[1], count=1,
        dtype="float32", crs=f"EPSG:{epsg}", transform=transform,
        nodata=np.nan, compress="deflate", predictor=3, tiled=True,
    ) as dst:
        dst.write(arr.astype(np.float32), 1)
        dst.set_band_description(1, description)


def export(tag: str) -> None:
    d = load_pair(tag)
    if d is None:
        print(f"{tag}: npz missing, skip")
        return
    los, coh = corrected_los(d)
    ref = reference_farfield(los, d["x"], d["y"])
    epsg = int(d["epsg"])
    OUT.mkdir(parents=True, exist_ok=True)
    write_band(OUT / f"{tag}.los_disp_m.tif", ref, d["x"], d["y"], epsg,
               "LOS displacement, m, + = toward satellite "
               "(ionosphere-corrected, far-field-referenced)")
    # Keep the coherence values only where the analysis LOS is valid, so the
    # exported pair cannot accidentally reintroduce invalid subswath-edge pixels.
    screened_coherence = np.where(np.isfinite(los), coh, np.nan)
    write_band(OUT / f"{tag}.coherence.tif", screened_coherence,
               d["x"], d["y"], epsg,
               "coherence magnitude (analysis validity mask applied)")
    print(f"{tag}: wrote {tag}.los_disp_m.tif + {tag}.coherence.tif")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--all", action="store_true",
                    help="export every .npz in data/nisar, not just the event pairs")
    args = ap.parse_args()
    tags = ([p.stem for p in sorted(DATA_NISAR.glob("*.npz"))]
            if args.all else list(PAIRS))
    for tag in tags:
        export(tag)


if __name__ == "__main__":
    main()
