"""Export the July 2026 opening model on the Hamling computational grid.

Ian Hamling can extend his depth-averaged opening figure with the 2026
dike if he has its opening distribution in his own file convention
(2026-09-17). This module samples this repository's fixed-plane fit at
the patch centres of his `plot_3d_dyke_only.dat` grid and writes one
opening value per geometry row, in that file's row order and numeric
format, exactly like a `dyke_only_<epoch>_ex.txt` delivery file.

Sampling is piecewise-constant: each Hamling patch centre is projected
onto the Wright-plane s axis (mh26.hamling.register) and takes the
value of the 2026 model patch containing that (s, depth) position —
the block model as fitted, not an interpolation. Positions outside the
2026 plane (s beyond its SE edge at ~71 km, depth rows below its 10-km
bottom) are zero; the grid's 20-km depth ladder therefore carries the
model in its upper ten rows only. The sampled field integrates on the
Hamling column widths to within a recorded ratio of the native fit
volume (the two grids tile the plane differently); quoted volumes stay
those of the native fit.

Report + file:  python -m mh26.export_2026_opening           (GUNW chain)
                python -m mh26.export_2026_opening --gslc    (30-m GSLC chain)
Outputs: results/opening_2026_hamling_grid[_gslc].{txt,json}.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .distributed import _display_source_path, load_geometry_source
from .hamling import (DATA_DIR, GEOMETRY_FILE, N_COLUMNS, N_DEPTH,
                      _find_fit, _sha256, load_openings, register)
from .params import RESULTS

# The delivery's opening files carry one "%.7e"-formatted value per
# geometry row; the export must be drop-in compatible.
VALUE_FORMAT = "%.7e"
VOLUME_RATIO_BOUNDS = (0.9, 1.1)


def compute_export(chain: str, fit_path: Path, data_dir: Path) -> tuple:
    """Sample the fit on the Hamling grid; return (values, payload)."""
    delivery = load_openings(Path(data_dir).resolve())
    fit_path = Path(fit_path).resolve()
    fit = json.loads(fit_path.read_text())
    if (fit.get("nisar_source") or "gunw") != chain:
        raise SystemExit(f"{fit_path} is not a {chain}-chain fit")
    geometry, fit_prov = load_geometry_source(fit_path)
    if "wright2006_provenance" not in geometry:
        raise SystemExit(f"{fit_path} was not fitted on the Wright plane")
    registration = register(delivery, geometry)

    nz, ns = geometry["dike_shape"]
    opening = np.asarray(fit["dislocation_magnitudes_m"])[:nz * ns]
    opening = opening.reshape(nz, ns)
    s_edges_km = np.asarray(geometry["dike_s_edges_m"], dtype=float) / 1e3
    z_edges_km = np.asarray(geometry["dike_z_edges_m"], dtype=float) / 1e3

    s_km = np.asarray(registration["s_km"])  # one per column, file order
    depth_centres = np.arange(N_DEPTH) + 0.5
    values = np.zeros((N_COLUMNS, N_DEPTH))
    s_index = np.searchsorted(s_edges_km, s_km, side="right") - 1
    z_index = np.searchsorted(z_edges_km, depth_centres, side="right") - 1
    for i in range(N_COLUMNS):
        if not 0 <= s_index[i] < ns:
            continue
        for j in range(N_DEPTH):
            if 0 <= z_index[j] < nz:
                values[i, j] = opening[z_index[j], s_index[i]]

    widths = delivery["columns_width_km"]
    sampled_volume = float((values * widths[:, None]).sum()) * 1e-3
    native_volume = float(fit["dike"]["dike_volume_km3"])
    ratio = sampled_volume / native_volume
    if not VOLUME_RATIO_BOUNDS[0] < ratio < VOLUME_RATIO_BOUNDS[1]:
        raise ValueError(
            f"sampled volume {sampled_volume:.4f} km3 is {ratio:.3f} of "
            f"the native fit volume {native_volume:.4f} km3; the grids "
            "no longer register")

    suffix = "_gslc" if chain == "gslc" else ""
    payload = {
        "generated_by": "mh26.export_2026_opening",
        "nisar_source": chain,
        "file": f"opening_2026_hamling_grid{suffix}.txt",
        "format": "one opening value (m) per row of the geometry file, "
                  f"in its row order, formatted {VALUE_FORMAT!r} — the "
                  "dyke_only_<epoch>_ex.txt convention",
        "geometry_file": GEOMETRY_FILE,
        "geometry_file_sha256": delivery["files"][GEOMETRY_FILE],
        "reference_fit": _display_source_path(fit_path),
        "reference_fit_sha256": _sha256(fit_path),
        "geometry_sha256": fit_prov["geometry_sha256"],
        "sampling": "piecewise-constant: each patch centre takes the "
                    "value of the 2026 model patch containing its "
                    "projected (s, depth) position; zero outside the "
                    "model plane (s past its SE edge, depth below its "
                    "bottom)",
        "plane_span_km": registration["plane_span_km"],
        "plane_bottom_km": float(z_edges_km[-1]),
        "columns_inside_plane": int(
            ((s_index >= 0) & (s_index < ns)).sum()),
        "n_nonzero": int((values > 0).sum()),
        "max_opening_m": float(values.max()),
        "sampled_volume_km3": sampled_volume,
        "native_fit_volume_km3": native_volume,
        "sampled_to_native_ratio": ratio,
    }
    return values, payload


def run(chain: str, fit_path: Path | None, data_dir: Path) -> None:
    fit_path = Path(fit_path or _find_fit(chain)).resolve()
    values, payload = compute_export(chain, fit_path, data_dir)
    suffix = "_gslc" if chain == "gslc" else ""
    txt_path = RESULTS / f"opening_2026_hamling_grid{suffix}.txt"
    np.savetxt(txt_path, values.reshape(-1), fmt=VALUE_FORMAT)
    json_path = RESULTS / f"opening_2026_hamling_grid{suffix}.json"
    json_path.write_text(json.dumps(payload, indent=1, sort_keys=True))
    print(f"wrote {txt_path.name}, {json_path.name}")
    print(f"{payload['n_nonzero']} nonzero of {values.size} values, "
          f"max {payload['max_opening_m']:.2f} m; sampled volume "
          f"{payload['sampled_volume_km3']:.4f} km3 = "
          f"{payload['sampled_to_native_ratio']:.3f} x native "
          f"{payload['native_fit_volume_km3']:.4f}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gslc", action="store_true",
                        help="use the 30-m GSLC-chain fixed-plane fit")
    parser.add_argument("--fit", type=Path,
                        help="explicit Wright-plane fit result JSON")
    parser.add_argument("--dir", type=Path, default=DATA_DIR,
                        help="directory with the Hamling opening files")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    run("gslc" if args.gslc else "gunw", args.fit, args.dir)
