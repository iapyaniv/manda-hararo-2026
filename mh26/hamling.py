"""Three-epoch opening comparison from the Hamling dike files.

Ian Hamling (GNS) provided `dyke_files.zip` (linked by Tim Wright in the
paper thread on 2026-09-16), whose `opening_files/` holds the published
opening distributions of the September 2005 megadike and fifteen later
(2006-2010) opening files on the ~74 km Hamling computational grid
(Hamling et al. 2009, 2010; the published sequence has 13 post-megadike
intrusions and the data flag two files as duplicate content — see
`version_diagnostics`, `PRIMARY_EXCLUDED`, and the report).
This module decodes that delivery, registers the grid against
the Wright et al. (2006) plane the 2026 fit uses (mh26.wright2006), and
renders the three-epoch U_perp comparison the plane was built for:
September 2005, the cumulative 2006-2010 sequence, and July 2026 on one
along-plane axis.

Decoded format (verified by the structural checks in `load_openings`):

- `plot_3d_dyke_only.dat`: 1500 rows x 9 columns, one fault patch per row:
  easting (m), northing (m) [UTM 37N], strike (deg), dip (90), rake-slot
  (180), example opening (m, unused here), along-strike patch length (km),
  patch top (km), patch bottom (km). Rows are depth-fastest: 20 one-km
  depth rows (0-20 km) per along-strike column, 75 columns ordered SE->NW,
  in six constant-strike segments (329.7-355.7 deg, ~0.94-1.14 km column
  widths, ~73.9 km trace).
- `dyke_only_<epoch>_ex.txt`: 1500 opening values in metres, one per patch.
  The files carry no coordinates; they are taken to follow the
  `plot_3d_dyke_only.dat` row order. That assumption rests on the equal row
  count and on the physically coherent profiles it produces (smooth depth
  profiles; the late small dikes vanish below ~8 km) and is flagged for
  Ian's confirmation in LIMITATIONS.md.
- `t300_*_wrp_dyke.grd`: GMT classic-netCDF grids (1024x1024, 90 m, UTM
  37N) of wrapped Envisat track-300 LOS in centimetres (z spans one fringe,
  0-2.81 cm = lambda/2); context only, not used here.

The openings are the published Envisat-era inversions exactly as provided:
no uncertainties were delivered, and their regularization and greens
functions differ from the 2026 fit's. The comparison is therefore of
published models on a shared axis, not a joint inversion.

Registration: each Hamling column centre is projected into the chained
segment frames of the 2026 contract (mh26.wright2006._project_s), giving a
monotone along-plane coordinate s. The two traces are close but not
identical surfaces (median off-plane offset ~0.1 km, up to ~1.4 km where
the Wright segmentation and the Hamling grid disagree mid-plane), so the
shared s axis compares openings by along-axis position, not patch by patch.

Report + figure:  python -m mh26.hamling           (GUNW-chain 2026 fit)
                  python -m mh26.hamling --gslc    (30-m GSLC-chain fit)
Outputs: results/uperp_three_epoch[_gslc].{json,md,png,pdf}.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from .distributed import _display_source_path, load_geometry_source
from .params import REPO, RESULTS
from .wright2006 import _chain_frames, _project_s

DATA_DIR = REPO / "data" / "ihamling" / "opening_files"
GEOMETRY_FILE = "plot_3d_dyke_only.dat"
SOURCE_URL = "https://sftp.gns.cri.nz/pub/ihamling/dyke_files.zip"
RECEIVED = ("linked by Tim Wright in the paper thread on 2026-09-16 "
            "(files produced by Ian Hamling, GNS)")

N_COLUMNS = 75
N_DEPTH = 20

# Filename month tags in event order. The delivery carries no exact dates;
# event dates must come from the Hamling et al. (2009, 2010) catalogue when
# captioning, so only the tags are recorded here.
EPOCH_TAGS = ("sept05", "june06", "july06", "sept06", "dec06", "jan07",
              "aug07", "nov07", "apr08", "july08", "oct08", "feb09",
              "june09", "july09", "apr10", "may10")

# The delivery holds 15 post-2005 files against the 13 published
# post-megadike intrusions, and two files duplicate others' content:
# may10 equals apr10 with every cell below 0.5 m set to zero (an exact,
# deterministic relation), and july09's positive cells all lie inside
# june09's with a 0.994 correlation, against at most 0.44 between any
# other pair — two estimates of one event, not two events. The primary
# cumulative therefore excludes the two files below; the report quotes the
# other candidate sums, and `version_diagnostics` re-derives the evidence
# from the data so a changed delivery cannot silently keep this decision.
# The file-to-event mapping still awaits Ian Hamling's confirmation.
PRIMARY_EXCLUDED = ("july09", "may10")

# The 13 post-2005 intrusions of the published sequence, used only for the
# alternative catalogue-count sum in the report (note it represents the
# May 2010 event by the filtered may10 field).
CATALOGUE_13 = ("june06", "july06", "sept06", "dec06", "jan07", "aug07",
                "nov07", "apr08", "july08", "oct08", "feb09", "june09",
                "may10")

OPENING_CMAP = LinearSegmentedColormap.from_list(
    "hamling_opening",
    [
        (0.00, "#ffffff"),
        (0.12, "#fff7bc"),
        (0.38, "#fec44f"),
        (0.68, "#f03b20"),
        (1.00, "#b1006e"),
    ],
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_openings(data_dir: Path = DATA_DIR) -> dict:
    """Parse and structurally verify the delivery; raise on any mismatch."""
    geometry_path = data_dir / GEOMETRY_FILE
    table = np.loadtxt(geometry_path)
    if table.shape != (N_COLUMNS * N_DEPTH, 9):
        raise ValueError(f"{geometry_path}: expected "
                         f"{N_COLUMNS * N_DEPTH} rows x 9 columns, "
                         f"got {table.shape}")
    blocks = table.reshape(N_COLUMNS, N_DEPTH, 9)
    if not np.allclose(blocks[:, :, :2], blocks[:, :1, :2]):
        raise ValueError(f"{geometry_path}: rows are not depth-fastest "
                         "(easting/northing varies inside a 20-row block)")
    tops, bottoms = blocks[0, :, 7], blocks[0, :, 8]
    if not (np.array_equal(tops, np.arange(N_DEPTH, dtype=float))
            and np.array_equal(bottoms, tops + 1.0)):
        raise ValueError(f"{geometry_path}: depth rows are not the "
                         "one-kilometre 0-20 km ladder")
    if not (np.allclose(blocks[:, :, 3], 90.0)
            and np.allclose(blocks[:, :, 7:9], blocks[:1, :, 7:9])):
        raise ValueError(f"{geometry_path}: dips or depth ladders vary "
                         "between columns")
    columns = blocks[:, 0, :]  # one row per along-strike column
    if not np.all(np.diff(columns[:, 1]) > 0):
        raise ValueError(f"{geometry_path}: columns are not ordered S->N")

    files = {GEOMETRY_FILE: _sha256(geometry_path)}
    openings = {}
    for tag in EPOCH_TAGS:
        path = data_dir / f"dyke_only_{tag}_ex.txt"
        values = np.loadtxt(path)
        if values.shape != (N_COLUMNS * N_DEPTH,):
            raise ValueError(f"{path}: expected {N_COLUMNS * N_DEPTH} "
                             f"values, got {values.shape}")
        if values.min() < 0.0:
            raise ValueError(f"{path}: negative opening")
        openings[tag] = values.reshape(N_COLUMNS, N_DEPTH)
        files[path.name] = _sha256(path)
    extra = sorted(p.name for p in data_dir.glob("dyke_only_*_ex.txt")
                   if p.name not in files)
    if extra:
        raise ValueError(f"unrecognized epoch files {extra}; extend "
                         "EPOCH_TAGS after checking their dates")

    return {
        "columns_e_m": columns[:, 0],
        "columns_n_m": columns[:, 1],
        "columns_strike_deg": columns[:, 2],
        "columns_width_km": columns[:, 6],
        "depth_edges_km": np.arange(N_DEPTH + 1, dtype=float),
        "openings_m": openings,
        "files": files,
    }


def register(delivery: dict, geometry: dict) -> dict:
    """Project the Hamling columns into the 2026 contract's s axis (km)."""
    prov = geometry["wright2006_provenance"]
    ext = prov["extension_se"]
    frames = _chain_frames(prov["segments_nw_to_se"],
                           ext["n_patches"] * ext["patch_width_m"])
    s_km = np.empty(N_COLUMNS)
    offplane_m = np.empty(N_COLUMNS)
    for i, (e, n) in enumerate(zip(delivery["columns_e_m"],
                                   delivery["columns_n_m"])):
        s, perp = _project_s(frames, np.array([n, e]))
        s_km[i] = s / 1e3
        offplane_m[i] = perp
    if not np.all(np.diff(s_km) < 0):
        raise ValueError("projected s is not monotone along the Hamling "
                         "grid; the traces do not share one axis")
    plane_span = (float(geometry["dike_s_edges_m"][0]) / 1e3,
                  float(geometry["dike_s_edges_m"][-1]) / 1e3)
    inside = (s_km >= plane_span[0]) & (s_km <= plane_span[1])
    return {
        "frame": "along-plane s of the mh26.wright2006 contract "
                 "(km, 0 at the NW end of the original plane)",
        "s_km": s_km,
        "offplane_m": offplane_m,
        "median_offplane_m": float(np.median(offplane_m)),
        "max_offplane_m": float(offplane_m.max()),
        "span_km": [float(s_km.min()), float(s_km.max())],
        "plane_span_km": list(plane_span),
        "columns_within_plane": int(inside.sum()),
        "se_limit_2005_km": float(
            geometry["dike_s_edges_m"][
                geometry["dike_shape"][1] - ext["n_patches"]]) / 1e3,
    }


def _epoch_stats(opening: np.ndarray, width_km: np.ndarray) -> dict:
    area_km2 = width_km[:, None] * 1.0  # one-km depth rows
    volume = float((opening * area_km2).sum()) * 1e-3
    above10 = float((opening[:, :10] * area_km2[:, :10]).sum()) * 1e-3
    return {
        "volume_km3": round(volume, 4),
        "volume_above_10km_km3": round(above10, 4),
        "max_opening_m": round(float(opening.max()), 3),
    }


def version_diagnostics(openings: dict) -> list[dict]:
    """Post-2005 file pairs whose content duplicates another's.

    Two independent intrusions on the same plane share little: no
    distinct-event pair in this delivery correlates above |r| = 0.44.
    A pair is flagged when the smaller file's positive cells all lie
    inside the larger's and either the shared values are exactly equal
    (a thresholded copy — the threshold is recovered and recorded; the
    removed small cells depress the correlation, so exactness alone
    suffices) or the whole fields correlate at r >= 0.98 (a re-estimate
    of the same event). Summing a flagged pair counts one event twice.
    """
    pairs = []
    for tag in EPOCH_TAGS[1:]:
        contained = openings[tag]
        nz = contained > 0
        if not nz.any():
            continue
        for other in EPOCH_TAGS[1:]:
            if other == tag:
                continue
            container = openings[other]
            nz_other = container > 0
            if nz_other.sum() <= nz.sum() or not np.all(nz_other[nz]):
                continue
            corr = float(np.corrcoef(contained.ravel(),
                                     container.ravel())[0, 1])
            exact = bool(np.array_equal(contained[nz], container[nz]))
            if not exact and corr < 0.98:
                continue
            threshold = None
            if exact:
                floor = float(contained[nz].min())
                if np.array_equal(contained,
                                  np.where(container >= floor, container,
                                           0.0)):
                    threshold = floor
            pairs.append({"contained": tag, "container": other,
                          "n_cells": int(nz.sum()),
                          "correlation": round(corr, 4),
                          "exact_value_subset": exact,
                          "threshold_m": threshold})
    return pairs


def _fit_opening_stats(fit: dict, geometry: dict) -> dict:
    """Integrate the fitted opening on its own grid, split at 10 km.

    A depth row straddling 10 km contributes its fractional overlap. The
    0.5-m-floor volume is a uniform-floor sensitivity, mirroring the one
    floor recovered from the delivered may10/apr10 pair; the papers
    filtered by local uncertainty, so it is not a matched reconstruction
    of their filtering.
    """
    nz, ns = geometry["dike_shape"]
    opening = np.asarray(fit["dislocation_magnitudes_m"])[:nz * ns]
    opening = opening.reshape(nz, ns)
    s_edges = np.asarray(geometry["dike_s_edges_m"], dtype=float)
    z_edges = np.asarray(geometry["dike_z_edges_m"], dtype=float)
    areas_km2 = (np.diff(s_edges)[None, :] / 1e3
                 * np.diff(z_edges)[:, None] / 1e3)
    z0, z1 = z_edges[:-1] / 1e3, z_edges[1:] / 1e3
    above_frac = np.clip((10.0 - z0) / (z1 - z0), 0.0, 1.0)
    floored = np.where(opening >= 0.5, opening, 0.0)
    return {
        "volume_km3": round(fit["dike"]["dike_volume_km3"], 4),
        "volume_above_10km_km3": round(
            float((opening * areas_km2 * above_frac[:, None]).sum())
            * 1e-3, 4),
        "volume_at_0p5m_floor_km3": round(
            float((floored * areas_km2).sum()) * 1e-3, 4),
        "max_opening_m": round(fit["dike"]["max_opening_m"], 3),
        "plane_bottom_km": float(z_edges[-1]) / 1e3,
    }


def _panel(ax, s_centres, z_edges, opening, vmax, line_levels, title,
           note):
    """One filled U_perp cross-section in the mh26.wright2006 panel style."""
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    s_lo = s_centres[0] - 0.5 * (s_centres[1] - s_centres[0])
    s_hi = s_centres[-1] + 0.5 * (s_centres[-1] - s_centres[-2])
    s_plot = np.concatenate(([s_lo], s_centres, [s_hi]))
    z_plot = np.concatenate(([z_edges[0]], z_centres, [z_edges[-1]]))
    padded = np.pad(opening, 1, mode="edge")
    filled = ax.contourf(s_plot, z_plot, padded,
                         levels=np.linspace(0.0, vmax, 21),
                         cmap=OPENING_CMAP,
                         extend="max" if opening.max() > vmax else "neither",
                         antialiased=True)
    filled.set_rasterized(True)
    levels = [lv for lv in line_levels if lv < float(opening.max())]
    if levels:
        ax.contour(s_plot, z_plot, padded, levels=levels, colors="0.15",
                   linewidths=0.55)
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", pad=5)
    ax.set_title(note, loc="right", fontsize=8, color="0.30", pad=6)
    ax.set_ylabel("Depth (km)")
    return filled


def three_epoch_figure(delivery: dict, registration: dict, fit: dict,
                       geometry: dict, chain: str, epochs: dict,
                       cumulative: np.ndarray) -> Path:
    s_ham = registration["s_km"][::-1]  # NW->SE for plotting
    panels = [
        ("September 2005", delivery["openings_m"]["sept05"][::-1].T,
         s_ham, np.arange(21.0), 13.0, [2, 4, 6, 8, 10, 12],
         "Hamling dike files, as provided",
         epochs["sept05"]),
        ("2006–2010 (cumulative, provisional set)", cumulative[::-1].T,
         s_ham, np.arange(21.0), 9.0, [1, 2, 4, 6, 8],
         "Hamling dike files; two-file mapping pending",
         epochs["cumulative_2006_2010"]),
    ]

    nz, ns = geometry["dike_shape"]
    x = np.asarray(fit["dislocation_magnitudes_m"])
    opening_2026 = x[:nz * ns].reshape(nz, ns)
    s_edges = geometry["dike_s_edges_m"] / 1e3
    z_edges = geometry["dike_z_edges_m"] / 1e3
    s_2026 = 0.5 * (s_edges[:-1] + s_edges[1:])
    source_note = ("this study, 30-m GSLC chain" if chain == "gslc"
                   else "this study, GUNW chain")
    panels.append(("July 2026", opening_2026, s_2026, z_edges, 2.5,
                   [0.5, 1.0, 1.5, 2.0], source_note, epochs["epoch_2026"]))

    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.2), sharex=True,
                             constrained_layout=True)
    xlim = (0.0, max(registration["span_km"][1] + 0.6,
                     registration["plane_span_km"][1]))
    se_2005 = registration["se_limit_2005_km"]
    for ax, (title, opening, s_centres, z_ax, vmax, lines, note,
             stats) in zip(axes, panels):
        filled = _panel(ax, s_centres, z_ax, opening, vmax, lines, title,
                        note)
        ax.set_xlim(*xlim)
        ax.set_ylim(0.0, 20.0)
        ax.invert_yaxis()
        ax.axvline(se_2005, color="0.25", lw=0.9, ls=(0, (4, 3)), zorder=4)
        label = (f"{stats['volume_km3']:.2f} km$^3$"
                 + (f" ({stats['volume_above_10km_km3']:.2f} above 10 km)"
                    if stats["volume_km3"]
                    - stats["volume_above_10km_km3"] > 5e-3 else "")
                 + f" · max {stats['max_opening_m']:.1f} m")
        ax.text(0.012, 0.06, label, transform=ax.transAxes, ha="left",
                va="bottom", fontsize=8.5,
                bbox={"facecolor": "white", "edgecolor": "none",
                      "alpha": 0.8, "pad": 1.6}, zorder=6)
        ticks = np.arange(0.0, vmax + 0.01, max(1.0, round(vmax / 5)))
        if vmax <= 3.0:
            ticks = np.arange(0.0, vmax + 0.01, 0.5)
        colorbar = fig.colorbar(filled, ax=ax, label="U$_\\perp$ (m)",
                                ticks=ticks, shrink=0.94, pad=0.015)
        colorbar.ax.tick_params(labelsize=8)

    ax_2026 = axes[2]
    ax_2026.axhspan(z_edges[-1], 20.0, color="0.88", zorder=1)
    if xlim[1] > s_edges[-1]:
        # Beyond the plane's SE edge there is no 2026 model; grey it so it
        # does not read as zero opening.
        ax_2026.axvspan(s_edges[-1], xlim[1], color="0.88", zorder=1)
    ax_2026.text(0.5 * sum(xlim), 0.5 * (z_edges[-1] + 20.0),
                 f"below the 2026 model plane ({z_edges[-1]:g} km)",
                 ha="center", va="center", fontsize=8, color="0.35",
                 zorder=2)
    ax_2026.text(se_2005 + 0.4, 0.06, "2005 SE limit",
                 transform=ax_2026.get_xaxis_transform(), rotation=90,
                 ha="left", va="top", fontsize=7.5, color="0.20",
                 bbox={"facecolor": "white", "edgecolor": "none",
                       "alpha": 0.78, "pad": 1.2}, zorder=5)
    ax_2026.set_xlabel("Distance along the 2005 dike plane, NW→SE (km)")

    suffix = "_gslc" if chain == "gslc" else ""
    out = RESULTS / f"uperp_three_epoch{suffix}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=250)
    plt.close(fig)
    return Path(f"{out}.png")


def _find_fit(chain: str) -> Path:
    matches = []
    for path in sorted(RESULTS.glob("distributed_fit_geometry_*.json")):
        is_gslc = path.stem.endswith("_gslc")
        if (chain == "gslc") != is_gslc:
            continue
        payload = json.loads(path.read_text())
        if "wright2006_provenance" in payload.get("geometry", {}):
            matches.append(path)
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one Wright-plane fit for the "
                         f"{chain} chain, found {matches}")
    return matches[0]


def run(chain: str, fit_path: Path | None, data_dir: Path) -> None:
    data_dir = Path(data_dir).resolve()
    delivery = load_openings(data_dir)
    fit_path = Path(fit_path or _find_fit(chain)).resolve()
    fit = json.loads(fit_path.read_text())
    if (fit.get("nisar_source") or "gunw") != chain:
        raise SystemExit(f"{fit_path} is not a {chain}-chain fit")
    geometry, fit_prov = load_geometry_source(fit_path)
    if "wright2006_provenance" not in geometry:
        raise SystemExit(f"{fit_path} was not fitted on the Wright plane")
    registration = register(delivery, geometry)

    epochs = {tag: _epoch_stats(delivery["openings_m"][tag],
                                delivery["columns_width_km"])
              for tag in EPOCH_TAGS}
    diagnostics = version_diagnostics(delivery["openings_m"])
    flagged = sorted({d["contained"] for d in diagnostics})
    if flagged != sorted(PRIMARY_EXCLUDED):
        raise ValueError(
            f"data-derived duplicate flags {flagged} no longer match the "
            f"reviewed exclusions {sorted(PRIMARY_EXCLUDED)}; the delivery "
            "changed — revisit PRIMARY_EXCLUDED before summing")
    cum_tags = [t for t in EPOCH_TAGS[1:] if t not in PRIMARY_EXCLUDED]
    cum = sum(delivery["openings_m"][t] for t in cum_tags)
    widths = delivery["columns_width_km"]
    epochs["cumulative_2006_2010"] = _epoch_stats(cum, widths)
    floored_cum = sum(np.where(delivery["openings_m"][t] >= 0.5,
                               delivery["openings_m"][t], 0.0)
                      for t in cum_tags)
    epochs["cumulative_2006_2010"]["volume_at_0p5m_floor_km3"] = \
        _epoch_stats(floored_cum, widths)["volume_km3"]
    epochs["epoch_2026"] = _fit_opening_stats(fit, geometry)

    figure = three_epoch_figure(delivery, registration, fit, geometry,
                                chain, epochs, cum)

    suffix = "_gslc" if chain == "gslc" else ""
    payload = {
        "generated_by": "mh26.hamling",
        "nisar_source": chain,
        "source": {
            "url": SOURCE_URL,
            "received": RECEIVED,
            "directory": _display_source_path(data_dir),
            "files": delivery["files"],
        },
        "layout": {
            "n_columns": N_COLUMNS,
            "n_depth_rows": N_DEPTH,
            "depth_range_km": [0.0, 20.0],
            "row_order": "depth-fastest (20 one-km rows per column), "
                         "columns S->N in the file",
            "opening_row_order_assumed_from_geometry_file": True,
            "column_width_km": [
                float(delivery["columns_width_km"].min()),
                float(delivery["columns_width_km"].max())],
            "crs": "UTM 37N metres",
        },
        "registration": {
            key: (value.tolist() if isinstance(value, np.ndarray) else value)
            for key, value in registration.items()
        },
        "epochs": {tag: epochs[tag] for tag in EPOCH_TAGS},
        "cumulative_2006_2010": {
            "files_summed": cum_tags,
            "primary_excluded": list(PRIMARY_EXCLUDED),
            "pair_diagnostics": diagnostics,
            **epochs["cumulative_2006_2010"],
            "alternatives": {
                "sum_of_all_delivered_files": _epoch_stats(
                    sum(delivery["openings_m"][t] for t in EPOCH_TAGS[1:]),
                    widths),
                "counted_once_14_files": _epoch_stats(
                    sum(delivery["openings_m"][t] for t in EPOCH_TAGS[1:]
                        if t != "may10"), widths),
                "keep_july09_drop_june09": _epoch_stats(
                    sum(delivery["openings_m"][t] for t in EPOCH_TAGS[1:]
                        if t not in ("june09", "may10")), widths),
                "drop_june09_and_apr10": _epoch_stats(
                    sum(delivery["openings_m"][t] for t in EPOCH_TAGS[1:]
                        if t not in ("june09", "apr10")), widths),
                "catalogue_13_events": {
                    "files": list(CATALOGUE_13),
                    **_epoch_stats(
                        sum(delivery["openings_m"][t] for t in CATALOGUE_13),
                        widths),
                },
            },
        },
        "epoch_2026": {
            "reference_fit": _display_source_path(fit_path),
            "reference_fit_sha256": _sha256(fit_path),
            "geometry_sha256": fit_prov["geometry_sha256"],
            **epochs["epoch_2026"],
        },
        "notes": [
            "openings are the published Envisat-era inversions exactly as "
            "provided (no uncertainties; different regularization and "
            "greens functions from the 2026 fit)",
            "the epoch-file row order is assumed to follow "
            "plot_3d_dyke_only.dat and awaits Ian Hamling's confirmation",
            "the delivery has 15 post-2005 files against the 13 published "
            "post-megadike intrusions; may10 is apr10 with cells below "
            "0.5 m removed, and july09's positive cells all lie inside "
            "june09's at r=0.994 (all other pairs |r|<=0.44), so the "
            "primary cumulative excludes both and the report quotes the "
            "candidate sums; the file-to-event mapping awaits Ian "
            "Hamling's confirmation",
            "the delivered historical fields were filtered below local "
            "uncertainty per the papers (the recovered 0.5-m floor is "
            "specific to the may10/apr10 pair), while the 2026 totals "
            "keep all fitted opening; the volume_at_0p5m_floor_km3 "
            "values apply one uniform 0.5-m floor to the primary "
            "cumulative and to 2026 as a sensitivity, not as a matched "
            "reconstruction of that filtering",
            "the shared s axis compares openings by along-axis position; "
            "the two traces differ by up to "
            f"{registration['max_offplane_m']:.0f} m mid-plane",
        ],
        "figure": figure.name,
    }
    json_path = RESULTS / f"uperp_three_epoch{suffix}.json"
    json_path.write_text(json.dumps(payload, indent=1, sort_keys=True))

    lines = [
        "# Three-epoch opening comparison (Hamling dike files)",
        "",
        f"NISAR chain of the 2026 fit: {chain} "
        f"({payload['epoch_2026']['reference_fit']}).",
        "",
        "Openings of the September 2005 megadike and the fifteen "
        "post-2005 opening files are the published inversions provided by "
        "Ian Hamling (GNS) on the ~74 km Hamling computational grid "
        "(75 columns x 20 one-km depth rows, UTM 37N). The 2026 "
        "distribution is this repository's fixed-plane fit on the Wright "
        "et al. (2006) plane. The Hamling columns project monotonically "
        "into that plane's s axis (median off-plane offset "
        f"{registration['median_offplane_m']:.0f} m, maximum "
        f"{registration['max_offplane_m']:.0f} m mid-plane), so "
        f"`{figure.name}` compares the three epochs by along-axis "
        "position on one axis.",
        "",
        "| epoch | volume (km3) | above 10 km (km3) | max opening (m) |",
        "|---|---|---|---|",
    ]
    for tag in EPOCH_TAGS:
        stats = epochs[tag]
        lines.append(f"| {tag} | {stats['volume_km3']:.3f} | "
                     f"{stats['volume_above_10km_km3']:.3f} | "
                     f"{stats['max_opening_m']:.2f} |")
    stats = epochs["cumulative_2006_2010"]
    alternatives = payload["cumulative_2006_2010"]["alternatives"]
    lines.append(f"| 2006-2010 cumulative ({len(cum_tags)} files, primary) "
                 f"| {stats['volume_km3']:.3f} | "
                 f"{stats['volume_above_10km_km3']:.3f} | "
                 f"{stats['max_opening_m']:.2f} |")
    for label, alt in (
            (f"all {len(EPOCH_TAGS) - 1} delivered files",
             alternatives["sum_of_all_delivered_files"]),
            ("14 files (only may10 dropped)",
             alternatives["counted_once_14_files"]),
            ("keep july09, drop june09",
             alternatives["keep_july09_drop_june09"]),
            ("keep july09 + may10, drop june09 + apr10",
             alternatives["drop_june09_and_apr10"]),
            ("catalogue 13 events (filtered may10)",
             alternatives["catalogue_13_events"])):
        lines.append(f"| — alternative: {label} | "
                     f"{alt['volume_km3']:.3f} | "
                     f"{alt['volume_above_10km_km3']:.3f} | "
                     f"{alt['max_opening_m']:.2f} |")
    fit_stats = epochs["epoch_2026"]
    lines.append(f"| July 2026 ({chain}) | {fit_stats['volume_km3']:.3f} | "
                 f"{fit_stats['volume_above_10km_km3']:.3f} | "
                 f"{fit_stats['max_opening_m']:.2f} |")
    bottom = fit_stats["plane_bottom_km"]
    lines += [
        "",
        "The delivery has 15 post-2005 files against the 13 published "
        "post-megadike intrusions, and the data flag two as duplicate "
        "content: may10 equals apr10 with every cell below 0.5 m set to "
        "zero (38% of that field's volume), and july09's positive cells "
        "all lie inside june09's at correlation 0.994 while no other "
        "pair exceeds |r| = 0.44 — two estimates of one event. The "
        "primary cumulative therefore sums 13 files (drops may10 and "
        "july09); the alternatives above bracket the mapping ambiguity "
        f"({min(v['volume_km3'] for v in [stats, *alternatives.values()]):.3f}"
        f"-{max(v['volume_km3'] for v in [stats, *alternatives.values()]):.3f}"
        " km3), which Ian Hamling has been asked to resolve.",
        "",
        "Comparability caveats: the delivered historical fields were "
        "filtered below local uncertainty per the papers (the recovered "
        "0.5-m floor is specific to the may10/apr10 pair), while the "
        "2026 totals keep all fitted opening. As a uniform-floor "
        "sensitivity — not a matched reconstruction — a 0.5-m floor "
        "gives "
        f"{stats['volume_at_0p5m_floor_km3']:.3f} km3 for the primary "
        f"cumulative and {fit_stats['volume_at_0p5m_floor_km3']:.3f} km3 "
        "for 2026. The 2005 model "
        "carries "
        f"{epochs['sept05']['volume_km3'] - epochs['sept05']['volume_above_10km_km3']:.2f}"
        " km3 below 10 km depth (a deep lobe peaking at 13-20 km in its "
        f"northern half); the 2026 model plane ends at {bottom:g} km, so "
        "like-for-like volume comparisons should quote the above-10-km "
        "column. See also the notes block in the JSON and LIMITATIONS.md. "
        "Exact event dates are not part of the delivery; take them from "
        "the Hamling et al. (2009, 2010) catalogue when captioning.",
        "",
    ]
    md_path = RESULTS / f"uperp_three_epoch{suffix}.md"
    md_path.write_text("\n".join(lines))

    print(f"wrote {figure}, {json_path.name}, {md_path.name}")
    print(f"registration: {registration['columns_within_plane']}/75 columns "
          f"inside the plane span, s {registration['span_km'][0]:.2f}-"
          f"{registration['span_km'][1]:.2f} km, off-plane median "
          f"{registration['median_offplane_m']:.0f} m / max "
          f"{registration['max_offplane_m']:.0f} m")
    print(f"cumulative: {len(cum_tags)} files summed; excluded "
          f"{list(PRIMARY_EXCLUDED)} per the recorded pair diagnostics")
    print(f"volumes (km3): 2005 {epochs['sept05']['volume_km3']:.2f} "
          f"({epochs['sept05']['volume_above_10km_km3']:.2f} above 10 km), "
          f"2006-2010 {stats['volume_km3']:.2f} "
          f"({stats['volume_above_10km_km3']:.2f} above 10 km), "
          f"2026 {fit_stats['volume_km3']:.3f}")


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
