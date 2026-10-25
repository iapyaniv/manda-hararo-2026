"""Three-epoch opening comparison from the Hamling dike files.

Ian Hamling (GNS) provided `dyke_files.zip` (via Tim Wright,
2026-09-16), whose `opening_files/` holds opening
distributions of the September 2005 megadike and fifteen later
(2006-2010) opening files on the ~74 km Hamling computational grid
(Hamling et al. 2009, 2010; the published sequence has 13 post-megadike
intrusions). His 2026-09-17 update (tracked in
`data/ihamling/2026-09-17/`) resolved the open
questions: the canonical 13 post-2005 files keep july09 and may10 and
drop the defunct earlier versions june09 and apr10 (see
`PRIMARY_EXCLUDED`); the 20-km sept05 file is the raw inversion whose
deep opening is ill-constrained (his `model_slip_sept05.txt`
is the updated 10-km model); the comparison figure should show the
upper 10 km only; and, per his layout note of the same day, the
updated model's dyke section is lines 171-920 of its 931 values (10
one-km depth rows per column at the same 75 horizontal positions; the
untrimmed ends are fault sections), so it is decoded on the shared
columns. This module decodes the delivery, registers the grid
against the Wright et al. (2006) plane the 2026 fit uses
(mh26.wright2006), and renders the three-epoch U_perp comparison the
plane was built for: September 2005, the cumulative 2006-2010
sequence, and July 2026 on one along-plane axis.

Decoded format (verified by the structural checks in `load_openings`):

- `plot_3d_dyke_only.dat`: 1500 rows x 9 columns, one fault patch per row:
  easting (m), northing (m) [UTM 37N], strike (deg), dip (90), rake-slot
  (180), example opening (m, unused here), along-strike patch length (km),
  patch top (km), patch bottom (km). Rows are depth-fastest: 20 one-km
  depth rows (0-20 km) per along-strike column, 75 columns ordered SE->NW,
  in six constant-strike segments (329.7-355.7 deg, ~0.94-1.14 km column
  widths, ~73.9 km trace).
- `dyke_only_<epoch>_ex.txt`: 1500 opening values in metres, one per patch,
  in the `plot_3d_dyke_only.dat` row order. Ian's own plotting script
  (`doplot3d.m`, provided 2026-09-17) indexes an opening file by the
  geometry file's row number, confirming the row-order reading that was
  first adopted on the equal row count and the physically coherent
  profiles it produces (smooth depth profiles; the late small dikes
  vanish below ~8 km).
- `t300_*_wrp_dyke.grd`: GMT classic-netCDF grids (1024x1024, 90 m, UTM
  37N) of wrapped Envisat track-300 LOS in centimetres (z spans one fringe,
  0-2.81 cm = lambda/2); context only, not used here.
- `2026-09-17/model_slip_sept05.txt`: 931 values; lines 171-920 are the
  dyke section - 10 one-km depth rows (0-10 km) per along-strike column,
  depth-fastest, at the plot_3d_dyke_only.dat horizontal positions - and
  the untrimmed ends are fault sections (Ian Hamling, 2026-09-17;
  verified in `load_updated_sept05`).

The openings are Envisat-era inversions exactly as provided: no
uncertainties were delivered, and their regularization and greens
functions differ from the 2026 fit's. Every post-2005 file carries a
positive floor of 0.29-0.72 m, consistent with the papers' discard of
opening below the local 2-sigma uncertainty (typically 0.3-0.4 m). The
sept05 file is unfiltered (values down to 4e-4 m), integrates to
3.55 km3 against the published ~2.5 km3, and holds deep lobes (up to
12.9 m at 16-17 km); Ian confirmed (2026-09-17) it is the raw
inversion regridded to 20 km with the ill-constrained deep opening
kept, so the figure shows the upper 10 km and comparisons quote the
above-10-km column. The September 2005 panel draws his updated 10-km
model (`model_slip_sept05.txt`): its decoded dyke section integrates
to 2.23 km3 (maximum 16.1 m, all above 10 km), near the published
~2.5 km3, and correlates at r = 0.69 with the raw field's upper 10 km
(two inversions of the same event); the raw field's numbers stay in
the report table. The comparison is of delivered models on a shared
axis, not a joint inversion.

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
import matplotlib.path as mpath
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from .distributed import _display_source_path, load_geometry_source
from .params import REPO, RESULTS
from .wright2006 import _chain_frames, _project_s

DATA_DIR = REPO / "data" / "ihamling" / "opening_files"
GEOMETRY_FILE = "plot_3d_dyke_only.dat"
SOURCE_URL = "https://sftp.gns.cri.nz/pub/ihamling/dyke_files.zip"
RECEIVED = ("provided by Ian Hamling (GNS) via Tim Wright, 2026-09-16")

N_COLUMNS = 75
N_DEPTH = 20

# Filename month tags in event order. The delivery carries no exact dates;
# event dates must come from the Hamling et al. (2009, 2010) catalogue when
# captioning, so only the tags are recorded here.
EPOCH_TAGS = ("sept05", "june06", "july06", "sept06", "dec06", "jan07",
              "aug07", "nov07", "apr08", "july08", "oct08", "feb09",
              "june09", "july09", "apr10", "may10")

# The delivery holds 15 post-2005 files against the 13 published
# post-megadike intrusions. Ian Hamling's 2026-09-17 update
# lists the canonical 13 files and calls june09 and apr10 "defunct"
# earlier versions: the set keeps july09 (the re-estimate of the June
# 2009 event) and may10 (apr10 with every cell below 0.5 m set to
# zero — for consistency with the other, similarly floored files).
# `version_diagnostics` still re-derives the two version pairs from the
# data, and `compute_report` requires them to match DUPLICATE_PAIRS
# with exactly one member of each pair excluded, so a changed delivery
# cannot silently keep this selection.
PRIMARY_EXCLUDED = ("june09", "apr10")
DUPLICATE_PAIRS = (frozenset({"apr10", "may10"}),
                   frozenset({"june09", "july09"}))

# Ian Hamling's updated September 2005 model (provided 2026-09-17):
# opening limited to 10 km, consistent with the Wright et
# al. (2006) original. His layout note of the same day: the dyke
# section is lines 171-920 (1-indexed) of the 931 values - 10 one-km
# depth rows per column, depth-fastest, at the plot_3d_dyke_only.dat
# horizontal positions - and the untrimmed ends are fault sections.
UPDATED_SEPT05 = (REPO / "data" / "ihamling" / "2026-09-17"
                  / "model_slip_sept05.txt")
UPDATED_N_VALUES = 931
UPDATED_DYKE_SLICE = slice(170, 920)  # lines 171-920, 1-indexed
UPDATED_N_DEPTH = 10

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


def load_updated_sept05(delivery: dict) -> tuple[np.ndarray, dict]:
    """Decode the updated 10-km September 2005 model; raise on mismatch.

    The adopted reading (depth-fastest, geometry-file column order) must
    correlate with the raw regridded field's upper 10 km — two inversions
    of the same event, r = 0.69, while the transposed reading gives
    r = 0.0 — so a changed file or a wrong orientation fails loudly
    instead of drawing a scrambled field.
    """
    values = np.loadtxt(UPDATED_SEPT05)
    if values.shape != (UPDATED_N_VALUES,) or values.min() < 0.0:
        raise ValueError(f"{UPDATED_SEPT05}: expected {UPDATED_N_VALUES} "
                         f"non-negative values, got shape {values.shape}")
    dyke = values[UPDATED_DYKE_SLICE]
    field = dyke.reshape(N_COLUMNS, UPDATED_N_DEPTH)
    raw_upper = delivery["openings_m"]["sept05"][:, :UPDATED_N_DEPTH]
    adopted = float(np.corrcoef(field.ravel(), raw_upper.ravel())[0, 1])
    transposed = float(np.corrcoef(
        dyke.reshape(UPDATED_N_DEPTH, N_COLUMNS).T.ravel(),
        raw_upper.ravel())[0, 1])
    if not (adopted > 0.5 and adopted > transposed):
        raise ValueError(
            f"{UPDATED_SEPT05}: the depth-fastest dyke-section reading "
            f"correlates at r={adopted:.3f} with the raw upper-10-km field "
            f"(transposed reading r={transposed:.3f}); the file no longer "
            "matches the 2026-09-17 layout note")
    meta = {
        "file": _display_source_path(UPDATED_SEPT05),
        "sha256": _sha256(UPDATED_SEPT05),
        "n_values": int(values.size),
        "layout": "dyke section lines 171-920 (1-indexed): 10 one-km "
                  "depth rows (0-10 km) per column, depth-fastest, at "
                  "the plot_3d_dyke_only.dat horizontal positions; "
                  "lines 1-170 and 921-931 are untrimmed fault "
                  "sections, excluded (Ian Hamling, 2026-09-17)",
        "decode_correlation_vs_raw_upper10km": adopted,
        "trimmed_fault_value_sums_m": [
            float(values[:UPDATED_DYKE_SLICE.start].sum()),
            float(values[UPDATED_DYKE_SLICE.stop:].sum())],
    }
    return field, meta


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
    # Full precision: the renderers round exactly once at format time
    # (an intermediate round() shifted third-decimal displays by 1).
    return {
        "volume_km3": volume,
        "volume_above_10km_km3": above10,
        "max_opening_m": float(opening.max()),
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
                          "correlation": corr,
                          "exact_value_subset": exact,
                          "threshold_m": threshold})
    return pairs


def _fit_opening_stats(fit: dict, geometry: dict) -> dict:
    """Integrate the fitted opening on its own grid, split at 10 km.

    A depth row straddling 10 km contributes its fractional overlap. The
    0.5-m-floor volume is a uniform-floor sensitivity, mirroring the one
    floor recovered from the delivered may10/apr10 pair; the papers
    describe filtering by local uncertainty, but its application to the
    delivered files is unconfirmed, so this is not a matched
    reconstruction of their filtering.
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
    # Derive volume and max from the patch vector; the fit's summary
    # block is only cross-checked (a summary edited apart from its
    # vector once passed every validator).
    volume_km3 = float((opening * areas_km2).sum()) * 1e-3
    max_opening_m = float(opening.max())
    for key, derived in (("dike_volume_km3", volume_km3),
                         ("max_opening_m", max_opening_m)):
        if not np.isclose(fit["dike"][key], derived, rtol=1e-12, atol=0.0):
            raise ValueError(
                f"fit summary {key}={fit['dike'][key]} does not match the "
                f"value derived from the opening vector ({derived})")
    return {
        "volume_km3": volume_km3,
        "volume_above_10km_km3":
            float((opening * areas_km2 * above_frac[:, None]).sum()) * 1e-3,
        "volume_at_0p5m_floor_km3":
            float((floored * areas_km2).sum()) * 1e-3,
        "max_opening_m": max_opening_m,
        "plane_bottom_km": float(z_edges[-1]) / 1e3,
    }


def _path_pieces(path) -> list:
    """The connected pieces of a contour path (one open or closed line each)."""
    vertices, codes = path.vertices, path.codes
    if codes is None:
        return [vertices]
    starts = list(np.flatnonzero(codes == path.MOVETO)) + [len(vertices)]
    return [vertices[a:b] for a, b in zip(starts[:-1], starts[1:]) if b - a > 1]


def _flatness(xy: np.ndarray, half_window: int = 3) -> np.ndarray:
    """Per vertex, how far the neighbourhood departs from a straight line
    (smaller is flatter: an inline label breaks the contour least there);
    the end vertices get infinity."""
    n = len(xy)
    scores = np.full(n, np.inf)
    if n <= 2 * half_window + 1:
        scores[n // 2] = 0.0
        return scores
    for i in range(half_window, n - half_window):
        a, b = xy[i - half_window], xy[i + half_window]
        chord = b - a
        length = float(np.hypot(*chord))
        if length == 0.0:
            continue
        local = xy[i - half_window:i + half_window + 1] - a
        scores[i] = float(np.abs(chord[0] * local[:, 1] - chord[1] * local[:, 0]).sum()
                          / length)
    return scores


def _outside_box(xy: np.ndarray, box) -> list:
    """The parts of the polyline ``xy`` outside the rectangle ``box``
    (``(x0, y0, x1, y1)``, same coordinates), cut exactly at its edge."""
    x0, y0, x1, y1 = box
    p0, d = xy[:-1], np.diff(xy, axis=0)
    t_in, t_out = np.zeros(len(d)), np.ones(len(d))
    for k, lo, hi in ((0, x0, x1), (1, y0, y1)):
        with np.errstate(divide="ignore", invalid="ignore"):
            ta, tb = (lo - p0[:, k]) / d[:, k], (hi - p0[:, k]) / d[:, k]
        t_lo, t_hi = np.minimum(ta, tb), np.maximum(ta, tb)
        parallel = d[:, k] == 0.0
        within = (p0[:, k] >= lo) & (p0[:, k] <= hi)
        t_lo = np.where(parallel, np.where(within, -np.inf, np.inf), t_lo)
        t_hi = np.where(parallel, np.where(within, np.inf, -np.inf), t_hi)
        t_in, t_out = np.maximum(t_in, t_lo), np.minimum(t_out, t_hi)
    pieces, piece = [], [xy[0]]
    for i in range(len(d)):
        if t_in[i] >= t_out[i]:
            piece.append(xy[i + 1])
            continue
        if t_in[i] > 0.0:
            piece.append(p0[i] + t_in[i] * d[i])
        pieces.append(piece)
        piece = [p0[i] + t_out[i] * d[i]] if t_out[i] < 1.0 else []
        if piece:
            piece.append(xy[i + 1])
    pieces.append(piece)
    return [np.array(p) for p in pieces if len(p) > 1]


def _cut_path(path, transform, box):
    """``path`` (data coordinates) less its parts inside ``box`` (screen
    pixels)."""
    pieces = []
    for xy in _path_pieces(path):
        pieces += _outside_box(transform.transform(xy), box)
    return mpath.Path.make_compound_path(
        *(mpath.Path(transform.inverted().transform(xy)) for xy in pieces))


def _label_contours(ax, lines, *, fontsize=6.0, min_width_mm=4.0,
                    min_separation_mm=5.0, margin_km=0.6, gap_mm=0.5):
    """One upright label per contour level, on the widest line of that
    level, at the flattest stretch that stays clear of the panel edges and
    at least ``min_separation_mm`` from the labels already placed (lower
    levels are placed first); the line is cut ``gap_mm`` round the label's
    box. A line narrower than ``min_width_mm`` in either direction (the
    innermost ring round a peak) stays unlabelled. Call it once the axis
    limits are final: the labels follow the drawn lines."""
    px_per_mm = ax.figure.dpi / 25.4
    renderer = ax.figure.canvas.get_renderer()
    (x_lo, x_hi), (y_lo, y_hi) = sorted(ax.get_xlim()), sorted(ax.get_ylim())
    positions, levels, placed_px = [], [], []
    for level, path in enumerate(lines.get_paths()):
        pieces = _path_pieces(path)
        if not pieces:
            continue
        widest = max(pieces, key=lambda xy: np.ptp(xy, axis=0).max())
        if np.ptp(ax.transData.transform(widest), axis=0).min() < min_width_mm * px_per_mm:
            continue
        inside = widest[(widest[:, 0] > x_lo + margin_km) & (widest[:, 0] < x_hi - margin_km)
                        & (widest[:, 1] > y_lo + margin_km) & (widest[:, 1] < y_hi - margin_km)]
        if len(inside) < 2:
            continue
        inside_px = ax.transData.transform(inside)
        scores = _flatness(inside)
        for i in np.argsort(scores):
            if not np.isfinite(scores[i]):
                break
            if all(np.hypot(*(inside_px[i] - p)) >= min_separation_mm * px_per_mm
                   for p in placed_px):
                positions.append(tuple(inside[i]))
                levels.append(level)
                placed_px.append(inside_px[i])
                break
    if not positions:
        return
    paths = list(lines.get_paths())
    for level, text in zip(levels, lines.clabel(manual=positions, inline=False,
                                                fontsize=fontsize, fmt="%g")):
        text.set_rotation(0)
        text.set_path_effects([path_effects.Stroke(linewidth=1.6, foreground="white"),
                               path_effects.Normal()])
        # The digits have no descender: the box is trimmed by the font's.
        extent = text.get_window_extent(renderer)
        descent = renderer.get_text_width_height_descent(
            text.get_text(), text.get_fontproperties(), ismath=False)[2]
        gap = gap_mm * px_per_mm
        paths[level] = _cut_path(paths[level], ax.transData,
                                 (extent.x0 - gap, extent.y0 + descent - gap,
                                  extent.x1 + gap, extent.y1 + gap))
    lines.set_paths(paths)


def _panel(ax, s_centres, z_edges, opening, vmax, line_levels, title,
           note, title_size=10.5):
    """One filled U_perp cross-section in the mh26.wright2006 panel style.
    Returns the filled set and the contour-line set (None without lines)."""
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
    lines = None
    if levels:
        lines = ax.contour(s_plot, z_plot, padded, levels=levels, colors="0.15",
                           linewidths=0.55)
    ax.set_title(title, loc="left", fontsize=title_size, fontweight="bold",
                 pad=5 if note else 3)
    if note:
        ax.set_title(note, loc="right", fontsize=8, color="0.30", pad=6)
    ax.set_ylabel("Depth (km)")
    return filled, lines


# Contour lines of the three sections (metres): September 2005, the
# 2006-2010 sum, July 2026. The manuscript caption quotes them.
SECTION_CONTOURS = ((3, 6, 9, 12, 15), (1, 2, 4, 6), (0.5, 1.0, 1.5, 2.0))


def draw_epoch_sections(fig, axes, registration: dict, fit: dict,
                        geometry: dict, chain: str, epochs: dict,
                        cumulative: np.ndarray,
                        updated_field: np.ndarray,
                        common_vmax: float | None = None,
                        manuscript: bool = False, cax=None) -> tuple:
    """Draw the three U_perp sections onto ``axes``; returns their xlim.

    ``common_vmax`` replaces the per-section colour scales with that
    single scale and one shared colorbar, so the three amplitudes
    compare directly. The caller sets the x label: the July 2026
    section need not be the figure's bottom row. ``manuscript`` keeps
    only what a caption cannot carry: panel letters and titles, the
    sections, their contours, the plane-edge line and the grey band;
    the source notes, volume labels and edge labels are left out. The
    shared colorbar then goes into ``cax``.
    """
    s_ham = registration["s_km"][::-1]  # NW->SE for plotting
    # Displayed depth range is 0-10 km on every panel, per Ian Hamling
    # (2026-09-17): the 2005 panel is his updated model, which stops at
    # 10 km, and the raw regridded files keep ill-constrained deep
    # opening that is not shown.
    panels = [
        # The note names only the studies: "updated" would be
        # delivery-relative — a manuscript reader knows no earlier
        # model. The raw-vs-updated distinction stays in the md/JSON.
        ("September 2005", updated_field[::-1].T,
         s_ham, np.arange(11.0), 16.5, list(SECTION_CONTOURS[0]),
         "Hamling et al. (2009, 2010) model",
         epochs["sept05_updated"]),
        ("2006–2010 (cumulative, 13 events)", cumulative[::-1].T[:10],
         s_ham, np.arange(11.0), 7.0, list(SECTION_CONTOURS[1]),
         "sum of the 13 dike models of Hamling et al. (2009, 2010)",
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
                   list(SECTION_CONTOURS[2]), source_note, epochs["epoch_2026"]))

    xlim = (0.0, max(registration["span_km"][1] + 0.6,
                     registration["plane_span_km"][1]))
    se_2005 = registration["se_limit_2005_km"]
    if manuscript:
        letters = ("(a) September 2005", "(b) 2006–2010, thirteen dikes",
                   "(c) July 2026")
        panels = [(letter, *rest[1:6], "", rest[7])
                  for letter, rest in zip(letters, panels)]
    filled = None
    for ax, (title, opening, s_centres, z_ax, vmax, lines, note,
             stats) in zip(axes, panels):
        if common_vmax is not None:
            vmax = common_vmax
        filled, contour_lines = _panel(ax, s_centres, z_ax, opening, vmax, lines,
                                       title, note,
                                       title_size=8 if manuscript else 10.5)
        ax.set_xlim(*xlim)
        ax.set_ylim(0.0, 10.0)
        ax.invert_yaxis()
        ax.axvline(se_2005, color="0.25", lw=0.9, ls=(0, (4, 3)), zorder=4)
        if manuscript:
            if contour_lines is not None:
                _label_contours(ax, contour_lines)
            continue
        # The label quotes the displayed (above-10-km) field; the full
        # file integral stays in parentheses where it differs.
        if stats["volume_km3"] - stats["volume_above_10km_km3"] > 5e-3:
            label = (f"{stats['volume_above_10km_km3']:.2f} km$^3$ above "
                     f"10 km (file total {stats['volume_km3']:.2f}) "
                     f"· max {float(np.max(opening)):.1f} m")
        else:
            label = (f"{stats['volume_km3']:.2f} km$^3$ "
                     f"· max {float(np.max(opening)):.1f} m")
        ax.text(0.012, 0.06, label, transform=ax.transAxes, ha="left",
                va="bottom", fontsize=8.5,
                bbox={"facecolor": "white", "edgecolor": "none",
                      "alpha": 0.8, "pad": 1.6}, zorder=6)
        if common_vmax is None:
            ticks = np.arange(0.0, vmax + 0.01, max(1.0, round(vmax / 5)))
            if vmax <= 3.0:
                ticks = np.arange(0.0, vmax + 0.01, 0.5)
            colorbar = fig.colorbar(filled, ax=ax, label="U$_\\perp$ (m)",
                                    ticks=ticks, shrink=0.94, pad=0.015)
            colorbar.ax.tick_params(labelsize=8)
    if common_vmax is not None:
        ticks = np.arange(0.0, common_vmax + 0.01,
                          max(1.0, round(common_vmax / 5)))
        if cax is not None:
            colorbar = fig.colorbar(filled, cax=cax, label="Opening (m)",
                                    ticks=ticks)
            colorbar.ax.tick_params(labelsize=7, length=2, width=0.5)
            colorbar.outline.set_linewidth(0.5)
        else:
            colorbar = fig.colorbar(filled, ax=list(axes[:3]),
                                    label="U$_\\perp$ (m)", ticks=ticks,
                                    shrink=0.97, pad=0.015)
            colorbar.ax.tick_params(labelsize=8)

    ax_2026 = axes[2]
    if z_edges[-1] < 10.0:
        ax_2026.axhspan(z_edges[-1], 10.0, color="0.88", zorder=1)
        if not manuscript:
            ax_2026.text(0.5 * sum(xlim), 0.5 * (z_edges[-1] + 10.0),
                         f"below the 2026 model plane ({z_edges[-1]:g} km)",
                         ha="center", va="center", fontsize=8, color="0.35",
                         zorder=2)
    if xlim[1] > s_edges[-1]:
        # Beyond the plane's SE edge there is no 2026 model; grey it so it
        # does not read as zero opening.
        ax_2026.axvspan(s_edges[-1], xlim[1], color="0.88", zorder=1)
    if manuscript:
        return xlim
    # The dashed line is the Wright et al. (2006) model-plane SE edge.
    # The drawn (updated) 2005 model holds no opening southeast of it;
    # the raw regridded field holds ~0.22 km3 there (the payload's
    # volume_se_of_wright_plane_edge_km3 block records both).
    ax_2026.text(se_2005 + 0.4, 0.5, "Wright 2006 plane edge",
                 transform=ax_2026.get_xaxis_transform(), rotation=90,
                 ha="left", va="center", fontsize=7.5, color="0.20",
                 bbox={"facecolor": "white", "edgecolor": "none",
                       "alpha": 0.78, "pad": 1.2}, zorder=5)
    return xlim


def three_epoch_figure(delivery: dict, registration: dict, fit: dict,
                       geometry: dict, chain: str, epochs: dict,
                       cumulative: np.ndarray,
                       updated_field: np.ndarray) -> Path:
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.2), sharex=True,
                             constrained_layout=True)
    draw_epoch_sections(fig, axes, registration, fit, geometry, chain,
                        epochs, cumulative, updated_field)
    axes[2].set_xlabel("Distance along the 2005 dike plane, NW→SE (km)")

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


def compute_report(chain: str, fit_path: Path, data_dir: Path) -> tuple:
    """Rebuild the full report payload deterministically.

    mh26.validate calls this against the canonical fit and the tracked
    delivery and requires exact equality with the stored JSON, so every
    number in the report is re-derived, not just hash-checked.
    """
    data_dir = Path(data_dir).resolve()
    delivery = load_openings(data_dir)
    fit_path = Path(fit_path).resolve()
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
    pairs = {frozenset({d["contained"], d["container"]}) for d in diagnostics}
    if pairs != set(DUPLICATE_PAIRS):
        raise ValueError(
            f"data-derived version pairs {sorted(map(sorted, pairs))} no "
            f"longer match {sorted(map(sorted, DUPLICATE_PAIRS))}; the "
            "delivery changed — revisit PRIMARY_EXCLUDED (Ian Hamling's "
            "2026-09-17 mapping) before summing")
    for pair in DUPLICATE_PAIRS:
        if len(pair & set(PRIMARY_EXCLUDED)) != 1:
            raise ValueError(
                f"PRIMARY_EXCLUDED {PRIMARY_EXCLUDED} must drop exactly "
                f"one member of the version pair {sorted(pair)}")
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

    updated_field, updated_meta = load_updated_sept05(delivery)
    epochs["sept05_updated"] = _epoch_stats(updated_field, widths)
    raw05 = delivery["openings_m"]["sept05"]
    se_mask = registration["s_km"] > registration["se_limit_2005_km"]
    se_volumes = {
        "sept05_raw": float(
            (raw05[se_mask] * widths[se_mask, None]).sum()) * 1e-3,
        "sept05_raw_above_10km": float(
            (raw05[se_mask, :10] * widths[se_mask, None]).sum()) * 1e-3,
        "sept05_updated": float(
            (updated_field[se_mask] * widths[se_mask, None]).sum()) * 1e-3,
    }
    # Why the drawn 2005 opening sits further NW than the raw field's
    # upper 10 km: the raw model's below-10-km opening is NW-weighted,
    # and the 10-km inversion must carry that moment in shallow patches.
    centroids = {}
    for key, field in (("raw_upper_10km", raw05[:, :10]),
                       ("raw_below_10km", raw05[:, 10:]),
                       ("updated_model", updated_field)):
        col_volumes = field.sum(axis=1) * widths
        centroids[key] = float((registration["s_km"] * col_volumes).sum()
                               / col_volumes.sum())
    deep_volumes = raw05[:, 10:].sum(axis=1) * widths
    along_strike = {
        "volume_centroid_km_from_nw": centroids,
        "raw_below_10km_fraction_nw_of_25km": float(
            deep_volumes[registration["s_km"] < 25.0].sum()
            / deep_volumes.sum()),
    }
    updated_sept05 = {
        **updated_meta,
        **epochs["sept05_updated"],
        "volume_km3_at_1km2_patches": float(updated_field.sum()) * 1e-3,
        "note": "Ian Hamling's updated 10-km September 2005 model "
                "(provided 2026-09-17), consistent with the Wright et al. "
                "(2006) original; decoded per his layout note of the "
                "same day and drawn as the September 2005 panel — the "
                "trimmed ends are fault sections, so the whole-file sum "
                "(2.478 km3 at 1-km2 patches) overstates the dyke",
    }

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
            "mapping": "canonical 13 post-2005 files per Ian Hamling "
                       "(2026-09-17; june09 and apr10 are "
                       "defunct earlier versions)",
            "pair_diagnostics": diagnostics,
            **epochs["cumulative_2006_2010"],
            "alternatives": {
                "sum_of_all_delivered_files": _epoch_stats(
                    sum(delivery["openings_m"][t] for t in EPOCH_TAGS[1:]),
                    widths),
                "superseded_primary_drop_july09_may10": _epoch_stats(
                    sum(delivery["openings_m"][t] for t in EPOCH_TAGS[1:]
                        if t not in ("july09", "may10")), widths),
                "untrimmed_apr10_for_may10": _epoch_stats(
                    sum(delivery["openings_m"][t] for t in EPOCH_TAGS[1:]
                        if t not in ("june09", "may10")), widths),
            },
        },
        "sept05_updated_model": updated_sept05,
        "volume_se_of_wright_plane_edge_km3": se_volumes,
        "sept05_along_strike": along_strike,
        "epoch_2026": {
            "reference_fit": _display_source_path(fit_path),
            "reference_fit_sha256": _sha256(fit_path),
            "geometry_sha256": fit_prov["geometry_sha256"],
            **epochs["epoch_2026"],
        },
        "notes": [
            "openings are Envisat-era inversions exactly as provided "
            "(no uncertainties; different regularization and greens "
            "functions from the 2026 fit)",
            "the epoch-file row order follows plot_3d_dyke_only.dat, "
            "confirmed by Ian Hamling's own plotting script doplot3d.m "
            "(provided 2026-09-17), which indexes an opening file by the "
            "geometry file's row number",
            "the delivered sept05 file is the raw inversion regridded "
            "to 20 km (Ian Hamling, 2026-09-17): it integrates to "
            "3.55 km3 against the published ~2.5 km3 and its deep "
            "opening (up to 12.9 m at 16-17 km; 1.74 km3 below 10 km) "
            "is ill-constrained, so the figure shows the upper 10 km "
            "and like-for-like comparisons quote the above-10-km "
            "column; his updated 10-km model (sept05_updated_model) is "
            "decoded per his 2026-09-17 layout note — dyke section "
            "lines 171-920, 10 one-km rows per column at the same 75 "
            "horizontal positions, the trimmed ends fault sections — "
            "and drawn as the September 2005 panel",
            "the updated model shifts opening toward the NW relative "
            "to the raw field's upper 10 km because the raw "
            "below-10-km opening is NW-weighted and the 10-km "
            "inversion carries that moment in shallow patches; "
            "sept05_along_strike records the volume centroids",
            "the delivery has 15 post-2005 files against the 13 "
            "published post-megadike intrusions; Ian Hamling's "
            "2026-09-17 update lists the canonical 13 (keep july09 and "
            "may10, drop the defunct june09 and apr10), matching the "
            "data-derived version pairs (may10 is apr10 with cells "
            "below 0.5 m removed; july09's positive cells all lie "
            "inside june09's at r=0.994, all other pairs |r|<=0.44); "
            "the alternatives record the superseded selection and the "
            "untrimmed-apr10 variant",
            "every post-2005 file carries a positive floor of "
            "0.29-0.72 m, consistent with the papers' discard of "
            "opening below the local 2-sigma uncertainty (typically "
            "0.3-0.4 m; only may10's 0.5-m cut is proven exactly, "
            "against apr10), while sept05 is unfiltered (cells down to "
            "4e-4 m) and the 2026 totals keep all fitted opening; the "
            "volume_at_0p5m_floor_km3 values apply one uniform 0.5-m "
            "floor to the primary cumulative and to 2026 as a "
            "sensitivity, not as a matched reconstruction of any "
            "delivered filtering",
            "the shared s axis compares openings by along-axis position; "
            "the two traces differ by up to "
            f"{registration['max_offplane_m']:.0f} m mid-plane",
        ],
        "figure": f"uperp_three_epoch{suffix}.png",
    }
    context = {"delivery": delivery, "registration": registration,
               "fit": fit, "geometry": geometry, "epochs": epochs,
               "cum": cum, "cum_tags": cum_tags,
               "updated_field": updated_field}
    return payload, context


def run(chain: str, fit_path: Path | None, data_dir: Path) -> None:
    fit_path = Path(fit_path or _find_fit(chain)).resolve()
    payload, context = compute_report(chain, fit_path, data_dir)
    delivery = context["delivery"]
    registration = context["registration"]
    fit, geometry = context["fit"], context["geometry"]
    epochs, cum = context["epochs"], context["cum"]
    cum_tags = context["cum_tags"]

    figure = three_epoch_figure(delivery, registration, fit, geometry,
                                chain, epochs, cum,
                                context["updated_field"])

    suffix = "_gslc" if chain == "gslc" else ""
    json_path = RESULTS / f"uperp_three_epoch{suffix}.json"
    json_path.write_text(json.dumps(payload, indent=1, sort_keys=True))

    lines = [
        "# Three-epoch opening comparison (Hamling dike files)",
        "",
        f"NISAR chain of the 2026 fit: {chain} "
        f"({payload['epoch_2026']['reference_fit']}).",
        "",
        "Openings of the September 2005 megadike and the fifteen "
        "post-2005 opening files are the inversions provided by "
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
        if tag == "sept05":
            upd_stats = epochs["sept05_updated"]
            lines.append(f"| sept05 — updated 10-km model (drawn in the "
                         f"figure) | {upd_stats['volume_km3']:.3f} | "
                         f"{upd_stats['volume_above_10km_km3']:.3f} | "
                         f"{upd_stats['max_opening_m']:.2f} |")
    stats = epochs["cumulative_2006_2010"]
    alternatives = payload["cumulative_2006_2010"]["alternatives"]
    lines.append(f"| 2006-2010 cumulative ({len(cum_tags)} files, primary) "
                 f"| {stats['volume_km3']:.3f} | "
                 f"{stats['volume_above_10km_km3']:.3f} | "
                 f"{stats['max_opening_m']:.2f} |")
    for label, alt in (
            (f"all {len(EPOCH_TAGS) - 1} delivered files",
             alternatives["sum_of_all_delivered_files"]),
            ("superseded primary (kept june09 + apr10)",
             alternatives["superseded_primary_drop_july09_may10"]),
            ("untrimmed apr10 in place of may10",
             alternatives["untrimmed_apr10_for_may10"])):
        lines.append(f"| — alternative: {label} | "
                     f"{alt['volume_km3']:.3f} | "
                     f"{alt['volume_above_10km_km3']:.3f} | "
                     f"{alt['max_opening_m']:.2f} |")
    fit_stats = epochs["epoch_2026"]
    lines.append(f"| July 2026 ({chain}) | {fit_stats['volume_km3']:.3f} | "
                 f"{fit_stats['volume_above_10km_km3']:.3f} | "
                 f"{fit_stats['max_opening_m']:.2f} |")
    bottom = fit_stats["plane_bottom_km"]
    updated = payload["sept05_updated_model"]
    lines += [
        "",
        "The delivery has 15 post-2005 files against the 13 published "
        "post-megadike intrusions. Ian Hamling's update of "
        "17 September 2026 lists the canonical 13 files and calls "
        "june09 and apr10 defunct earlier versions, so the primary "
        "cumulative keeps july09 (the re-estimate of the June 2009 "
        "event) and may10 (apr10 with every cell below 0.5 m set to "
        "zero, 38% of that field's volume — kept for consistency with "
        "the other floored files). The data-derived pair diagnostics "
        "match that identification: july09's positive cells all lie "
        "inside june09's at correlation 0.994 while no other pair "
        "exceeds |r| = 0.44. The alternatives above record the "
        "superseded selection and the untrimmed-apr10 variant.",
        "",
        "The row order of the opening files follows "
        "plot_3d_dyke_only.dat; Ian Hamling's own plotting script "
        "(doplot3d.m, provided 17 September 2026, tracked in "
        "data/ihamling/2026-09-17/) indexes an opening file by "
        "the geometry file's row number, confirming that reading.",
        "",
        "Comparability caveats: every post-2005 file carries a "
        "positive floor of 0.29-0.72 m, consistent with the papers' "
        "discard of opening below the local 2-sigma uncertainty "
        "(typically 0.3-0.4 m), while the delivered sept05 file is the "
        "raw inversion regridded to 20 km (Ian Hamling, 17 September "
        "2026): it is unfiltered (cells down to 4e-4 m), integrates to "
        f"{epochs['sept05']['volume_km3']:.2f} km3 against the "
        "published ~2.5 km3, and carries "
        f"{epochs['sept05']['volume_km3'] - epochs['sept05']['volume_above_10km_km3']:.2f}"
        " km3 of ill-constrained opening below 10 km (deep lobes "
        f"peaking at 13-20 km, up to "
        f"{epochs['sept05']['max_opening_m']:.1f} m). The figure "
        "therefore shows the upper 10 km of every panel, and its "
        "September 2005 panel draws his updated 10-km model "
        f"(`{Path(updated['file']).name}`, provided 17 September 2026), "
        "decoded per his layout note of the same day: the dyke section "
        "is lines 171-920 of the 931 values — 10 one-km depth rows per "
        "column, depth-fastest, at the same 75 horizontal positions — "
        "and the trimmed ends are fault sections. The decoded field "
        "correlates at r = "
        f"{updated['decode_correlation_vs_raw_upper10km']:.2f} with the "
        "raw field's upper 10 km (two inversions of the same event) "
        f"and integrates to {updated['volume_km3']:.3f} km3 with the "
        "geometry file's column widths (maximum "
        f"{updated['max_opening_m']:.2f} m, all above 10 km), near the "
        "published ~2.5 km3 — the earlier whole-file total of "
        "2.478 km3 at 1-km2 patches included the fault values now "
        "trimmed. Relative to the raw field's upper 10 km, the updated "
        "model shifts opening toward the NW (along-strike volume "
        "centroid "
        f"{payload['sept05_along_strike']['volume_centroid_km_from_nw']['raw_upper_10km']:.1f}"
        " -> "
        f"{payload['sept05_along_strike']['volume_centroid_km_from_nw']['updated_model']:.1f}"
        " km from the NW end): the raw model's below-10-km opening is "
        "itself NW-weighted ("
        f"{payload['sept05_along_strike']['raw_below_10km_fraction_nw_of_25km']:.0%}"
        " of it NW of s = 25 km, centroid "
        f"{payload['sept05_along_strike']['volume_centroid_km_from_nw']['raw_below_10km']:.1f}"
        " km), and limited to 10 km the inversion carries that moment "
        "in shallow NW patches instead — hence also the larger drawn "
        "maximum. The updated model holds "
        f"{payload['volume_se_of_wright_plane_edge_km3']['sept05_updated']:.3f}"
        " km3 southeast of the Wright 2006 plane edge (the raw field "
        "holds "
        f"{payload['volume_se_of_wright_plane_edge_km3']['sept05_raw']:.2f}"
        " km3 there, "
        f"{payload['volume_se_of_wright_plane_edge_km3']['sept05_raw_above_10km']:.2f}"
        " km3 of it above 10 km), so the dashed edge line also bounds "
        "the drawn 2005 opening. "
        "Like-for-like volume comparisons should quote the "
        "above-10-km column. The 2026 totals keep "
        "all fitted opening. As a uniform-floor sensitivity — not a "
        "matched reconstruction — a 0.5-m floor gives "
        f"{stats['volume_at_0p5m_floor_km3']:.3f} km3 for the primary "
        f"cumulative and {fit_stats['volume_at_0p5m_floor_km3']:.3f} km3 "
        "for 2026. "
        f"The 2026 model plane ends at {bottom:g} km. "
        "See also the notes block in the JSON and LIMITATIONS.md. "
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
          f"{list(PRIMARY_EXCLUDED)} per Ian Hamling's 2026-09-17 mapping")
    print(f"volumes (km3): 2005 raw {epochs['sept05']['volume_km3']:.2f} "
          f"({epochs['sept05']['volume_above_10km_km3']:.2f} above 10 km), "
          f"2005 updated {epochs['sept05_updated']['volume_km3']:.2f} "
          "(drawn), "
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
