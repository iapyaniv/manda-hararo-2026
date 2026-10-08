"""Dike-plane geometry from the Wright et al. (2006) source model.

Builds a complete geometry contract (mh26.distributed --geometry) from the
collaborator-provided file `data/complex_geom2_revised_test.dat`, so the 2026
opening distribution can be inverted on the same plane as the 2005 intrusion
and the 2006-2010 sequence (Wright et al. 2006; Hamling et al. 2009, 2010).

The file is an ASCII table, one source per row, 12 columns, UTM 37N
kilometres. The column meaning is not documented in the file; it was inferred
from the values and verified by the segment-chaining check below:

    1 strike (deg)   2 dip (deg)      3 rake (deg)    4 easting (km)
    5 northing (km)  6 length (km)    7 top (km)      8 bottom (km)
    9 n along-strike patches          10 n down-dip patches
    11 group tag     12 type (1 fault, 2 dike, 3 point source)

For type-3 point rows the shared top/bottom names do not apply: column 7 is
the source depth (5 km, consistent with the Wright et al. 2006 Mogi sources)
and column 8 is unused (0).

Coordinates are segment CENTRES: the five type-2 rows chain end-to-end into
one continuous 61 km trace once each centre is extended by +-length/2 along
its strike (junction residuals are below 0.1 m, recorded in the output).
Type-1 rows are the Dabbahu graben-bounding faults and type-3 rows the two
point sources of the 2005 model; both lie on the northern half of the plane,
are not part of the 2026 problem, and are stored as metadata only.

This is the Wright et al. (2006) source geometry. It is NOT the later
Hamling et al. (2009) computational grid (~75 km, 1 km x 1 km patches, whose
northern part follows this geometry); do not present results on this plane
as being on that grid.

The 2026 uniform-fit dike is aligned with the southernmost segment (strike
differs by ~4 deg; its endpoints lie 0.35 and 0.88 km off the segment's
line) and overruns the segment's SE end by ~3 km, so the plane is extended
SE along that segment's strike (--extend-se) and the extension is recorded
in the output. The along-plane coordinate s starts at 0 at the NW end of the
ORIGINAL plane and increases SE. The delivered 2005/2006-2010 opening
distributions project monotonically onto this axis (mh26.hamling), so s
needs no remap.

The plane is segmented, so no single straight strike-frame reproduces it.
The scalar `strike`/`trace_centre` contract fields are set to the chord of
the chain (origin at s=0, azimuth toward increasing s); reconstructing patch
positions from them is approximate, with the maximum chord error recorded in
the provenance. Consumers that need positions must use the saved patch
arrays. mh26.stress refuses fits on this geometry for that reason, and the
grafted fault's own straight-frame axis metadata is preserved under
`wright2006_provenance.fault_frame`.

By default the trace-supported east graben fault of the primary 2026 fit is
grafted in unchanged (absolute patch positions from
results/distributed_fit_mapped.json), so shallow graben slip is not forced
into spurious shallow opening; --no-faults writes a dike-only contract.

Build the contract:   python -m mh26.wright2006
Then invert:          python -m mh26.distributed --geometry results/geometry_wright2006.json
U_perp panel:         python -m mh26.wright2006 --panel results/distributed_fit_geometry_<hash>.json

A --gslc fit keeps a _gslc panel name, so the two chains cannot overwrite
each other. The publication-style version of this panel is
``mh26.summary_figures --only wright_uperp``.
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

from .distributed import (
    _display_source_path,
    _patch_areas,
    _strike_vectors,
    geometry_sha256,
    geometry_to_json,
    load_geometry_source,
)
from .params import REPO, RESULTS

DAT_DEFAULT = REPO / "data" / "complex_geom2_revised_test.dat"
MAPPED_FIT = RESULTS / "distributed_fit_mapped.json"
OUT_GEOMETRY = RESULTS / "geometry_wright2006.json"

# Junctions in the source model are not exactly closed; residuals beyond this
# would mean the column interpretation is wrong, not that the model is sloppy.
MAX_JUNCTION_GAP_M = 400.0

COLUMNS = ("strike_deg", "dip_deg", "rake_deg", "east_km", "north_km",
           "length_km", "top_km", "bottom_km", "n_along", "n_downdip",
           "group", "type")


def parse_dat(path: Path) -> dict:
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        values = [float(v) for v in line.split()]
        if len(values) != len(COLUMNS):
            raise ValueError(f"{path}: expected {len(COLUMNS)} columns, "
                             f"got {len(values)}")
        rows.append(dict(zip(COLUMNS, values)))
    by_type = {1.0: "fault", 2.0: "dike", 3.0: "point"}
    out = {"fault": [], "dike": [], "point": []}
    for row in rows:
        kind = by_type.get(row["type"])
        if kind is None:
            raise ValueError(f"{path}: unknown source type {row['type']}")
        out[kind].append(row)
    if not out["dike"]:
        raise ValueError(f"{path}: no dike segments (type 2) found")
    return out


def _point_record(row: dict) -> dict:
    # For type-3 rows the shared top/bottom names would record an inverted
    # depth interval (top 5, bottom 0); column 7 is the source depth and
    # column 8 is unused.
    record = {k: row[k] for k in COLUMNS
              if k not in ("top_km", "bottom_km")}
    record["depth_km"] = row["top_km"]
    record["column8_unused"] = row["bottom_km"]
    return record


def _segment_ends(row: dict) -> tuple[np.ndarray, np.ndarray]:
    """(NW end, SE end) in metres (north, east)."""
    along, _ = _strike_vectors(row["strike_deg"])
    centre = 1e3 * np.array([row["north_km"], row["east_km"]])
    half = 500.0 * row["length_km"]
    a, b = centre + along * half, centre - along * half
    return (a, b) if a[0] >= b[0] else (b, a)


def chain_segments(dike_rows: list[dict]) -> tuple[list[dict], list[float]]:
    """Order segments NW->SE and measure the junction closure residuals."""
    segments = sorted(dike_rows, key=lambda r: -r["north_km"])
    residuals = []
    for prev, nxt in zip(segments[:-1], segments[1:]):
        gap = float(np.linalg.norm(_segment_ends(prev)[1]
                                   - _segment_ends(nxt)[0]))
        residuals.append(gap)
        if gap > MAX_JUNCTION_GAP_M:
            raise ValueError(
                f"segment junction gap {gap:.0f} m exceeds "
                f"{MAX_JUNCTION_GAP_M:.0f} m; column interpretation is suspect")
    return segments, residuals


def build_contract(dat_path: Path, extend_se: int,
                   with_faults: bool) -> dict:
    sources = parse_dat(dat_path)
    segments, junction_gaps = chain_segments(sources["dike"])

    tops = {row["top_km"] for row in segments}
    bottoms = {row["bottom_km"] for row in segments}
    n_downs = {row["n_downdip"] for row in segments}
    if len(tops) != 1 or len(bottoms) != 1 or len(n_downs) != 1:
        raise ValueError("dike segments disagree on depth extent or rows; "
                         "a single (z, s) grid no longer applies")
    dips = {row["dip_deg"] for row in segments}
    if dips != {90.0}:
        raise ValueError(
            f"dike segment dips {sorted(dips)} are not supported: patch "
            "construction assumes a vertical plane, so a non-90 dip in a "
            "custom --dat table would be silently discarded")
    z_top = 1e3 * tops.pop()
    z_bot = 1e3 * bottoms.pop()
    nz = int(n_downs.pop())
    z_edges = np.linspace(z_top, z_bot, nz + 1)

    # Column descriptors NW->SE: (centre_n, centre_e, strike, width_m).
    columns = []
    for row in segments:
        along, _ = _strike_vectors(row["strike_deg"])
        nw_end, _ = _segment_ends(row)
        n_along = int(row["n_along"])
        width = 1e3 * row["length_km"] / n_along
        # The NW end may be either +along/2 or -along/2 from the centre;
        # step from the NW end toward the SE end regardless.
        step = -along if along[0] >= 0 else along
        for j in range(n_along):
            centre = nw_end + step * width * (j + 0.5)
            columns.append((centre[0], centre[1], row["strike_deg"], width))
    se_row = segments[-1]
    along_se, _ = _strike_vectors(se_row["strike_deg"])
    step_se = -along_se if along_se[0] >= 0 else along_se
    se_end = _segment_ends(se_row)[1]
    ext_width = 1e3 * se_row["length_km"] / int(se_row["n_along"])
    for j in range(extend_se):
        centre = se_end + step_se * ext_width * (j + 0.5)
        columns.append((centre[0], centre[1], se_row["strike_deg"], ext_width))

    s_edges = np.concatenate([[0.0], np.cumsum([w for *_, w in columns])])
    ns = len(columns)

    patches, disl, kinds = [], [], []
    for z0, z1 in zip(z_edges[:-1], z_edges[1:]):
        for c_n, c_e, strike, width in columns:
            patches.append([c_n, c_e, 0.5 * (z0 + z1), strike, 90.0,
                            -width / 2, width / 2,
                            -(z1 - z0) / 2, (z1 - z0) / 2])
            disl.append([0.0, 0.0, 1.0])
            kinds.append("dike")

    # Scalar frame: the chord of the chained trace, origin at s=0 and azimuth
    # toward increasing s, so `trace_centre + along(strike) * s` is at least
    # directionally and origin-consistent with the s axis. It remains a chord
    # approximation of a bent plane; the maximum error against the true patch
    # centres is recorded in the provenance, and consumers that need positions
    # must use the patch arrays.
    origin = _segment_ends(segments[0])[0]
    chord = _segment_ends(segments[-1])[1] - origin
    chord_strike = float(np.degrees(np.arctan2(chord[1], chord[0])) % 360.0)
    along_chord, _ = _strike_vectors(chord_strike)
    s_centres = 0.5 * (s_edges[:-1] + s_edges[1:])
    column_centres = np.array([[c_n, c_e] for c_n, c_e, *_ in columns])
    chord_error = float(np.linalg.norm(
        column_centres - (origin + np.outer(s_centres, along_chord)),
        axis=1).max())

    geometry = {
        "geometry_contract_version": 1,
        "dike_shape": (nz, ns),
        "dike_s_edges_m": s_edges,
        "dike_z_edges_m": z_edges,
        "dike_s_centres": s_centres,
        "dike_z_centres": 0.5 * (z_edges[:-1] + z_edges[1:]),
        "strike": chord_strike,
        "trace_centre": (float(origin[1]), float(origin[0])),
    }

    fault_source_note = "none (dike-only Wright plane)"
    fault_frame = None
    if with_faults:
        mapped, mapped_prov = load_geometry_source(MAPPED_FIT)
        # The mapped fit's fault-axis metadata (column offsets, s axis) is
        # defined in ITS straight frame, which this contract's scalars no
        # longer describe; keep that frame so the metadata stays decodable.
        fault_frame = {
            "strike_deg": float(mapped["strike"]),
            "trace_centre_e_n_m": [float(v) for v in mapped["trace_centre"]],
        }
        sel = mapped["kinds"] != "dike"
        for row, d, kind in zip(mapped["patches"][sel], mapped["disl"][sel],
                                mapped["kinds"][sel]):
            patches.append(list(row))
            disl.append(list(d))
            kinds.append(str(kind))
        for key in ("fault_shape", "fault_s_edges_m", "fault_s_centres_m",
                    "fault_z_edges_m", "min_fault_patch_thickness_m",
                    "fault_offset", "fault_dip", "fault_trace_source",
                    "enabled_fault_sides", "fault_column_offsets_m",
                    "fault_masks"):
            geometry[key] = mapped[key]
        fault_source_note = (
            "east graben fault grafted unchanged from the primary "
            f"trace-supported fit ({mapped_prov['path']}, geometry "
            f"{mapped_prov['geometry_sha256'][:12]}); its along-strike axis "
            "metadata refers to the uniform-fit frame preserved in "
            "wright2006_provenance.fault_frame, not the s axis above")
    else:
        geometry.update({
            "fault_shape": (1, 1),
            "fault_s_edges_m": np.array([0.0, 1.0]),
            "fault_s_centres_m": np.array([0.5]),
            "fault_z_edges_m": np.array([0.0, 1.0]),
            "min_fault_patch_thickness_m": 0.0,
            "fault_offset": 0.0,
            "fault_dip": 0.0,
            "fault_trace_source": "none",
            "enabled_fault_sides": [],
            "fault_column_offsets_m": {"east": [0.0], "west": [0.0]},
            "fault_masks": {"east": np.zeros((1, 1), dtype=bool),
                            "west": np.zeros((1, 1), dtype=bool)},
        })

    geometry["patches"] = np.array(patches, dtype=np.float64)
    geometry["disl"] = np.array(disl, dtype=np.float64)
    geometry["kinds"] = np.array(kinds)
    geometry["areas"] = _patch_areas(geometry["patches"])

    geometry["wright2006_provenance"] = {
        "dat_path": _display_source_path(dat_path),
        "dat_sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "received": "provided by Juliet Biggs, 2026-08-26 (geometry of "
                    "Wright et al. 2006 and the Hamling et al. papers)",
        "column_interpretation": " | ".join(COLUMNS),
        "column_note": "top_km/bottom_km apply to fault and dike rows; for "
                       "point rows column 7 is the source depth and column 8 "
                       "is unused (recorded as depth_km/column8_unused)",
        "units": "UTM 37N kilometres, converted to metres here",
        "segments_nw_to_se": [
            {k: row[k] for k in COLUMNS} for row in segments],
        "junction_residuals_m": [round(g, 1) for g in junction_gaps],
        "s_axis": "s=0 at the NW end of the original plane, increasing SE; "
                  "the SE extension occupies s > original length",
        "straight_frame_max_error_m": round(chord_error, 1),
        "frame_note": "strike/trace_centre describe the chord of the "
                      "segmented trace (origin at s=0, azimuth toward "
                      "increasing s); position reconstruction from them is "
                      "approximate — use the patch arrays",
        "fault_frame": fault_frame,
        "extension_se": {
            "n_patches": extend_se,
            "patch_width_m": ext_width,
            "reason": "the 2026 dike is aligned with the southernmost "
                      "segment (~4 deg strike difference, endpoints 0.35 and "
                      "0.88 km off its line) and overruns its SE end; the "
                      "extension lets the data bound the opening instead of "
                      "the plane edge",
        },
        "unmodeled_rows": {
            "graben_faults_2005": [
                {k: row[k] for k in COLUMNS} for row in sources["fault"]],
            "point_sources_2005": [
                _point_record(row) for row in sources["point"]],
        },
        "fault_treatment": fault_source_note,
    }
    return geometry


def _chain_frames(segments: list[dict], ext_len_m: float) -> list[dict]:
    """Per-segment straight frames NW->SE; the last one carries the extension."""
    frames, s0 = [], 0.0
    for i, row in enumerate(segments):
        along, _ = _strike_vectors(row["strike_deg"])
        d = -along if along[0] >= 0 else along
        length = 1e3 * row["length_km"]
        if i == len(segments) - 1:
            length += ext_len_m
        frames.append({"s0": s0, "nw": _segment_ends(row)[0], "d": d,
                       "L": length})
        s0 += length
    return frames


def _project_s(frames: list[dict], point: np.ndarray) -> tuple[float, float]:
    """Signed along-chain coordinate of a point and its off-plane distance.

    The terminal segment lines extend to infinity, so a point beyond either
    end of the plane gets an s outside the plane edges instead of being
    clamped onto them; interior segments clamp to their own extent so joints
    resolve to the nearer side.
    """
    best = None
    for i, f in enumerate(frames):
        t = float((point - f["nw"]) @ f["d"])
        lo = -np.inf if i == 0 else 0.0
        hi = np.inf if i == len(frames) - 1 else f["L"]
        t = min(max(t, lo), hi)
        perp = float(np.linalg.norm(point - (f["nw"] + f["d"] * t)))
        if best is None or perp < best[1]:
            best = (f["s0"] + t, perp)
    return best


def coverage_report(geometry: dict) -> dict:
    """Along-plane span of the 2026 uniform-fit trace, in the s coordinate."""
    fit_path = RESULTS / "okada_fit.json"
    model = json.loads(fit_path.read_text())["values"]
    e_c, n_c, strike, length = model[0], model[1], model[2], model[3]
    along, _ = _strike_vectors(strike)
    prov = geometry["wright2006_provenance"]
    ext = prov["extension_se"]
    frames = _chain_frames(prov["segments_nw_to_se"],
                           ext["n_patches"] * ext["patch_width_m"])
    spans, offsets = [], []
    for sign in (+1.0, -1.0):
        s, perp = _project_s(frames,
                             np.array([n_c, e_c]) + sign * along * length / 2)
        spans.append(s)
        offsets.append(perp)
    s_lo, s_hi = sorted(spans)
    return {
        "uniform_trace_span_s_km": [s_lo / 1e3, s_hi / 1e3],
        "plane_span_s_km": [float(geometry["dike_s_edges_m"][0]) / 1e3,
                            float(geometry["dike_s_edges_m"][-1]) / 1e3],
        "margin_nw_km": (s_lo - float(geometry["dike_s_edges_m"][0])) / 1e3,
        "margin_se_km": (float(geometry["dike_s_edges_m"][-1]) - s_hi) / 1e3,
        "endpoint_offplane_offsets_m": [round(v, 1) for v in offsets],
    }


def write_contract(args: argparse.Namespace) -> None:
    if args.extend_se < 0:
        raise SystemExit("--extend-se must be >= 0")
    geometry = build_contract(Path(args.dat).resolve(), args.extend_se,
                              with_faults=not args.no_faults)
    coverage = coverage_report(geometry)
    geometry["wright2006_provenance"]["coverage_2026_uniform_fit"] = coverage
    if coverage["margin_se_km"] < 0 or coverage["margin_nw_km"] < 0:
        raise SystemExit(
            f"2026 trace (s {coverage['uniform_trace_span_s_km'][0]:.1f}-"
            f"{coverage['uniform_trace_span_s_km'][1]:.1f} km) is not "
            "contained in the plane; increase --extend-se")
    payload = geometry_to_json(geometry)
    OUT_GEOMETRY.write_text(json.dumps(payload, indent=1, sort_keys=True))

    nz, ns = geometry["dike_shape"]
    prov = geometry["wright2006_provenance"]
    print(f"wrote {OUT_GEOMETRY}")
    print(f"geometry sha256: {geometry_sha256(geometry)}")
    print(f"dike plane: {ns} columns x {nz} rows "
          f"({geometry['dike_s_edges_m'][-1] / 1e3:.1f} km x "
          f"{geometry['dike_z_edges_m'][-1] / 1e3:.0f} km), "
          f"{len(prov['segments_nw_to_se'])} segments + "
          f"{args.extend_se} SE extension patches")
    print(f"junction residuals (m): {prov['junction_residuals_m']}")
    print(f"straight-frame (chord) max error: "
          f"{prov['straight_frame_max_error_m']:.0f} m — use patch arrays "
          "for positions")
    print(f"2026 uniform trace spans s = "
          f"{coverage['uniform_trace_span_s_km'][0]:.1f}-"
          f"{coverage['uniform_trace_span_s_km'][1]:.1f} km "
          f"(margins NW {coverage['margin_nw_km']:.1f} km, "
          f"SE {coverage['margin_se_km']:.1f} km; endpoint off-plane "
          f"offsets {coverage['endpoint_offplane_offsets_m']} m)")
    print(f"faults: {prov['fault_treatment']}")


def uperp_panel(result_path: Path) -> None:
    """Render one opening cross-section using Hamling's U_perp panel view."""
    payload = json.loads(Path(result_path).read_text())
    geometry, _ = load_geometry_source(result_path)
    prov = geometry.get("wright2006_provenance")
    if prov is None:
        raise SystemExit(f"{result_path} was not fitted on the Wright plane")
    x = np.asarray(payload["dislocation_magnitudes_m"])
    nz, ns = geometry["dike_shape"]
    opening = x[:nz * ns].reshape(nz, ns)

    s_edges = geometry["dike_s_edges_m"] / 1e3
    z_edges = geometry["dike_z_edges_m"] / 1e3
    s_centres = 0.5 * (s_edges[:-1] + s_edges[1:])
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    # Repeat the boundary patch values at their edges so the contours cover
    # the saved model domain without extending beyond it.
    s_plot = np.concatenate(([s_edges[0]], s_centres, [s_edges[-1]]))
    z_plot = np.concatenate(([z_edges[0]], z_centres, [z_edges[-1]]))
    opening_plot = np.pad(opening, 1, mode="edge")

    vmax = 2.5
    fill_levels = np.linspace(0.0, vmax, 21)
    opening_cmap = LinearSegmentedColormap.from_list(
        "hamling_opening",
        [
            (0.00, "#ffffff"),
            (0.12, "#fff7bc"),
            (0.38, "#fec44f"),
            (0.68, "#f03b20"),
            (1.00, "#b1006e"),
        ],
    )
    extend = "max" if float(opening.max()) > vmax else "neither"

    fig, ax = plt.subplots(figsize=(9.0, 3.0), constrained_layout=True)
    filled = ax.contourf(
        s_plot,
        z_plot,
        opening_plot,
        levels=fill_levels,
        cmap=opening_cmap,
        extend=extend,
        antialiased=True,
    )
    filled.set_rasterized(True)
    line_levels = [level for level in (0.5, 1.0, 1.5, 2.0)
                   if level < float(opening.max())]
    if line_levels:
        ax.contour(
            s_plot,
            z_plot,
            opening_plot,
            levels=line_levels,
            colors="0.15",
            linewidths=0.55,
        )

    n_ext = prov["extension_se"]["n_patches"]
    if n_ext:
        original_se = s_edges[ns - n_ext]
        ax.axvspan(original_se, s_edges[-1], color="0.65", alpha=0.12,
                   zorder=2)
        ax.axvline(original_se, color="0.25", lw=0.9, ls=(0, (4, 3)),
                   zorder=4)
        ax.text(
            original_se + 0.35,
            0.5,
            "2005 SE limit",
            transform=ax.get_xaxis_transform(),
            rotation=90,
            ha="left",
            va="center",
            fontsize=7.5,
            color="0.20",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 1.2},
            zorder=5,
        )
        ax.text(
            0.5 * (original_se + s_edges[-1]),
            0.035,
            f"{s_edges[-1] - original_se:.0f} km SE\nextension",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=7.0,
            color="0.20",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 1.2},
            zorder=5,
        )

    ax.set_xlim(s_edges[0], s_edges[-1])
    ax.set_ylim(z_edges[0], z_edges[-1])
    ax.set_xlabel("Distance along the 2005 dike plane, NW→SE (km)")
    ax.set_ylabel("Modelled depth (km)")
    ax.invert_yaxis()
    ax.set_title("July 2026", loc="left", fontsize=11, fontweight="bold",
                 pad=7)
    ax.set_title(
        f"Wright et al. (2006) plane · 0–{z_edges[-1]:g} km model depth",
        loc="right",
        fontsize=8.5,
        color="0.30",
        pad=8,
    )
    ax.text(0.985, 0.92, "U$_\\perp$", transform=ax.transAxes,
            ha="right", va="top", fontsize=11, fontstyle="italic")
    colorbar = fig.colorbar(
        filled,
        ax=ax,
        label="U$_\\perp$ (m)",
        ticks=np.arange(0.0, vmax + 0.01, 0.5),
        shrink=0.92,
        pad=0.025,
    )
    colorbar.ax.tick_params(labelsize=8)
    chain = "_gslc" if payload.get("nisar_source") == "gslc" else ""
    out = RESULTS / f"uperp_2026_wright2006{chain}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=250)
    plt.close(fig)

    # Edge-independence diagnostic: opening pressed against the SE plane
    # edge, or a large share of volume in the invented extension, means the
    # plane still truncates the source.
    areas = (np.diff(geometry["dike_s_edges_m"])[None, :]
             * np.diff(geometry["dike_z_edges_m"])[:, None])
    volume = float((opening * areas).sum())
    ext_fraction = 0.0
    if n_ext and volume > 0.0:
        ext_fraction = float((opening[:, ns - n_ext:]
                              * areas[:, ns - n_ext:]).sum()) / volume
    print(f"wrote {out}.png and .pdf "
          f"(max opening {opening.max():.2f} m; SE-edge column max "
          f"{opening[:, -1].max():.2f} m; extension holds "
          f"{100 * ext_fraction:.1f}% of the fitted volume)")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dat", type=Path, default=DAT_DEFAULT,
                        help="collaborator geometry table")
    parser.add_argument("--extend-se", type=int, default=5, metavar="N",
                        help="patches appended SE of the southernmost segment")
    parser.add_argument("--no-faults", action="store_true",
                        help="omit the grafted east graben fault")
    parser.add_argument("--panel", type=Path, metavar="RESULT",
                        help="render the U_perp panel from a saved "
                             "distributed fit instead of building the contract")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    if args.panel is not None:
        uperp_panel(args.panel)
    else:
        write_contract(args)
