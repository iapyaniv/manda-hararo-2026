"""distributed opening + graben-fault inversion on four geometries.

Extends the uniform-opening model (mh26.invert) in the two directions the
uniform fit cannot follow:

1. the dike plane is discretized into patches with a free opening
   distribution (positivity + Laplacian smoothing; weight selected by
   held-out prediction);
2. two inward-dipping graben-bounding normal faults flank the dike, with a
   free normal-slip distribution, because the uniform-fit residuals
   concentrate near the axis where graben faulting is unmodeled.

The data are the two NISAR GUNW geometries used by mh26.invert plus the two
LiCSAR co-event pairs (block means with per-block LOS unit vectors). The four
windows end 0–5 days after 18 July, so each dataset is normalized to the
co-event-through-18-July amplitude with the Sentinel-1 decay summary of
mh26.s1_tracks before the joint fit; the uncorrected variant is reported as a
sensitivity.

The dike-plane position, strike and the graben geometry are FIXED a priori
(strike/trace from the uniform fit; fault offset/dip from the residual graben
width); only dislocation magnitudes are inverted, so the linear problem stays
convex. Geometry variants are reported as sensitivities. Comparison with the
2005–2010 dike-plane geometry files will follow once those files are
available.

Run after `mh26.invert` and `mh26.s1_tracks`:  python -m mh26.distributed
Self-checks only:                              python -m mh26.distributed --selftest

Writes results/distributed_fit.{json,md}, results/distributed_planes.png,
and results/distributed_fit_panels.png.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyrocko.modelling import okada_ext
from scipy.optimize import lsq_linear

from .gunw import (
    block_average,
    corrected_los,
    load_pair,
    los_unit_neu,
    pair_window,
    reference_farfield,
)
from .okada import forward
from .params import (
    CENTER_E,
    CENTER_N,
    FARFIELD_RADIUS,
    INVERSION_PAIRS,
    LAME_LAMBDA,
    REPO,
    RESULTS,
    SHEAR_MODULUS,
    S1_DECOMPOSITION_SETS,
)
from .s1_tracks import (
    MIN_COEVENT_ABS_R,
    block_fields,
    coevent_pixel_stats,
    load_frame,
    load_pair_los,
    reference,
)

# --- fixed source geometry -------------------------------------------------
# Dike plane: strike and surface trace from the uniform fit; the patch plane
# extends beyond the uniform rectangle so the data, not the parameterisation,
# bound the opening distribution.
DIKE_HALF_LENGTH = 14_000.0     # m along strike each side of the centroid
DIKE_PATCH_ALONG = 1_750.0      # m
DIKE_DEPTH_TOP = 0.0
DIKE_DEPTH_BOTTOM = 10_000.0
DIKE_PATCH_DEPTH = 1_000.0      # m

# Graben-bounding normal faults: antithetic pair, dipping toward the axis.
# The co-event residual trough and the reprocessed-fringe graben are ~3-4 km
# wide, so the traces sit 2 km either side of the dike trace.
FAULT_OFFSET = 2_000.0          # m horizontal, each side of the dike trace
FAULT_DIP = 65.0                # degrees, toward the axis
FAULT_HALF_LENGTH = 12_250.0    # m along strike (7 patches of 1.75 km per side)
FAULT_PATCH_ALONG = 1_750.0     # m
FAULT_DEPTH_BOTTOM = 4_000.0
FAULT_PATCH_DEPTH = 4_000.0 / 3
MIN_FAULT_PATCH_THICKNESS = 100.0   # m vertical; thinner slivers are dropped

# --- inversion controls ------------------------------------------------------
OPENING_MAX = 15.0              # m, generous physical cap
SLIP_MAX = 10.0                 # m
OFFSET_BOUND = 0.2              # m, per-dataset datum nuisance
# The smoothing weight is chosen by held-out prediction (10-km checkerboard),
# not by L-curve curvature: on this problem the discrete curvature maximum is
# not bracketed by any moderate grid (it keeps climbing into the heavily
# over-damped branch), so a curvature knee is a boundary artifact. The grid
# must bracket the CV minimum; selection fails loudly at a grid edge.
LAMBDA_GRID = np.logspace(-1.5, 3.0, 19)
CV_CELL_M = 10_000.0            # checkerboard cell for the selection criterion
DECAY_REFERENCE = dt.date(2026, 7, 18)

S1_COEVENT = {frame: pair for frame, pair in S1_DECOMPOSITION_SETS["co-event"]}
S1_GEOMETRY = {"014A_07688_131313": "S1 asc 07-12→07-18",
               "079D_07694_131313": "S1 desc 07-11→07-23"}
NISAR_LABEL = {"172_008_A": "NISAR asc 07-08→07-20",
               "165_082_D": "NISAR desc 07-07→07-19"}


def _strike_vectors(strike_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Unit (north, east) vectors along strike and perpendicular-left (ENE)."""
    s = np.radians(strike_deg)
    along = np.array([np.cos(s), np.sin(s)])
    perp = np.array([np.cos(s - np.pi / 2), np.sin(s - np.pi / 2)])
    return along, perp


def _source_path(path: str | Path | None) -> Path:
    return RESULTS / "fault_traces.json" if path is None else Path(path)


def load_fault_trace_source(path: str | Path | None = None) -> tuple[dict, dict]:
    """Per-fault-column trace offsets from a fault-trace JSON file.

    Default source is results/fault_traces.json; an external file with the
    same schema (s_centres_m plus east/west offset_m arrays, dike-frame
    metres, ENE positive) can be supplied instead, e.g. a converted
    2005-2010 fault geometry. Offsets are interpolated at the fault-column
    centres from the trace bins; columns beyond the mapped extent take the
    side's median so the patch grids stay rectangular for the
    regularisation.
    """
    source = _source_path(path)
    raw = source.read_bytes()
    payload = json.loads(raw)
    s_centres = np.asarray(payload["s_centres_m"], dtype=float)
    if s_centres.ndim != 1 or len(s_centres) < 2 \
            or not np.all(np.diff(s_centres) > 0):
        raise ValueError(f"{source}: s_centres_m must be a strictly increasing array")
    f_edges = np.arange(-FAULT_HALF_LENGTH, FAULT_HALF_LENGTH + 1,
                        FAULT_PATCH_ALONG)
    columns = 0.5 * (f_edges[:-1] + f_edges[1:])
    out = {}
    for side in ("east", "west"):
        offsets = np.asarray(payload[side]["offset_m"], dtype=float)
        if offsets.shape != s_centres.shape:
            raise ValueError(f"{source}: {side} offsets do not match s_centres_m")
        ok = np.isfinite(offsets)
        if ok.sum() < 5:
            raise ValueError(f"{source}: too few {side} trace bins")
        filled = np.interp(columns, s_centres[ok], offsets[ok],
                           left=np.nan, right=np.nan)
        filled[~np.isfinite(filled)] = np.median(offsets[ok])
        out[side] = filled
    resolved = source.resolve()
    try:
        display_path = str(resolved.relative_to(REPO.resolve()))
    except ValueError:
        display_path = str(source)
    provenance = {
        "path": display_path,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "trace_support": payload.get("trace_support"),
    }
    return out, provenance


def load_fault_traces(path: str | Path | None = None) -> dict:
    """Load interpolated trace offsets without the provenance record."""
    return load_fault_trace_source(path)[0]


def build_geometry(model: np.ndarray,
                   fault_offset: float = FAULT_OFFSET,
                   fault_dip: float = FAULT_DIP,
                   fault_trace_offsets: dict | None = None,
                   fault_trace_source: str | None = None,
                   enabled_fault_sides: tuple[str, ...] = ("east", "west")) -> dict:
    """Patch tables for the dike plane and the two graben faults.

    Returns pyrocko patch rows [north, east, depth, strike, dip, al1, al2,
    aw1, aw2] (patch centres), unit-dislocation rows, and index/shape
    bookkeeping for the regularisation and the plots.
    """
    e_c, n_c, strike = model[0], model[1], model[2]
    along, perp = _strike_vectors(strike)

    patches, disl, kinds = [], [], []

    # Dike: vertical opening patches.
    s_edges = np.arange(-DIKE_HALF_LENGTH, DIKE_HALF_LENGTH + 1, DIKE_PATCH_ALONG)
    z_edges = np.arange(DIKE_DEPTH_TOP, DIKE_DEPTH_BOTTOM + 1, DIKE_PATCH_DEPTH)
    dike_shape = (len(z_edges) - 1, len(s_edges) - 1)
    dike_s_centres = 0.5 * (s_edges[:-1] + s_edges[1:])
    dike_z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    for z0, z1 in zip(z_edges[:-1], z_edges[1:]):
        for s0, s1 in zip(s_edges[:-1], s_edges[1:]):
            sc = 0.5 * (s0 + s1)
            centre = np.array([n_c, e_c]) + along * sc
            patches.append([centre[0], centre[1], 0.5 * (z0 + z1), strike, 90.0,
                            -(s1 - s0) / 2, (s1 - s0) / 2,
                            -(z1 - z0) / 2, (z1 - z0) / 2])
            disl.append([0.0, 0.0, 1.0])
            kinds.append("dike")

    # Graben faults: normal slip = hanging wall down-dip. pyrocko's positive
    # dip-slip is reverse (verified in selftest), so unit normal slip enters
    # as dip-slip -1 and the inverted parameter stays non-negative.
    #
    # Each inward-dipping column terminates where it meets the vertical dike
    # plane (depth |offset|*tan(dip)): patches must not continue through the
    # dike onto the far side of the rift. The idealized ±2 km geometry is
    # unaffected (its columns meet the dike below the 4 km base), but mapped
    # trace offsets under ~1.9 km truncate their columns. The last patch of a
    # truncated column is shortened so its down-dip edge touches the dike.
    f_edges = np.arange(-FAULT_HALF_LENGTH, FAULT_HALF_LENGTH + 1, FAULT_PATCH_ALONG)
    fz_edges = np.arange(0.0, FAULT_DEPTH_BOTTOM + 1, FAULT_PATCH_DEPTH)
    fault_shape = (len(fz_edges) - 1, len(f_edges) - 1)
    dip_r = np.radians(fault_dip)
    invalid_sides = set(enabled_fault_sides) - {"east", "west"}
    if invalid_sides:
        raise ValueError(f"unknown fault sides: {sorted(invalid_sides)}")
    fault_masks = {}
    fault_z_cut = {}
    column_offsets_by_side = {}
    for side, sign in (("east", +1.0), ("west", -1.0)):
        # pyrocko dips to the right of strike; flip the west fault's strike so
        # both faults dip toward the axis.
        f_strike = strike if side == "east" else strike + 180.0
        if fault_trace_offsets is not None:
            column_offsets = np.asarray(fault_trace_offsets[side], dtype=float)
        else:
            column_offsets = np.full(len(f_edges) - 1, sign * fault_offset)
        if column_offsets.shape != (len(f_edges) - 1,):
            raise ValueError(
                f"{side} trace has {len(column_offsets)} columns; "
                f"expected {len(f_edges) - 1}"
            )
        column_offsets_by_side[side] = column_offsets.copy()
        z_cut = np.minimum(FAULT_DEPTH_BOTTOM,
                           np.abs(column_offsets) * np.tan(dip_r))
        mask = np.zeros(fault_shape, dtype=bool)
        if side not in enabled_fault_sides:
            fault_masks[side] = mask
            fault_z_cut[side] = z_cut
            continue
        for i, (z0, z1) in enumerate(zip(fz_edges[:-1], fz_edges[1:])):
            for j, ((s0, s1), p_off) in enumerate(zip(zip(f_edges[:-1],
                                                          f_edges[1:]),
                                                      column_offsets)):
                z1_eff = min(z1, z_cut[j])
                if z1_eff - z0 < MIN_FAULT_PATCH_THICKNESS:
                    continue
                aw = 0.5 * (z1_eff - z0) / np.sin(dip_r)
                # Centre depth from aw so the top edge is exactly z0 (>= 0);
                # 0.5*(z0+z1_eff) can round the top a hair above the surface,
                # which the Okada kernel rejects.
                zc = z0 + aw * np.sin(dip_r)
                sc = 0.5 * (s0 + s1)
                trace = np.array([n_c, e_c]) + along * sc + perp * p_off
                centre = trace - perp * sign * (zc / np.tan(dip_r))
                patches.append([centre[0], centre[1], zc, f_strike, fault_dip,
                                -(s1 - s0) / 2, (s1 - s0) / 2, -aw, aw])
                disl.append([0.0, -1.0, 0.0])
                kinds.append(f"fault_{side}")
                mask[i, j] = True
        fault_masks[side] = mask
        fault_z_cut[side] = z_cut

    kinds = np.array(kinds)
    return {
        "patches": np.array(patches, dtype=np.float64),
        "disl": np.array(disl, dtype=np.float64),
        "kinds": kinds,
        "areas": _patch_areas(np.array(patches)),
        "dike_shape": dike_shape,
        "fault_shape": fault_shape,
        "dike_s_edges_m": s_edges,
        "dike_z_edges_m": z_edges,
        "dike_s_centres": dike_s_centres,
        "dike_z_centres": dike_z_centres,
        "fault_s_edges_m": f_edges,
        "fault_s_centres_m": 0.5 * (f_edges[:-1] + f_edges[1:]),
        "fault_z_edges_m": fz_edges,
        "min_fault_patch_thickness_m": MIN_FAULT_PATCH_THICKNESS,
        "strike": float(strike),
        "trace_centre": (float(e_c), float(n_c)),
        "fault_offset": float(fault_offset),
        "fault_dip": float(fault_dip),
        "fault_trace_source": (fault_trace_source
                               or ("mapped" if fault_trace_offsets is not None
                                   else "idealized")),
        "enabled_fault_sides": list(enabled_fault_sides),
        "fault_column_offsets_m": {
            side: [float(v) for v in column_offsets_by_side[side]]
            for side in ("east", "west")},
        "fault_masks": fault_masks,
        "fault_termination": "columns end at the dike plane "
                             "(depth = |offset| * tan(dip)) or at "
                             f"{FAULT_DEPTH_BOTTOM / 1e3:.0f} km, whichever "
                             "is shallower",
        "fault_z_cut_km": {side: [float(v) / 1e3 for v in fault_z_cut[side]]
                           for side in ("east", "west")},
        "fault_retained_patches": {side: int(fault_masks[side].sum())
                                   for side in ("east", "west")},
    }


def _json_geometry_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _json_geometry_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_geometry_value(item) for item in value]
    return value


def geometry_to_json(geometry: dict) -> dict:
    """Serialize the complete geometry needed to interpret a saved solution."""
    out = {key: _json_geometry_value(value) for key, value in geometry.items()}
    out["geometry_contract_version"] = 1
    return out


def geometry_from_json(payload: dict) -> dict:
    """Restore and validate a geometry saved by :func:`geometry_to_json`."""
    if payload.get("geometry_contract_version") != 1:
        raise ValueError("distributed result has no complete geometry contract; rerun it")
    geometry = dict(payload)
    for key in ("patches", "disl", "areas", "dike_s_edges_m",
                "dike_z_edges_m", "dike_s_centres", "dike_z_centres",
                "fault_s_edges_m", "fault_s_centres_m", "fault_z_edges_m"):
        geometry[key] = np.asarray(geometry[key], dtype=float)
    geometry["kinds"] = np.asarray(geometry["kinds"], dtype=str)
    for key in ("dike_shape", "fault_shape"):
        geometry[key] = tuple(int(value) for value in geometry[key])
    geometry["trace_centre"] = tuple(float(value)
                                     for value in geometry["trace_centre"])
    geometry["fault_masks"] = {
        side: np.asarray(mask, dtype=bool)
        for side, mask in geometry["fault_masks"].items()
    }
    patches = geometry["patches"]
    disl = geometry["disl"]
    kinds = geometry["kinds"]
    areas = geometry["areas"]
    if patches.ndim != 2 or patches.shape[1] != 9:
        raise ValueError("saved patches must have shape (n, 9)")
    if disl.shape != (len(patches), 3) or kinds.shape != (len(patches),) \
            or areas.shape != (len(patches),):
        raise ValueError("saved geometry arrays have inconsistent patch counts")
    if not np.allclose(areas, _patch_areas(patches), rtol=1e-12, atol=1e-6):
        raise ValueError("saved patch areas do not match the saved patch bounds")
    for side in ("east", "west"):
        mask = geometry["fault_masks"][side]
        if mask.shape != geometry["fault_shape"]:
            raise ValueError(f"saved {side} fault mask has the wrong shape")
        if int(mask.sum()) != int(np.sum(kinds == f"fault_{side}")):
            raise ValueError(f"saved {side} fault mask does not match patch kinds")
    return geometry


def geometry_sha256(geometry: dict) -> str:
    """Stable identity for the complete serialized geometry."""
    encoded = json.dumps(geometry_to_json(geometry), sort_keys=True,
                         separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def _greens_geometry_key(geometry: dict) -> str:
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(geometry["patches"]).tobytes())
    digest.update(np.ascontiguousarray(geometry["disl"]).tobytes())
    return digest.hexdigest()


def _patch_areas(patches: np.ndarray) -> np.ndarray:
    return ((patches[:, 6] - patches[:, 5]) * (patches[:, 8] - patches[:, 7]))


def greens_los(geometry: dict, pts: np.ndarray, units: np.ndarray) -> np.ndarray:
    """LOS Green's matrix (n_pts, n_patches) for unit dislocation per patch."""
    rec = np.zeros((len(pts), 3))
    rec[:, 0] = pts[:, 0]
    rec[:, 1] = pts[:, 1]
    res = okada_ext.okada(geometry["patches"], geometry["disl"], rec,
                          LAME_LAMBDA, SHEAR_MODULUS,
                          nthreads=0, rotate_sdn=0, stack_sources=0)
    # res[i, j, 0:3] = (north, east, down) at receiver j from unit patch i
    disp = np.stack([res[:, :, 0], res[:, :, 1], -res[:, :, 2]], axis=-1)
    return np.einsum("pjk,jk->jp", disp, units)


def laplacian(shape: tuple[int, int]) -> np.ndarray:
    """Second-difference roughness with implicit zero outside the plane edge."""
    nz, ns = shape
    n = nz * ns
    L = np.zeros((n, n))
    for i in range(nz):
        for j in range(ns):
            k = i * ns + j
            L[k, k] = 4.0
            for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < nz and 0 <= jj < ns:
                    L[k, ii * ns + jj] = -1.0
    return L


def decay_factor(window_end: dt.date, decay: dict) -> float:
    """Captured fraction of (co-event + transient) for a window ending at t2.

    All four co-event windows open before the onset bracket, so each captures
    the full co-event step plus the transient through its end date.
    """
    t2 = max(0.0, (window_end - DECAY_REFERENCE).days)
    a, tau = decay["amplitude_total_pattern_fraction"], decay["tau_days"]
    return 1.0 + a * (1.0 - np.exp(-t2 / tau))


def robust_sigma(values: np.ndarray) -> float:
    med = np.median(values)
    return float(1.4826 * np.median(np.abs(values - med)))


def load_datasets(model: np.ndarray, decay: dict, apply_decay: bool = True,
                  nisar_source: str = "gunw") -> list[dict]:
    """The four co-event geometries as weighted block observations.

    nisar_source "gunw" uses the product crops (80 m, product unwrapping);
    "gslc" uses the mh26.gslc reprocessed fields (30 m, own unwrapping),
    with LOS unit vectors still interpolated from the GUNW geometry cubes of
    the same acquisitions.
    """
    datasets = []

    for tag in INVERSION_PAIRS:
        d = load_pair(tag)
        if d is None:
            raise SystemExit(f"{tag}.npz missing from data/nisar — see data/README.md")
        window = pair_window(d)
        if nisar_source == "gslc":
            from .gslc import DATA_GSLC, UNWRAP_COHERENCE_MIN

            g = np.load(DATA_GSLC / f"{tag}_unw.npz")
            los = np.where(g["coherence"] > UNWRAP_COHERENCE_MIN,
                           g["los"].astype(np.float64), np.nan)
            # 1.28-km blocks: the reprocessed grid posting is 30 m.
            pts, vals, grid = block_average(los, g["x"], g["y"], block=42)
            label = NISAR_LABEL[tag].replace("NISAR", "NISAR/GSLC")
        else:
            los, _ = corrected_los(d)
            ref = reference_farfield(los, d["x"], d["y"])
            pts, vals, grid = block_average(ref, d["x"], d["y"])
            label = NISAR_LABEL[tag]
        units = los_unit_neu(d, pts)
        end = dt.date(*[int(v) for v in window.split("→")[1].split("-")])
        datasets.append(dict(name=label, tag=tag, pts=pts, vals=vals,
                             units=units, grid=grid, window_end=end,
                             kind="nisar"))

    for frame, pair in S1_COEVENT.items():
        fd = load_frame(frame)
        los = reference(fd, load_pair_los(fd, pair))
        stats = coevent_pixel_stats(fd, los, model)
        if abs(stats["correlation"]) < MIN_COEVENT_ABS_R:
            raise ValueError(f"{frame} {pair}: co-event correlation too weak "
                             f"({stats['correlation']:.3f}) to fix the frame sign")
        blocks = block_fields(fd, stats["applied_observation_sign"] * los)
        end = dt.date(int(pair[9:13]), int(pair[13:15]), int(pair[15:17]))
        datasets.append(dict(name=S1_GEOMETRY[frame], tag=f"{frame}:{pair}",
                             pts=blocks["pts"], vals=blocks["vals"],
                             units=blocks["units"], grid=blocks["grid"],
                             ok_mask=blocks["ok_mask"], window_end=end,
                             kind="s1"))

    for ds in datasets:
        factor = decay_factor(ds["window_end"], decay) if apply_decay else 1.0
        ds["decay_factor"] = factor
        ds["vals"] = ds["vals"] / factor
        radius = np.hypot(ds["pts"][:, 1] - CENTER_E, ds["pts"][:, 0] - CENTER_N)
        far = radius > FARFIELD_RADIUS
        if far.sum() < 50:
            raise ValueError(f"{ds['name']}: too few far-field blocks for a noise estimate")
        ds["sigma"] = robust_sigma(ds["vals"][far])
        print(f"{ds['name']}: {len(ds['vals'])} blocks, decay factor "
              f"{factor:.3f}, far-field sigma {100 * ds['sigma']:.1f} cm")
    return datasets


def assemble(datasets: list[dict], geometry: dict, with_faults: bool = True):
    """Weighted design matrix [G | dataset offsets] and data vector."""
    keep = np.ones(len(geometry["kinds"]), dtype=bool)
    if not with_faults:
        keep = geometry["kinds"] == "dike"
    G_rows, d_rows, w_rows = [], [], []
    n_sets = len(datasets)
    geometry_key = _greens_geometry_key(geometry)
    for i, ds in enumerate(datasets):
        if "_greens" not in ds or ds.get("_greens_geometry_key") != geometry_key:
            ds["_greens"] = greens_los(geometry, ds["pts"], ds["units"])
            ds["_greens_geometry_key"] = geometry_key
        G = ds["_greens"]
        G = G[:, keep]
        offsets = np.zeros((len(ds["vals"]), n_sets))
        offsets[:, i] = 1.0
        G_rows.append(np.hstack([G, offsets]))
        d_rows.append(ds["vals"])
        w_rows.append(np.full(len(ds["vals"]), 1.0 / ds["sigma"]))
    A = np.vstack(G_rows)
    b = np.concatenate(d_rows)
    w = np.concatenate(w_rows)
    return A * w[:, None], b * w, keep


def regularization_blocks(geometry: dict, keep: np.ndarray) -> list[np.ndarray]:
    """Laplacian blocks in parameter order: dike, then each retained fault.

    Truncated fault columns drop patches, so each fault Laplacian is
    restricted to the retained patches; a dropped neighbour contributes the
    same implicit zero as the plane edge.
    """
    kinds = geometry["kinds"][keep]
    blocks = [laplacian(geometry["dike_shape"])]
    if (kinds != "dike").any():
        full = laplacian(geometry["fault_shape"])
        for side in ("east", "west"):
            m = geometry["fault_masks"][side].ravel()
            blocks.append(full[np.ix_(m, m)])
    return blocks


def solve(A: np.ndarray, b: np.ndarray, keep: np.ndarray, geometry: dict,
          lam: float, n_sets: int) -> np.ndarray:
    kinds = geometry["kinds"][keep]
    blocks = regularization_blocks(geometry, keep)
    n_par = keep.sum()
    L = np.zeros((sum(B.shape[0] for B in blocks), n_par + n_sets))
    r = c = 0
    for B in blocks:
        L[r:r + B.shape[0], c:c + B.shape[1]] = B
        r += B.shape[0]
        c += B.shape[1]
    A_full = np.vstack([A, lam * L])
    b_full = np.concatenate([b, np.zeros(L.shape[0])])
    upper = np.where(kinds == "dike", OPENING_MAX, SLIP_MAX)
    lo = np.concatenate([np.zeros(n_par), -OFFSET_BOUND * np.ones(n_sets)])
    hi = np.concatenate([upper, OFFSET_BOUND * np.ones(n_sets)])
    res = lsq_linear(A_full, b_full, bounds=(lo, hi), method="trf",
                     max_iter=400, tol=1e-12)
    if not res.success:
        raise RuntimeError(f"lsq_linear failed at lambda={lam}: {res.message}")
    return res.x


def weighted_rms(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> float:
    return float(np.sqrt(np.mean((A @ x - b) ** 2)))


def dataset_rms(datasets: list[dict], geometry: dict, keep: np.ndarray,
                x: np.ndarray) -> dict:
    """Physical (unweighted) per-dataset RMS in metres."""
    out = {}
    for i, ds in enumerate(datasets):
        G = ds["_greens"][:, keep]
        pred = G @ x[:keep.sum()] + x[keep.sum() + i]
        out[ds["name"]] = float(np.sqrt(np.mean((pred - ds["vals"]) ** 2)))
        ds[f"_pred_{keep.sum()}"] = pred
    out["all"] = float(np.sqrt(np.mean(np.concatenate(
        [ds[f"_pred_{keep.sum()}"] - ds["vals"] for ds in datasets]) ** 2)))
    return out


def uniform_reference_rms(datasets: list[dict], model: np.ndarray) -> dict:
    """RMS of the fixed uniform-opening model on the same blocks (offsets refit)."""
    out = {}
    residuals = []
    for ds in datasets:
        pred = np.sum(forward(model, ds["pts"]) * ds["units"], axis=1)
        offset = float(np.mean(ds["vals"] - pred))
        r = ds["vals"] - pred - offset
        out[ds["name"]] = float(np.sqrt(np.mean(r ** 2)))
        residuals.append(r)
    out["all"] = float(np.sqrt(np.mean(np.concatenate(residuals) ** 2)))
    return out


def checkerboard_holdout(datasets: list[dict], geometry: dict, lam: float,
                         cell: float = 5_000.0) -> list[dict]:
    """Spatial 5-km checkerboard holdout: fit one half, score the other.

    In-sample RMS always falls when parameters are added; the holdout asks
    whether the fault basis explains structure that repeats across space.
    Both complementary splits are reported for both model variants.
    """
    n_sets = len(datasets)
    masks = []
    for ds in datasets:
        ix = np.floor((ds["pts"][:, 1] - CENTER_E) / cell).astype(int)
        iy = np.floor((ds["pts"][:, 0] - CENTER_N) / cell).astype(int)
        masks.append((ix + iy) % 2)
    out = []
    for parity in (0, 1):
        row = {"train_parity": int(parity)}
        for with_faults in (False, True):
            keep = np.ones(len(geometry["kinds"]), dtype=bool)
            if not with_faults:
                keep = geometry["kinds"] == "dike"
            A_rows, b_rows, test = [], [], []
            for i, ds in enumerate(datasets):
                G = ds["_greens"][:, keep]
                offsets = np.zeros((len(ds["vals"]), n_sets))
                offsets[:, i] = 1.0
                block = np.hstack([G, offsets])
                train = masks[i] == parity
                w = 1.0 / ds["sigma"]
                A_rows.append(block[train] * w)
                b_rows.append(ds["vals"][train] * w)
                test.append((block[~train], ds["vals"][~train]))
            x = solve(np.vstack(A_rows), np.concatenate(b_rows), keep,
                      geometry, lam, n_sets)
            resid = np.concatenate([b - Bt @ x for Bt, b in test])
            key = "with_faults" if with_faults else "opening_only"
            row[f"test_rms_m_{key}"] = float(np.sqrt(np.mean(resid ** 2)))
        out.append(row)
        print(f"holdout parity {parity}: opening-only "
              f"{100 * row['test_rms_m_opening_only']:.2f} cm, with faults "
              f"{100 * row['test_rms_m_with_faults']:.2f} cm")
    return out


def _cv_checkerboard_rms(datasets: list[dict], geometry: dict, lam: float,
                         cell: float) -> list[float]:
    """Held-out test RMS (m) of the with-faults model, both parities."""
    n_sets = len(datasets)
    keep = np.ones(len(geometry["kinds"]), dtype=bool)
    masks = []
    for ds in datasets:
        ix = np.floor((ds["pts"][:, 1] - CENTER_E) / cell).astype(int)
        iy = np.floor((ds["pts"][:, 0] - CENTER_N) / cell).astype(int)
        masks.append((ix + iy) % 2)
    out = []
    for parity in (0, 1):
        A_rows, b_rows, test = [], [], []
        for i, ds in enumerate(datasets):
            G = ds["_greens"][:, keep]
            offsets = np.zeros((len(ds["vals"]), n_sets))
            offsets[:, i] = 1.0
            block = np.hstack([G, offsets])
            train = masks[i] == parity
            w = 1.0 / ds["sigma"]
            A_rows.append(block[train] * w)
            b_rows.append(ds["vals"][train] * w)
            test.append((block[~train], ds["vals"][~train]))
        x = solve(np.vstack(A_rows), np.concatenate(b_rows), keep,
                  geometry, lam, n_sets)
        resid = np.concatenate([b - Bt @ x for Bt, b in test])
        out.append(float(np.sqrt(np.mean(resid ** 2))))
    return out


def sensor_transfer(datasets: list[dict], geometry: dict, lam: float) -> dict:
    """Cross-sensor prediction: fit one sensor, score the other (offsets refit).

    NISAR->S1 holds out the sensor that did not set the source and fault
    geometry. S1->NISAR is conditional because the NISAR test fields affected
    geometry selection before this fit. Both directions also include real
    cross-sensor amplitude differences.
    """
    out = {}
    for train_kind, label in (("nisar", "NISAR->S1"), ("s1", "S1->NISAR")):
        train = [ds for ds in datasets if ds["kind"] == train_kind]
        test = [ds for ds in datasets if ds["kind"] != train_kind]
        row = {}
        for with_faults, key in ((False, "opening_only_cm"),
                                 (True, "with_faults_cm")):
            A_t, b_t, keep_t = assemble(train, geometry,
                                        with_faults=with_faults)
            x_t = solve(A_t, b_t, keep_t, geometry, lam, len(train))
            m = x_t[:keep_t.sum()]
            resid = []
            for ds in test:
                pred = ds["_greens"][:, keep_t] @ m
                off = float(np.mean(ds["vals"] - pred))
                resid.append(ds["vals"] - pred - off)
            row[key] = float(100 * np.sqrt(np.mean(
                np.concatenate(resid) ** 2)))
        row["validation_role"] = (
            "held-out-sensor test; S1 did not set the saved geometry"
            if train_kind == "nisar" else
            "conditional reverse transfer; NISAR test fields affected "
            "geometry selection"
        )
        out[label] = row
    return out


def select_lambda(datasets: list[dict], geometry: dict, A: np.ndarray,
                  b: np.ndarray, keep: np.ndarray, n_sets: int) -> dict:
    """Smoothing weight by held-out prediction over LAMBDA_GRID.

    For each grid value: the full-data solution (kept for the report and the
    chosen fit), the 10-km checkerboard held-out RMS (selection criterion),
    and the discrete L-curve curvature (reported for reference only — see the
    LAMBDA_GRID note). Fails if the CV minimum touches the grid boundary.
    """
    rows, solutions = [], []
    for lam in LAMBDA_GRID:
        x = solve(A, b, keep, geometry, lam, n_sets)
        rough = 0.0
        base = 0
        for B in regularization_blocks(geometry, keep):
            seg = x[base:base + B.shape[1]]
            rough += float(np.sum((B @ seg) ** 2))
            base += B.shape[1]
        per = dataset_rms(datasets, geometry, keep, x)
        cv = _cv_checkerboard_rms(datasets, geometry, lam, CV_CELL_M)
        rows.append({"lambda": float(lam),
                     "weighted_rms": weighted_rms(A, b, x),
                     "rms_all": per["all"],
                     "roughness": rough,
                     "cv10_test_rms_m": cv,
                     "cv10_test_rms_mean_m": float(np.mean(cv))})
        solutions.append(x)
        print(f"lambda {lam:8.3f}: fit rms {100 * per['all']:.2f} cm, "
              f"10-km held-out rms {100 * np.mean(cv):.3f} cm")
    cv_mean = np.array([r["cv10_test_rms_mean_m"] for r in rows])
    chosen = int(np.argmin(cv_mean))
    if chosen in (0, len(LAMBDA_GRID) - 1):
        raise RuntimeError("lambda selection: CV minimum at the grid "
                           "boundary — extend LAMBDA_GRID")
    knee = lcurve_knee(np.array([r["lambda"] for r in rows]),
                       np.array([r["weighted_rms"] for r in rows]),
                       np.array([r["roughness"] for r in rows]))
    lam = float(LAMBDA_GRID[chosen])
    print(f"lambda chosen by 10-km held-out RMS: {lam:.3f} "
          f"(curvature knee for reference: {LAMBDA_GRID[knee]:.3f}"
          f"{', boundary-adjacent' if knee >= len(LAMBDA_GRID) - 2 or knee <= 1 else ''})")
    return {
        "rule": "minimum mean 10-km checkerboard held-out RMS (with-faults "
                "model, both parities); curvature knee reported for "
                "reference only",
        "grid": [float(v) for v in LAMBDA_GRID],
        "sweep": rows,
        "chosen_index": chosen,
        "lambda_chosen": lam,
        "curvature_knee_lambda": float(LAMBDA_GRID[knee]),
        "curvature_knee_boundary_adjacent":
            bool(knee >= len(LAMBDA_GRID) - 2 or knee <= 1),
        "solution": solutions[chosen],
        "sensor_transfer_at_chosen": sensor_transfer(datasets, geometry, lam),
    }


def lcurve_knee(lambdas: np.ndarray, chis: np.ndarray, roughs: np.ndarray) -> int:
    """Index of maximum curvature on the log-log misfit/roughness curve."""
    x = np.log10(np.asarray(chis))
    y = np.log10(np.maximum(np.asarray(roughs), 1e-30))
    curvature = np.full(len(x), -np.inf)
    for i in range(1, len(x) - 1):
        a = np.hypot(x[i] - x[i - 1], y[i] - y[i - 1])
        b = np.hypot(x[i + 1] - x[i], y[i + 1] - y[i])
        c = np.hypot(x[i + 1] - x[i - 1], y[i + 1] - y[i - 1])
        area = abs((x[i] - x[i - 1]) * (y[i + 1] - y[i - 1])
                   - (x[i + 1] - x[i - 1]) * (y[i] - y[i - 1]))
        if a * b * c > 0:
            curvature[i] = 2 * area / (a * b * c)
    return int(np.argmax(curvature))


def moment_summary(geometry: dict, keep: np.ndarray, x: np.ndarray) -> dict:
    kinds = geometry["kinds"][keep]
    areas = geometry["areas"][keep]
    m = x[:keep.sum()]
    dike = kinds == "dike"
    volume = float(np.sum(areas[dike] * m[dike]))
    dike_shape = geometry["dike_shape"]
    opening_grid = m[dike].reshape(dike_shape)
    weights = m[dike].reshape(dike_shape)
    total = weights.sum()
    s_centroid = float((weights.sum(axis=0) * geometry["dike_s_centres"]).sum() / total)
    z_centroid = float((weights.sum(axis=1) * geometry["dike_z_centres"]).sum() / total)
    out = {
        "dike_volume_km3": volume / 1e9,
        "dike_moment_Nm": SHEAR_MODULUS * volume,
        "max_opening_m": float(m[dike].max()),
        "opening_centroid_alongstrike_km": s_centroid / 1e3,
        "opening_centroid_depth_km": z_centroid / 1e3,
        "opening_grid": opening_grid,
    }
    for side in ("east", "west"):
        sel = kinds == f"fault_{side}"
        if sel.any():
            out[f"fault_{side}_max_slip_m"] = float(m[sel].max())
            out[f"fault_{side}_potency_m3"] = float(np.sum(areas[sel] * m[sel]))
            grid = np.full(geometry["fault_shape"], np.nan)
            grid[geometry["fault_masks"][side]] = m[sel]
            out[f"fault_{side}_grid"] = grid
    return out


# --- self-checks -------------------------------------------------------------

def selftest() -> bool:
    ok = True
    model = np.array([683225.36, 1350366.91, 152.58, 17150.0, 1494.0, 4989.0, 2.812])

    # 1. A patch grid tiling the uniform rectangle must reproduce okada.forward.
    e_c, n_c, strike, length, top, bottom, opening = model
    along, _ = _strike_vectors(strike)
    patches, disl = [], []
    for z0, z1 in zip(np.linspace(top, bottom, 5)[:-1], np.linspace(top, bottom, 5)[1:]):
        for s0, s1 in zip(np.linspace(-length / 2, length / 2, 9)[:-1],
                          np.linspace(-length / 2, length / 2, 9)[1:]):
            centre = np.array([n_c, e_c]) + along * 0.5 * (s0 + s1)
            patches.append([centre[0], centre[1], 0.5 * (z0 + z1), strike, 90.0,
                            -(s1 - s0) / 2, (s1 - s0) / 2, -(z1 - z0) / 2, (z1 - z0) / 2])
            disl.append([0.0, 0.0, opening])
    rng = np.random.default_rng(7)
    pts = np.column_stack([n_c + rng.uniform(-3e4, 3e4, 300),
                           e_c + rng.uniform(-3e4, 3e4, 300)])
    rec = np.zeros((len(pts), 3))
    rec[:, 0], rec[:, 1] = pts[:, 0], pts[:, 1]
    res = okada_ext.okada(np.array(patches), np.array(disl), rec,
                          LAME_LAMBDA, SHEAR_MODULUS, nthreads=0,
                          rotate_sdn=0, stack_sources=1)
    tiled = np.column_stack([res[:, 0], res[:, 1], -res[:, 2]])
    direct = forward(model, pts)
    err = np.abs(tiled - direct).max()
    # Superposition is analytically exact; the residual is cancellation noise
    # in the Okada edge terms. 0.1 mm is 1000x below the data noise floor,
    # while a patch-placement error would appear at the 0.1-1 m scale.
    ok_tile = err < 1e-4
    print(f"selftest: tiled-patch vs single-rectangle max |diff| = {err:.2e} "
          f"({'PASS' if ok_tile else 'FAIL'})")
    ok &= ok_tile

    # 2. Unit normal slip on the east graben fault must drop the axis side.
    geometry = build_geometry(model)
    east = geometry["kinds"] == "fault_east"
    x = np.zeros(east.sum())
    x[:] = 1.0
    axis_pt = np.array([[n_c, e_c]])
    _, perp = _strike_vectors(strike)
    outside_pt = np.array([[n_c + 4000 * perp[0], e_c + 4000 * perp[1]]])
    up_axis, up_out = [], []
    for pt, box in ((axis_pt, up_axis), (outside_pt, up_out)):
        rec = np.zeros((1, 3))
        rec[0, :2] = pt[0]
        res = okada_ext.okada(geometry["patches"][east], geometry["disl"][east],
                              rec, LAME_LAMBDA, SHEAR_MODULUS, nthreads=0,
                              rotate_sdn=0, stack_sources=1)
        box.append(-res[0, 2])
    ok_fault = up_axis[0] < 0
    print(f"selftest: normal slip on east fault -> axis up {up_axis[0]:+.3f} m, "
          f"outside up {up_out[0]:+.3f} m ({'PASS (axis subsides)' if ok_fault else 'FAIL'})")
    ok &= ok_fault
    return ok


# --- reporting ---------------------------------------------------------------

def plot_planes(geometry: dict, summary: dict, lambda_label: str,
                suffix: str = "") -> None:
    fig = plt.figure(figsize=(14, 9))
    grid_spec = fig.add_gridspec(2, 2, height_ratios=[1.3, 1.0])

    ax = fig.add_subplot(grid_spec[0, :])
    op = summary["opening_grid"]
    extent = [geometry["dike_s_edges_m"][0] / 1e3,
              geometry["dike_s_edges_m"][-1] / 1e3,
              geometry["dike_z_edges_m"][-1] / 1e3,
              geometry["dike_z_edges_m"][0] / 1e3]
    im = ax.imshow(op, cmap="magma_r", extent=extent, aspect="auto",
                   vmin=0, vmax=max(1e-3, op.max()))
    ax.set_xlabel("along strike from centroid, NNW → SSE (km)")
    ax.set_ylabel("depth (km)")
    ax.set_title(f"dike opening distribution (max {summary['max_opening_m']:.2f} m, "
                 f"volume {summary['dike_volume_km3']:.3f} km³)")
    plt.colorbar(im, ax=ax, label="opening (m)", shrink=0.9)

    fault_grids = [summary[f"fault_{side}_grid"]
                   for side in ("east", "west")
                   if f"fault_{side}_grid" in summary]
    fault_vmax = max(1e-3, *(np.nanmax(grid) for grid in fault_grids))
    for i, side in enumerate(("west", "east")):
        ax = fig.add_subplot(grid_spec[1, i])
        key = f"fault_{side}_grid"
        if key not in summary:
            ax.axis("off")
            continue
        sl = summary[key]
        f_extent = [geometry["fault_s_edges_m"][0] / 1e3,
                    geometry["fault_s_edges_m"][-1] / 1e3,
                    geometry["fault_z_edges_m"][-1] / 1e3,
                    geometry["fault_z_edges_m"][0] / 1e3]
        im = ax.imshow(sl, cmap="viridis_r", extent=f_extent, aspect="auto",
                       vmin=0, vmax=fault_vmax)
        ax.set_title(f"{side} graben fault normal slip "
                     f"(max {summary[f'fault_{side}_max_slip_m']:.2f} m)", fontsize=10)
        ax.set_xlabel("along strike (km)")
        ax.set_ylabel("depth (km)")
        plt.colorbar(im, ax=ax, label="slip (m)", shrink=0.9)

    if geometry["fault_trace_source"] != "idealized":
        enabled = geometry["enabled_fault_sides"]
        medians = "/".join(
            f"{side} {np.median(geometry['fault_column_offsets_m'][side])/1e3:+.1f}"
            for side in enabled
        )
        fault_text = (f"faults at trace columns ({medians} km median, "
                      "truncated at the dike)")
    else:
        fault_text = f"faults ±{geometry['fault_offset']/1e3:.1f} km"
    fig.suptitle("Distributed dislocation model — dike opening + graben faults\n"
                 f"(strike {geometry['strike']:.1f}°, {fault_text} "
                 f"dip {geometry['fault_dip']:.0f}° inward; {lambda_label})",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = RESULTS / f"distributed_planes{suffix}.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"figure: {out}")


def plot_panels(datasets: list[dict], geometry: dict, keep: np.ndarray,
                x: np.ndarray, suffix: str = "") -> None:
    fig, axes = plt.subplots(4, 3, figsize=(15, 18))
    for i, ds in enumerate(datasets):
        pred = ds["_greens"][:, keep] @ x[:keep.sum()] + x[keep.sum() + i]
        if ds["kind"] == "nisar":
            xs, ys, obs_grid = ds["grid"]
            shape_ok = np.isfinite(obs_grid)
            pred_grid = np.full(obs_grid.shape, np.nan)
            pred_grid[shape_ok] = pred
            extent = [xs[0] / 1e3, xs[-1] / 1e3, ys[-1] / 1e3, ys[0] / 1e3]
        else:
            obs_grid = ds["grid"]
            pred_grid = np.full(obs_grid.shape, np.nan)
            pred_grid[ds["ok_mask"]] = pred
            extent = None
        obs_plot = obs_grid.copy()
        v = np.nanpercentile(np.abs(obs_plot), 99)
        for j, (F, title) in enumerate(((obs_plot, f"{ds['name']}\nobserved (decay-normalized)"),
                                        (pred_grid, "distributed model"),
                                        (obs_plot - pred_grid, "residual"))):
            ax = axes[i, j]
            im = ax.imshow(F * 100, cmap="RdBu_r", vmin=-100 * v, vmax=100 * v,
                           extent=extent)
            ax.set_title(title, fontsize=9)
            ax.set_xticks([]), ax.set_yticks([])
            plt.colorbar(im, ax=ax, shrink=0.75, label="LOS (cm)")
    fig.suptitle("Distributed model against the four co-event geometries", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = RESULTS / f"distributed_fit_panels{suffix}.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"figure: {out}")


def main(nisar_source: str = "gunw", fault_source: str = "idealized",
         traces_path: str | Path | None = None) -> None:
    RESULTS.mkdir(exist_ok=True)
    suffix = ("" if nisar_source == "gunw" else "_gslc")
    fit_path = RESULTS / "okada_fit.json"
    tracks_path = RESULTS / "s1_tracks.json"
    if not fit_path.exists() or not tracks_path.exists():
        raise SystemExit("run mh26.invert and mh26.s1_tracks first")
    model = np.asarray(json.loads(fit_path.read_text())["values"][:7])
    decay = json.loads(tracks_path.read_text())["post_event_decay"]

    if not selftest():
        raise SystemExit("selftest failed — do not trust the Green's functions")

    trace_offsets = None
    trace_provenance = None
    enabled_fault_sides = ("east", "west")
    if fault_source != "idealized":
        trace_offsets, trace_provenance = load_fault_trace_source(traces_path)
        if fault_source == "mapped":
            suffix = "_mapped" + suffix
            support = trace_provenance.get("trace_support") or {}
            enabled_fault_sides = tuple(
                side for side in ("east", "west")
                if support.get(side, {}).get("consistent_axis_side_down_step")
            )
            if not enabled_fault_sides:
                raise ValueError("no trace side passes the LOS-step support rule")
        else:
            suffix = f"_traces_{trace_provenance['sha256'][:12]}" + suffix
        trace_provenance["modeled_fault_sides"] = list(enabled_fault_sides)
    geometry = build_geometry(
        model,
        fault_trace_offsets=trace_offsets,
        fault_trace_source=fault_source,
        enabled_fault_sides=enabled_fault_sides,
    )
    datasets = load_datasets(model, decay, nisar_source=nisar_source)
    n_sets = len(datasets)

    A, b, keep = assemble(datasets, geometry, with_faults=True)
    # The smoothing weight is selected once, on the primary configuration
    # (idealized geometry, GUNW inputs); every variant reuses it so the
    # geometry and input comparisons are at matched regularization.
    primary = nisar_source == "gunw" and fault_source == "idealized"
    if primary:
        selection = select_lambda(datasets, geometry, A, b, keep, n_sets)
        lam = selection["lambda_chosen"]
        x = selection["solution"]
    else:
        ref_path = RESULTS / "distributed_fit.json"
        if not ref_path.exists():
            raise SystemExit("run the primary fit (mh26.distributed) first: "
                             "variants reuse its smoothing weight")
        lam = float(json.loads(ref_path.read_text())["lambda_chosen"])
        selection = {"rule": "smoothing weight reused from the primary fit "
                             "(distributed_fit.json) so variants stay "
                             "comparable",
                     "lambda_chosen": lam}
        x = solve(A, b, keep, geometry, lam, n_sets)
        print(f"lambda reused from the primary fit: {lam:.3f}")

    rms_full = dataset_rms(datasets, geometry, keep, x)
    summary = moment_summary(geometry, keep, x)

    # Opening-only variant at the same lambda quantifies what the faults add.
    A2, b2, keep2 = assemble(datasets, geometry, with_faults=False)
    x2 = solve(A2, b2, keep2, geometry, lam, n_sets)
    rms_opening = dataset_rms(datasets, geometry, keep2, x2)
    summary_opening = moment_summary(geometry, keep2, x2)

    rms_uniform = uniform_reference_rms(datasets, model)

    holdout = checkerboard_holdout(datasets, geometry, lam)

    # Sensitivities at the chosen lambda: the smoothing envelope (half a
    # decade each way, so the lambda-dependence of every quoted quantity is
    # visible), sensor subsets, geometry variants, and no decay
    # normalization. Reported as envelopes, not uncertainties.
    sensitivities = []
    slip_sides = "/".join(enabled_fault_sides)

    def slip_values(values: dict) -> str:
        return "/".join(
            f"{values.get(f'fault_{side}_max_slip_m', 0):.2f}"
            for side in enabled_fault_sides
        )

    for env_label, factor in (("smoothing λ×1/3.2", 10 ** -0.5),
                              ("smoothing λ×3.2", 10 ** 0.5)):
        x_e = solve(A, b, keep, geometry, lam * factor, n_sets)
        s_e = moment_summary(geometry, keep, x_e)
        r_e = dataset_rms(datasets, geometry, keep, x_e)
        sensitivities.append({
            "variant": env_label,
            "dike_volume_km3": s_e["dike_volume_km3"],
            "max_opening_m": s_e["max_opening_m"],
            "fault_east_max_slip_m": s_e.get("fault_east_max_slip_m", 0.0),
            "fault_west_max_slip_m": s_e.get("fault_west_max_slip_m", 0.0),
            "fault_potency_sum_m3": s_e.get("fault_east_potency_m3", 0.0)
            + s_e.get("fault_west_potency_m3", 0.0),
            "rms_all_cm": 100 * r_e["all"],
        })
        print(f"sensitivity {env_label}: volume {s_e['dike_volume_km3']:.3f} "
              f"km³, {slip_sides} max slip {slip_values(s_e)} m, "
              f"rms {100 * r_e['all']:.2f} cm")
    for label, kinds_kept in (("NISAR only", ("nisar",)),
                              ("Sentinel-1 only", ("s1",))):
        subset = [ds for ds in datasets if ds["kind"] in kinds_kept]
        A_s, b_s, keep_s = assemble(subset, geometry, with_faults=True)
        x_s = solve(A_s, b_s, keep_s, geometry, lam, len(subset))
        s_s = moment_summary(geometry, keep_s, x_s)
        r_s = dataset_rms(subset, geometry, keep_s, x_s)
        sensitivities.append({
            "variant": label,
            "dike_volume_km3": s_s["dike_volume_km3"],
            "max_opening_m": s_s["max_opening_m"],
            "fault_east_max_slip_m": s_s.get("fault_east_max_slip_m", 0.0),
            "fault_west_max_slip_m": s_s.get("fault_west_max_slip_m", 0.0),
            "fault_potency_sum_m3": s_s.get("fault_east_potency_m3", 0.0)
            + s_s.get("fault_west_potency_m3", 0.0),
            "rms_all_cm": 100 * r_s["all"],
        })
        print(f"sensitivity {label}: volume {s_s['dike_volume_km3']:.3f} km³, "
              f"{slip_sides} max slip {slip_values(s_s)} m")
    if fault_source != "idealized":
        def _shifted(delta: float) -> dict:
            return {"east": np.asarray(trace_offsets["east"]) + delta,
                    "west": np.asarray(trace_offsets["west"]) - delta}
        geometry_variants = [
            ("trace columns 250 m inward",
             {"fault_trace_offsets": _shifted(-250.0),
              "enabled_fault_sides": enabled_fault_sides}),
            ("trace columns 250 m outward",
             {"fault_trace_offsets": _shifted(+250.0),
              "enabled_fault_sides": enabled_fault_sides}),
            ("fault dip 55°",
             {"fault_trace_offsets": trace_offsets, "fault_dip": 55.0,
              "enabled_fault_sides": enabled_fault_sides}),
            ("fault dip 75°",
             {"fault_trace_offsets": trace_offsets, "fault_dip": 75.0,
              "enabled_fault_sides": enabled_fault_sides}),
        ]
    else:
        geometry_variants = [
            ("faults ±1.5 km", {"fault_offset": 1_500.0}),
            ("faults ±2.5 km", {"fault_offset": 2_500.0}),
            ("fault dip 55°", {"fault_dip": 55.0}),
            ("fault dip 75°", {"fault_dip": 75.0}),
        ]
    for label, kwargs in geometry_variants:
        g_v = build_geometry(model, **kwargs)
        for ds in datasets:
            ds.pop("_greens", None)
        A_v, b_v, keep_v = assemble(datasets, g_v, with_faults=True)
        x_v = solve(A_v, b_v, keep_v, g_v, lam, n_sets)
        s_v = moment_summary(g_v, keep_v, x_v)
        r_v = dataset_rms(datasets, g_v, keep_v, x_v)
        sensitivities.append({
            "variant": label,
            "dike_volume_km3": s_v["dike_volume_km3"],
            "max_opening_m": s_v["max_opening_m"],
            "fault_east_max_slip_m": s_v.get("fault_east_max_slip_m", 0.0),
            "fault_west_max_slip_m": s_v.get("fault_west_max_slip_m", 0.0),
            "fault_potency_sum_m3": s_v.get("fault_east_potency_m3", 0.0)
            + s_v.get("fault_west_potency_m3", 0.0),
            "rms_all_cm": 100 * r_v["all"],
        })
        print(f"sensitivity {label}: volume {s_v['dike_volume_km3']:.3f} km³, "
              f"rms {100 * r_v['all']:.2f} cm")

    # A saved dip sweep for trace-column geometry: the truncation depth scales
    # with tan(dip), so the sweep tests dip and fault area together.
    dip_sweep = None
    if fault_source != "idealized":
        dip_sweep = []
        for dip in (55.0, 60.0, 65.0, 70.0, 75.0):
            g_d = build_geometry(model, fault_trace_offsets=trace_offsets,
                                 fault_dip=dip,
                                 enabled_fault_sides=enabled_fault_sides)
            for ds in datasets:
                ds.pop("_greens", None)
            A_d, b_d, keep_d = assemble(datasets, g_d, with_faults=True)
            x_d = solve(A_d, b_d, keep_d, g_d, lam, n_sets)
            s_d = moment_summary(g_d, keep_d, x_d)
            r_d = dataset_rms(datasets, g_d, keep_d, x_d)
            pot_d = (s_d.get("fault_east_potency_m3", 0.0)
                     + s_d.get("fault_west_potency_m3", 0.0))
            dip_sweep.append({
                "dip_deg": dip,
                "rms_all_cm": 100 * r_d["all"],
                "fault_to_dike_potency_ratio":
                    pot_d / (s_d["dike_volume_km3"] * 1e9),
                "fault_east_max_slip_m": s_d.get("fault_east_max_slip_m", 0.0),
                "fault_west_max_slip_m": s_d.get("fault_west_max_slip_m", 0.0),
            })
            print(f"dip sweep {dip:.0f}°: rms {100 * r_d['all']:.4f} cm, "
                  f"potency ratio {dip_sweep[-1]['fault_to_dike_potency_ratio']:.3f}")

    for ds in datasets:
        ds.pop("_greens", None)
    datasets_nd = load_datasets(model, decay, apply_decay=False,
                                nisar_source=nisar_source)
    A_nd, b_nd, keep_nd = assemble(datasets_nd, geometry, with_faults=True)
    x_nd = solve(A_nd, b_nd, keep_nd, geometry, lam, n_sets)
    s_nd = moment_summary(geometry, keep_nd, x_nd)
    sensitivities.append({
        "variant": "no decay normalization",
        "dike_volume_km3": s_nd["dike_volume_km3"],
        "max_opening_m": s_nd["max_opening_m"],
        "fault_east_max_slip_m": s_nd.get("fault_east_max_slip_m", 0.0),
        "fault_west_max_slip_m": s_nd.get("fault_west_max_slip_m", 0.0),
        "fault_potency_sum_m3": s_nd.get("fault_east_potency_m3", 0.0)
        + s_nd.get("fault_west_potency_m3", 0.0),
        "rms_all_cm": 100 * dataset_rms(datasets_nd, geometry,
                                        keep_nd, x_nd)["all"],
    })

    fault_potency = (summary.get("fault_east_potency_m3", 0.0)
                     + summary.get("fault_west_potency_m3", 0.0))
    dike_potency = summary["dike_volume_km3"] * 1e9
    payload = {
        "geometry": geometry_to_json(geometry),
        "geometry_sha256": geometry_sha256(geometry),
        "dislocation_magnitudes_m": [float(value)
                                     for value in x[:keep.sum()]],
        "dataset_offsets_m": {
            ds["name"]: float(x[keep.sum() + index])
            for index, ds in enumerate(datasets)
        },
        "n_parameters": int(keep.sum()) + n_sets,
        "datasets": [{"name": ds["name"], "n_blocks": int(len(ds["vals"])),
                      "decay_factor": ds["decay_factor"],
                      "sigma_cm": 100 * ds["sigma"]} for ds in datasets],
        "lambda_selection": {k: v for k, v in selection.items()
                             if k != "solution"},
        "lambda_chosen": lam,
        "rms_m": {"uniform_reference": rms_uniform,
                  "distributed_opening_only": rms_opening,
                  "distributed_with_faults": rms_full},
        "dike": {k: (v if not isinstance(v, np.ndarray) else None)
                 for k, v in summary.items() if not k.endswith("_grid")},
        "opening_distribution": summary["opening_grid"].tolist(),
        # Dropped (truncated) patches are null, not NaN: the JSON stays RFC
        # 8259 parseable and the mask separates "no patch" from zero slip.
        "fault_slip": {side: [[float(v) if np.isfinite(v) else None
                               for v in row]
                              for row in summary[f"fault_{side}_grid"]]
                       for side in ("east", "west")
                       if f"fault_{side}_grid" in summary},
        "potency_partition": {
            "dike_m3": dike_potency,
            "faults_m3": fault_potency,
            "fault_to_dike_potency_ratio": fault_potency / dike_potency,
            "note": "sum(fault area x slip) / sum(dike area x opening); not "
                    "a moment-tensor ratio",
        },
        "holdout_5km_checkerboard": holdout,
        "sensitivities": sensitivities,
        "nisar_source": nisar_source,
        **({"fault_trace_provenance": trace_provenance}
           if trace_provenance is not None else {}),
        **({"trace_dip_sweep": dip_sweep} if dip_sweep is not None else {}),
        "notes": [
            "geometry fixed a priori; only dislocation magnitudes inverted",
            "datasets normalized to co-event-through-18-July with the S1 decay summary",
            "weights: per-dataset far-field robust sigma; offsets: per-dataset nuisance",
            "regularization: Laplacian per surface, single lambda selected by "
            "minimum mean 10-km checkerboard held-out RMS for the primary fit",
            "comparison with the 2005-2010 dike and fault geometries requires "
            "those geometry files",
        ],
    }
    (RESULTS / f"distributed_fit{suffix}.json").write_text(
        json.dumps(payload, indent=2, default=float, allow_nan=False) + "\n")

    excluded_sides = [side for side in ("east", "west")
                      if trace_provenance is not None
                      and side not in enabled_fault_sides]
    lines = [
        "# Distributed opening + graben-fault model (four geometries"
        + ("; NISAR from GSLC reprocessing" if nisar_source == "gslc" else "")
        + ("; fault-trace geometry" if fault_source != "idealized" else "")
        + ")", "",
        "- fault geometry: "
        + (f"anchored at trace columns from "
           f"{trace_provenance['path']} (SHA-256 "
           f"{trace_provenance['sha256']}); dip 65° inward; each column ends "
           "where it meets the dike plane"
           if fault_source != "idealized"
           else "idealized traces ±2 km from the dike, dip 65° inward"),
        *(["- trace evidence rule: " + ", ".join(excluded_sides)
           + " lacks a consistent LOS step across geometries and is excluded "
             "from this fit; the retained trace geometry remains a "
             "sensitivity, not a fault map"]
          if excluded_sides else []),
        f"- data: {', '.join(ds['name'] for ds in datasets)} "
        f"({sum(len(ds['vals']) for ds in datasets)} blocks)",
        "- decay normalization to 2026-07-18 (factors "
        + ", ".join(f"{ds['decay_factor']:.3f}" for ds in datasets) + ")",
        f"- smoothing: λ = {lam:.2f} "
        + (f"chosen by minimum 10-km checkerboard held-out RMS on the grid "
           f"{LAMBDA_GRID[0]:.2g}–{LAMBDA_GRID[-1]:.0f} (the discrete "
           "L-curve curvature maximum is not bracketed by any moderate grid "
           "on this problem and is reported for reference only)"
           if primary else
           "reused from the primary fit so variants compare at matched "
           "regularization"), "",
        "| model | RMS all | NISAR asc | NISAR desc | S1 asc | S1 desc |",
        "|---|---|---|---|---|---|",
    ]
    for label, rr in (("uniform rectangle (fixed, offsets refit)", rms_uniform),
                      ("distributed opening", rms_opening),
                      ("distributed opening + graben faults", rms_full)):
        cells = [rr["all"]] + [rr[ds["name"]] for ds in datasets]
        lines.append("| " + label + " | "
                     + " | ".join(f"{100 * c:.2f} cm" for c in cells) + " |")
    scope_fault_lines = ([
        "- The retained east-fault term improves both sensor-subset fits, but",
        "  its magnitude changes between them; treat its slip values as",
        "  regularized magnitudes, not resolved slip. The west trace is absent",
        "  because it fails the cross-geometry LOS-step support rule.",
    ] if enabled_fault_sides == ("east",) else [
        "- The total fault-like contribution is repeatable across data subsets,",
        "  but the east/west partition is not (see the sensor-subset rows);",
        "  treat per-fault slip values as regularized magnitudes, not resolved",
        "  slip. The west-fault maximum sits at the model-grid edge.",
    ])
    fault_maxima = ", ".join(
        f"{side} {summary.get(f'fault_{side}_max_slip_m', 0):.2f} m"
        for side in enabled_fault_sides
    )
    sensitivity_slips = "–".join(enabled_fault_sides) + " max slip m"
    dip_slips = "/".join(enabled_fault_sides) + " max slip"
    lines.extend([
        "",
        "- the fault term's in-sample test is the distributed-opening row "
        "against the with-faults row; the uniform row is a fixed reference "
        "with many fewer parameters, not a fault test",
        "",
        f"- dike volume **{summary['dike_volume_km3']:.3f} km³** "
        f"(opening-only variant {summary_opening['dike_volume_km3']:.3f}); "
        f"max opening **{summary['max_opening_m']:.2f} m** at "
        f"{summary['opening_centroid_depth_km']:.1f} km centroid depth",
        f"- {'graben faults' if len(enabled_fault_sides) > 1 else 'retained fault term'}: "
        f"regularized maximum normal slip {fault_maxima}; "
        f"fault/dike potency ratio "
        f"{payload['potency_partition']['fault_to_dike_potency_ratio']:.3f} "
        "(Σ area·slip / Σ area·opening; not a moment-tensor ratio)",
        "",
        "## Spatial holdout (5-km checkerboard)",
        "",
    ])
    for row in holdout:
        lines.append(
            f"- train parity {row['train_parity']}: held-out RMS "
            f"opening-only {100 * row['test_rms_m_opening_only']:.2f} cm → "
            f"with faults {100 * row['test_rms_m_with_faults']:.2f} cm")
    lines.extend([
        "",
        f"## Sensitivities (volume km³ / max opening m / {sensitivity_slips} / RMS cm)",
        "",
    ])
    for s in sensitivities:
        slip_values = "–".join(
            f"{s[f'fault_{side}_max_slip_m']:.2f}"
            for side in enabled_fault_sides
        )
        lines.append(f"- {s['variant']}: {s['dike_volume_km3']:.3f} / "
                     f"{s['max_opening_m']:.2f} / {slip_values} / "
                     f"{s['rms_all_cm']:.2f}")
    if dip_sweep is not None:
        lines.extend([
            "",
            "## Dip sweep (trace-column anchoring; truncation depth scales "
            "with tan dip)",
            "",
            f"| dip | RMS all | potency ratio | {dip_slips} |",
            "|---|---|---|---|",
        ])
        for row in dip_sweep:
            slip_values = " / ".join(
                f"{row[f'fault_{side}_max_slip_m']:.2f}"
                for side in enabled_fault_sides
            )
            lines.append(f"| {row['dip_deg']:.0f}° | "
                         f"{row['rms_all_cm']:.4f} cm | "
                         f"{row['fault_to_dike_potency_ratio']:.3f} | "
                         f"{slip_values} m |")
    if primary:
        lines.extend([
            "",
            "## Smoothing selection (10-km checkerboard held-out RMS)",
            "",
            "| λ | fit RMS | held-out RMS (parity mean) |",
            "|---|---|---|",
        ])
        for r in selection["sweep"]:
            marker = " ←" if r["lambda"] == lam else ""
            lines.append(f"| {r['lambda']:.3g} | {100 * r['rms_all']:.2f} cm "
                         f"| {100 * r['cv10_test_rms_mean_m']:.3f} cm"
                         f"{marker} |")
        st = selection["sensor_transfer_at_chosen"]
        lines.extend([
            "",
            "- curvature knee for reference: λ = "
            f"{selection['curvature_knee_lambda']:.3g}"
            + (" (boundary-adjacent — not usable as a selection)"
               if selection["curvature_knee_boundary_adjacent"] else ""),
            "- cross-sensor transfer at the chosen λ (offsets refit): "
            + "; ".join(
                f"{k}: opening-only {v['opening_only_cm']:.2f} → with faults "
                f"{v['with_faults_cm']:.2f} cm" for k, v in st.items()),
            "- validation role: NISAR→S1 is the held-out-sensor direction. "
            "S1→NISAR is conditional because the NISAR fields affected "
            "geometry selection.",
        ])
    lines.extend([
        "",
        "## Scope",
        "",
        "- The dike plane, strike, and fault geometry are fixed a priori; the",
        "  inversion solves only dislocation magnitudes (convex, positive).",
    ] + ([
        "- The fault surfaces keep one fixed strike while their column",
        "  offsets follow the supplied traces, so plan-view steps between",
        "  columns are parameterisation steps, not resolved geometry; the",
        "  smoothing couples stepped columns as neighbours. Columns end at",
        "  the dike plane. A direct test against the 2005-2010 dike and",
        "  fault geometries requires those geometry files.",
    ] if fault_source != "idealized" else [
        "- The graben-fault geometry is a two-plane idealization chosen from the",
        "  residual trough; a direct test against the 2005-2010 dike and fault",
        "  geometries requires those geometry files.",
    ]) + scope_fault_lines + [
        "- Formal resolution/uncertainty of the slip distribution is not",
        "  computed here; the sensitivity list brackets the main choices.",
    ])
    md = "\n".join(lines) + "\n"
    (RESULTS / f"distributed_fit{suffix}.md").write_text(md)
    print(md)

    plot_planes(geometry, summary,
                f"λ = {lam:.2f} by 10-km held-out RMS" if primary
                else f"λ = {lam:.2f}, reused from the primary fit",
                suffix=suffix)
    # Rebuild the Green's cache for the chosen geometry before plotting panels.
    for ds in datasets:
        ds.pop("_greens", None)
    A, b, keep = assemble(datasets, geometry, with_faults=True)
    plot_panels(datasets, geometry, keep, x, suffix=suffix)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gslc", action="store_true",
                        help="use the reprocessed NISAR fields")
    parser.add_argument("--selftest", action="store_true",
                        help="run Green-function checks and exit")
    traces = parser.add_mutually_exclusive_group()
    traces.add_argument("--mapped", action="store_true",
                        help="use results/fault_traces.json")
    traces.add_argument("--traces", type=Path, metavar="FILE",
                        help="use an external fault-trace JSON file; outputs "
                             "receive a content-hash suffix")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    if args.selftest:
        raise SystemExit(0 if selftest() else 1)
    source = "external" if args.traces is not None else (
        "mapped" if args.mapped else "idealized"
    )
    main(nisar_source="gslc" if args.gslc else "gunw",
         fault_source=source, traces_path=args.traces)
