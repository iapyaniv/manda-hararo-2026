"""joint NISAR + Sentinel-1 co-event decomposition (four geometries).

Solves, per cell of the shared 1.28-km UTM lattice, for two displacement
components from the four co-event LOS geometries (NISAR asc/desc GUNW, LiCSAR
asc/desc):

    d_i = (sin(az) e_i + cos(az) n_i) * h  +  u_i * w

where h is the horizontal displacement along the rift-perpendicular azimuth
az (ENE positive) and w is vertical. The rift-perpendicular constraint uses
the dike symmetry instead of the earlier zero-north assumption: with four
near-polar look vectors the north component is nearly unobservable, and
setting it to zero aliases real rift-normal motion into the other components.
The zero-north variant is reported alongside for comparison with the earlier
products.

Windows end 0-5 days after 18 July, so each dataset is first normalized to
the co-event-through-18-July amplitude with the Sentinel-1 decay summary
(mh26.s1_tracks); the joint solution otherwise mixes different transient
fractions. With two unknowns and up to four observations per cell the
weighted residual is a cross-sensor closure diagnostic, reported as a map.

Run after `mh26.invert` and `mh26.s1_tracks`:  python -m mh26.decompose3d

The free (east, north, up) solve is also rotated into the dike frame —
along strike (positive toward the near-north end of the axis),
perpendicular-left of strike (ENE, the same axis as h), and up — with the
east–north cross-covariance carried into the rotated formal sigmas, so the
rotation is exact and identical to a solve in that basis.

Writes results/decomposition_joint.{json,md} and
results/decomposition_joint.png.
"""

from __future__ import annotations

import datetime as dt
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .gunw import (
    corrected_los,
    load_pair,
    los_unit_neu,
    pair_window,
    reference_farfield,
    valid_sample_mask,
)
from .okada import forward
from .params import (
    CENTER_E,
    CENTER_N,
    FARFIELD_RADIUS,
    INVERSION_PAIRS,
    PROBES,
    RESULTS,
    S1_DECOMPOSITION_SETS,
)
from .s1_tracks import (
    AXIS_PROBE,
    MIN_COEVENT_ABS_R,
    coevent_pixel_stats,
    load_frame,
    load_pair_los,
    reference,
)
from .distributed import (
    NISAR_LABEL,
    S1_GEOMETRY,
    _strike_vectors,
    decay_factor,
    robust_sigma,
)

LATTICE_STEP = 1_280.0
LATTICE_HALF = 40_000.0
MIN_SAMPLES_NISAR = 60      # of 256 pixels per cell (80 m posting)
MIN_SAMPLES_S1 = 40         # of ~130 pixels per cell (LiCSAR grid)
# Two geometries give an exact (h, w) solve — near the axis only the two
# NISAR geometries survive C-band decorrelation. Redundancy (and the closure
# map) exists only where three or more geometries remain.
MIN_GEOMETRIES = 2
PROBE_RADIUS = 2_000.0
FRAME_FRACTION_MIN = 0.9    # NISAR cells must be 90% valid-sample-mask covered
MIN_SAMPLES_NISAR_GSLC = 420    # of ~1820 30-m pixels per cell (same fraction)

S1_COEVENT = dict(S1_DECOMPOSITION_SETS["co-event"])


def _crosscheck_ratio_range() -> str:
    """The measured S1/NISAR amplitude-ratio range, from the saved cross-checks."""
    checks = json.loads(
        (RESULTS / "s1_tracks.json").read_text())["coevent_crosschecks"]
    ratios = sorted(v["amplitude_ratio"] for v in checks.values())
    return f"{ratios[0]:.2f}\u2013{ratios[-1]:.2f}"


def lattice_axes() -> tuple[np.ndarray, np.ndarray]:
    n_cells = int(2 * LATTICE_HALF / LATTICE_STEP)
    centres_e = CENTER_E - LATTICE_HALF + LATTICE_STEP * (np.arange(n_cells) + 0.5)
    centres_n = CENTER_N - LATTICE_HALF + LATTICE_STEP * (np.arange(n_cells) + 0.5)
    return centres_e, centres_n


def bin_cells(X: np.ndarray, Y: np.ndarray, fields: dict[str, np.ndarray],
              valid: np.ndarray, min_samples: int) -> dict[str, np.ndarray]:
    """Cell means of every field over `valid` pixels on the shared lattice."""
    centres_e, centres_n = lattice_axes()
    n_cells = len(centres_e)
    xi = np.floor((X - (CENTER_E - LATTICE_HALF)) / LATTICE_STEP).astype(int)
    yi = np.floor((Y - (CENTER_N - LATTICE_HALF)) / LATTICE_STEP).astype(int)
    ok = valid & (xi >= 0) & (xi < n_cells) & (yi >= 0) & (yi < n_cells)
    flat = yi[ok] * n_cells + xi[ok]
    counts = np.bincount(flat, minlength=n_cells * n_cells)
    out = {}
    for name, field in fields.items():
        sums = np.bincount(flat, weights=field[ok], minlength=n_cells * n_cells)
        means = np.full(n_cells * n_cells, np.nan)
        np.divide(sums, counts, out=means, where=counts >= min_samples)
        out[name] = means.reshape(n_cells, n_cells)
    return out


def load_cells(model: np.ndarray, decay: dict,
               apply_decay: bool = True,
               nisar_source: str = "gunw") -> list[dict]:
    """Each co-event geometry as lattice-cell means of LOS and unit vectors.

    nisar_source "gunw" uses the 80-m product fields; "gslc" uses the
    30-m reprocessed unwrapped fields (mh26.gslc; iono already subtracted),
    with the LOS geometry, valid-sample mask, and window still taken from
    the GUNW crop of the same pair.
    """
    cells = []

    for tag in INVERSION_PAIRS:
        d = load_pair(tag)
        if d is None:
            raise SystemExit(f"{tag}.npz missing from data/nisar")
        if nisar_source == "gslc":
            from .gslc import DATA_GSLC, UNWRAP_COHERENCE_MIN
            f = np.load(DATA_GSLC / f"{tag}_unw.npz")
            los_g = np.where(f["coherence"] > UNWRAP_COHERENCE_MIN,
                             np.asarray(f["los"], dtype=np.float64), np.nan)
            ref = reference_farfield(los_g, f["x"], f["y"])
            Xg, Yg = np.meshgrid(f["x"], f["y"])
            binned = bin_cells(Xg, Yg,
                               {"los": np.where(np.isfinite(ref), ref, 0.0)},
                               np.isfinite(ref), MIN_SAMPLES_NISAR_GSLC)
        else:
            los, _ = corrected_los(d)
            ref = reference_farfield(los, d["x"], d["y"])
            Xg, Yg = np.meshgrid(d["x"], d["y"])
            binned = bin_cells(Xg, Yg,
                               {"los": np.where(np.isfinite(ref), ref, 0.0)},
                               np.isfinite(ref), MIN_SAMPLES_NISAR)
        X, Y = np.meshgrid(d["x"], d["y"])
        # Frame-coverage rule: cells that straddle the swath/subswath
        # boundary or a product-mask gap average partly-filled, edge-degraded
        # samples (such cells produced a spurious rim patch near the
        # ascending frame's western edge). The criterion is the product
        # valid-sample mask alone — value-independent, applied to every
        # NISAR geometry: a cell must be at least 90% mask-valid. Coherence
        # plays no role here, so sparse-but-in-frame axial cells survive.
        # LiCSAR frames have no swath boundary inside the lattice, so no
        # equivalent rule exists for Sentinel-1.
        in_frame = valid_sample_mask(d["mask"]).astype(float)
        frac = bin_cells(X, Y, {"in_frame": in_frame},
                         np.ones(X.shape, dtype=bool), 1)["in_frame"]
        binned["los"][~(frac >= FRAME_FRACTION_MIN)] = np.nan
        # Unit vectors from the product geometry cubes at the cell centres of
        # the valid cells (the cubes are smooth at the cell scale).
        centres_e, centres_n = lattice_axes()
        CE, CN = np.meshgrid(centres_e, centres_n)
        valid_cells = np.isfinite(binned["los"])
        pts = np.column_stack([CN[valid_cells], CE[valid_cells]])
        units = los_unit_neu(d, pts)
        n_grid = np.full(CE.shape, np.nan)
        e_grid = np.full(CE.shape, np.nan)
        u_grid = np.full(CE.shape, np.nan)
        n_grid[valid_cells] = units[:, 0]
        e_grid[valid_cells] = units[:, 1]
        u_grid[valid_cells] = units[:, 2]
        window = pair_window(d)
        end = dt.date(*[int(v) for v in window.split("→")[1].split("-")])
        cells.append(dict(name=NISAR_LABEL[tag], los=binned["los"], n=n_grid,
                          e=e_grid, u=u_grid, window_end=end))

    for frame, pair in S1_COEVENT.items():
        fd = load_frame(frame)
        los = reference(fd, load_pair_los(fd, pair))
        stats = coevent_pixel_stats(fd, los, model)
        if abs(stats["correlation"]) < MIN_COEVENT_ABS_R:
            raise ValueError(f"{frame} {pair}: co-event correlation too weak")
        aligned = stats["applied_observation_sign"] * los
        valid = np.isfinite(aligned)
        fields = {"los": np.where(valid, aligned, 0.0),
                  "n": fd["n"], "e": fd["e"], "u": fd["u"]}
        binned = bin_cells(fd["X"], fd["Y"], fields, valid, MIN_SAMPLES_S1)
        end = dt.date(int(pair[9:13]), int(pair[13:15]), int(pair[15:17]))
        cells.append(dict(name=S1_GEOMETRY[frame], window_end=end, **binned))

    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    R = np.hypot(CE - CENTER_E, CN - CENTER_N)
    for ds in cells:
        factor = decay_factor(ds["window_end"], decay) if apply_decay else 1.0
        ds["decay_factor"] = factor
        ds["los"] = ds["los"] / factor
        far = (R > FARFIELD_RADIUS) & np.isfinite(ds["los"])
        if far.sum() < 50:
            raise ValueError(f"{ds['name']}: too few far-field cells")
        ds["sigma"] = robust_sigma(ds["los"][far])
        print(f"{ds['name']}: {int(np.isfinite(ds['los']).sum())} cells, "
              f"decay factor {factor:.3f}, sigma {100 * ds['sigma']:.1f} cm")
    return cells


def solve_cells(cells: list[dict], perp_az_deg: float | None) -> dict:
    """Per-cell weighted least squares for (h, w) or (E, U) when az is None."""
    shape = cells[0]["los"].shape
    h = np.full(shape, np.nan)
    w = np.full(shape, np.nan)
    closure = np.full(shape, np.nan)
    n_geo = np.zeros(shape, dtype=int)
    sig_h = np.full(shape, np.nan)
    sig_w = np.full(shape, np.nan)

    if perp_az_deg is not None:
        az = np.radians(perp_az_deg)
        horiz = (np.sin(az), np.cos(az))  # (east, north) weights of h
    else:
        horiz = (1.0, 0.0)  # E with N = 0

    stack_los = np.stack([c["los"] for c in cells])
    stack_e = np.stack([c["e"] for c in cells])
    stack_n = np.stack([c["n"] for c in cells])
    stack_u = np.stack([c["u"] for c in cells])
    weights = np.array([1.0 / c["sigma"] for c in cells])

    finite = np.isfinite(stack_los) & np.isfinite(stack_e)
    count = finite.sum(axis=0)
    rows, cols = np.where(count >= MIN_GEOMETRIES)
    for r, c in zip(rows, cols):
        sel = finite[:, r, c]
        a_h = horiz[0] * stack_e[sel, r, c] + horiz[1] * stack_n[sel, r, c]
        A = np.column_stack([a_h, stack_u[sel, r, c]])
        d = stack_los[sel, r, c]
        wt = weights[sel]
        Aw = A * wt[:, None]
        dw = d * wt
        try:
            AtA = Aw.T @ Aw
            cov = np.linalg.inv(AtA)
        except np.linalg.LinAlgError:
            continue
        m = cov @ (Aw.T @ dw)
        resid = dw - Aw @ m
        h[r, c], w[r, c] = m
        n_geo[r, c] = int(sel.sum())
        # Weighted residual in sigma units; only redundant cells close.
        if sel.sum() > 2:
            closure[r, c] = float(np.sqrt(np.mean(resid ** 2)))
        sig_h[r, c], sig_w[r, c] = np.sqrt(np.diag(cov))
    return {"h": h, "w": w, "closure": closure, "n_geometries": n_geo,
            "sigma_h": sig_h, "sigma_w": sig_w}


def solve_cells_3d(cells: list[dict]) -> dict:
    """Free per-cell (east, north, up) solve on cells with >=3 geometries.

    Near-polar LOS geometries observe north weakly, so its formal sigma is
    several times the east/up sigmas; no damping is applied and the sigmas
    are reported so the north field can be weighed accordingly.
    """
    shape = cells[0]["los"].shape
    out = {k: np.full(shape, np.nan) for k in
           ("e", "n", "u", "closure", "sigma_e", "sigma_n", "sigma_u",
            "cov_en")}
    n_geo = np.zeros(shape, dtype=int)
    stack_los = np.stack([c["los"] for c in cells])
    stack_e = np.stack([c["e"] for c in cells])
    stack_n = np.stack([c["n"] for c in cells])
    stack_u = np.stack([c["u"] for c in cells])
    weights = np.array([1.0 / c["sigma"] for c in cells])
    finite = np.isfinite(stack_los) & np.isfinite(stack_e)
    count = finite.sum(axis=0)
    rows, cols = np.where(count >= 3)
    for r, c in zip(rows, cols):
        sel = finite[:, r, c]
        A = np.column_stack([stack_e[sel, r, c], stack_n[sel, r, c],
                             stack_u[sel, r, c]])
        wt = weights[sel]
        Aw = A * wt[:, None]
        dw = stack_los[sel, r, c] * wt
        try:
            cov = np.linalg.inv(Aw.T @ Aw)
        except np.linalg.LinAlgError:
            continue
        m = cov @ (Aw.T @ dw)
        resid = dw - Aw @ m
        out["e"][r, c], out["n"][r, c], out["u"][r, c] = m
        sig = np.sqrt(np.diag(cov))
        out["sigma_e"][r, c], out["sigma_n"][r, c], out["sigma_u"][r, c] = sig
        out["cov_en"][r, c] = cov[0, 1]
        n_geo[r, c] = int(sel.sum())
        if sel.sum() > 3:
            out["closure"][r, c] = float(np.sqrt(np.mean(resid ** 2)))
    out["n_geometries"] = n_geo
    return out


DIKE_FRAME_AXES = {
    "par": "along strike, positive toward strike + 180° (NNW here, the "
           "near-north direction)",
    "perp": "perpendicular-left of strike, positive toward strike − 90° "
            "(ENE here; the axis of the constrained solve's h)",
    "u": "up",
}


def dike_frame(fields: dict, strike_deg: float) -> dict:
    """Rotate a free (east, north, up) solve into the dike frame.

    Uses the `_strike_vectors` frame (along strike, perpendicular-left, up),
    so U⊥ is the same ENE-positive axis as h in `solve_cells`. U∥ is positive
    toward the near-north end of the axis (strike + 180°), so a map draws
    northward motion warm. A linear reparameterization leaves the
    least-squares solution unchanged; the rotated variances need the
    east–north cross-covariance, which the near-polar geometries make large,
    so it is carried explicitly.
    """
    along, perp = _strike_vectors(strike_deg)  # (north, east) unit vectors
    var_e, var_n = fields["sigma_e"] ** 2, fields["sigma_n"] ** 2
    cov = fields["cov_en"]
    out = {"u": fields["u"], "sigma_u": fields["sigma_u"],
           "closure": fields["closure"], "n_geometries": fields["n_geometries"]}
    for name, (v_n, v_e) in (("par", -along), ("perp", perp)):
        out[name] = v_e * fields["e"] + v_n * fields["n"]
        out[f"sigma_{name}"] = np.sqrt(
            v_e ** 2 * var_e + v_n ** 2 * var_n + 2.0 * v_e * v_n * cov)
    return out


def model_dike_frame(model: np.ndarray) -> dict[str, np.ndarray]:
    """Uniform-model (par, perp, up) on the lattice, for comparison."""
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    disp = forward(model, np.column_stack([CN.ravel(), CE.ravel()]))
    along, perp = _strike_vectors(float(model[2]))
    out = {}
    for name, (v_n, v_e) in (("par", -along), ("perp", perp)):
        out[name] = (v_e * disp[:, 1] + v_n * disp[:, 0]).reshape(CE.shape)
    out["u"] = disp[:, 2].reshape(CE.shape)
    return out


def wright_plane_se_strike(reference_strike_deg: float) -> float | None:
    """Strike of the southernmost Wright et al. (2006) dike segment.

    Taken from the saved `mh26.wright2006` contract when present, and
    expressed in the sense within 90° of the reference strike so that the
    along-strike sign convention is the same as for the fitted dike.
    """
    path = RESULTS / "geometry_wright2006.json"
    if not path.exists():
        return None
    provenance = json.loads(path.read_text())["wright2006_provenance"]
    strike = float(provenance["segments_nw_to_se"][-1]["strike_deg"])
    if abs(((strike - reference_strike_deg + 180.0) % 360.0) - 180.0) > 90.0:
        strike = (strike + 180.0) % 360.0
    return strike


def summarize_dike_frame(frame: dict, model_frame: dict, strike_deg: float,
                         fields_3d: dict, alt_strike_deg: float | None) -> dict:
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    near18 = (np.hypot(CE - CENTER_E, CN - CENTER_N) < 18_000) \
        & np.isfinite(frame["par"])
    probes = {comp: {"EAST": probe_mean(frame[comp], PROBES["EAST"]),
                     "WEST": probe_mean(frame[comp], PROBES["WEST"])}
              for comp in ("perp", "par", "u")}
    east, west = probes["perp"]["EAST"], probes["perp"]["WEST"]
    snr = np.abs(frame["par"][near18]) / frame["sigma_par"][near18]
    out = {
        "strike_deg": float(strike_deg),
        "axes": DIKE_FRAME_AXES,
        "flank_probe_cm": probes,
        "model_flank_probe_cm": {
            comp: {"EAST": probe_mean(model_frame[comp], PROBES["EAST"]),
                   "WEST": probe_mean(model_frame[comp], PROBES["WEST"])}
            for comp in ("perp", "par", "u")},
        "opening_divergence_cm": (east - west)
        if east is not None and west is not None else None,
        "median_sigma_cm": {
            comp: float(100 * np.nanmedian(frame[f"sigma_{comp}"]))
            for comp in ("perp", "par", "u")},
        "par_p01_p99_cm_within_18km": [
            float(v) for v in 100 * np.nanpercentile(frame["par"][near18], [1, 99])
        ] if near18.sum() >= 50 else None,
        "par_fraction_above_2sigma_within_18km": float(np.mean(snr > 2.0))
        if near18.sum() >= 50 else None,
    }
    if alt_strike_deg is not None:
        alt = dike_frame(fields_3d, alt_strike_deg)
        out["strike_sensitivity"] = {
            "alternative_strike_deg": float(alt_strike_deg),
            "source": "southernmost Wright et al. (2006) segment "
                      "(results/geometry_wright2006.json)",
            "flank_probe_cm": {
                comp: {"EAST": probe_mean(alt[comp], PROBES["EAST"]),
                       "WEST": probe_mean(alt[comp], PROBES["WEST"])}
                for comp in ("perp", "par")},
        }
    return out


def closure_regions(fields: dict) -> dict:
    """Median closure and >1-sigma fraction, whole scene and near field."""
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    R = np.hypot(CE - CENTER_E, CN - CENTER_N)
    cl = fields["closure"]
    out = {}
    for name, mask in (("all", np.isfinite(cl)),
                       ("within_25km", np.isfinite(cl) & (R < 25_000)),
                       ("within_15km", np.isfinite(cl) & (R < 15_000)),
                       ("beyond_28km", np.isfinite(cl) & (R > FARFIELD_RADIUS))):
        if mask.sum() < 10:
            continue
        vals = cl[mask]
        out[name] = {"n_cells": int(mask.sum()),
                     "median_sigma_units": float(np.median(vals)),
                     "fraction_above_1sigma": float(np.mean(vals > 1.0))}
    return out


def closure_probes(fields: dict) -> dict:
    """Closure inside the probe disks that define the reported values.

    The radial summaries hide local systematics; the probe disks are where
    the reported values live, so their redundant-cell closure is reported
    separately.
    """
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    cl = fields["closure"]
    out = {}
    for name, point in (("EAST", PROBES["EAST"]), ("WEST", PROBES["WEST"]),
                        ("AXIS", AXIS_PROBE)):
        disk = np.hypot(CE - point[0], CN - point[1]) < PROBE_RADIUS
        red = disk & np.isfinite(cl)
        entry = {"n_redundant_cells": int(red.sum())}
        if red.any():
            entry["median_sigma_units"] = float(np.median(cl[red]))
            entry["fraction_above_1sigma"] = float(np.mean(cl[red] > 1.0))
        out[name] = entry
    return out


def probe_geometry_counts(n_geometries: np.ndarray,
                          point: tuple[float, float]) -> dict:
    """How many looks feed the cells inside one probe disk."""
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    mask = (np.hypot(CE - point[0], CN - point[1]) < PROBE_RADIUS) \
        & (n_geometries > 0)
    if not mask.any():
        return {"n_cells": 0}
    vals = n_geometries[mask]
    return {"n_cells": int(mask.sum()),
            "min_geometries": int(vals.min()),
            "max_geometries": int(vals.max())}


def probe_mean(field: np.ndarray, point: tuple[float, float]) -> float | None:
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    mask = (np.hypot(CE - point[0], CN - point[1]) < PROBE_RADIUS) & np.isfinite(field)
    return float(np.nanmean(field[mask]) * 100) if mask.sum() else None


def model_fields(model: np.ndarray, perp_az_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Uniform-model (h, w) on the lattice, for the comparison row."""
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    disp = forward(model, np.column_stack([CN.ravel(), CE.ravel()]))
    az = np.radians(perp_az_deg)
    h = (np.sin(az) * disp[:, 1] + np.cos(az) * disp[:, 0]).reshape(CE.shape)
    w = disp[:, 2].reshape(CE.shape)
    return h, w


def fmt(value: float | None) -> str:
    return "decorrelated" if value is None else f"{value:+.1f} cm"


def summarize(fields: dict, label: str) -> dict:
    div = None
    east = probe_mean(fields["h"], PROBES["EAST"])
    west = probe_mean(fields["h"], PROBES["WEST"])
    if east is not None and west is not None:
        div = east - west
    return {
        "label": label,
        "n_cells": int(np.isfinite(fields["h"]).sum()),
        "n_cells_ge3": int((fields["n_geometries"] >= 3).sum()),
        "n_cells_3geo": int((fields["n_geometries"] == 3).sum()),
        "n_cells_4geo": int((fields["n_geometries"] == 4).sum()),
        "flank_h_probe_cm": {"EAST": east, "WEST": west},
        "flank_probe_geometries": {
            "EAST": probe_geometry_counts(fields["n_geometries"], PROBES["EAST"]),
            "WEST": probe_geometry_counts(fields["n_geometries"], PROBES["WEST"]),
        },
        "opening_divergence_cm": div,
        "axis_w_probe_cm": probe_mean(fields["w"], AXIS_PROBE),
        "axis_probe_geometries": probe_geometry_counts(fields["n_geometries"],
                                                       AXIS_PROBE),
        "median_sigma_h_cm": float(100 * np.nanmedian(fields["sigma_h"])),
        "median_sigma_w_cm": float(100 * np.nanmedian(fields["sigma_w"])),
        "closure": closure_regions(fields),
        "closure_probes": closure_probes(fields),
    }


def plot(fields: dict, model_h: np.ndarray, model_w: np.ndarray,
         perp_az: float, suffix: str = "") -> None:
    centres_e, centres_n = lattice_axes()
    extent = [centres_e[0] / 1e3, centres_e[-1] / 1e3,
              centres_n[0] / 1e3, centres_n[-1] / 1e3]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    v_h = np.nanpercentile(np.abs(fields["h"]), 99)
    v_w = np.nanpercentile(np.abs(fields["w"]), 99)
    panels = [
        (axes[0, 0], fields["h"], v_h, f"observed: rift-perpendicular horizontal\n(az {perp_az:.0f}°, ENE +)"),
        (axes[0, 1], fields["w"], v_w, "observed: vertical (up +)"),
        (axes[0, 2], fields["closure"], None, "cross-sensor closure (weighted RMS, σ units)"),
        (axes[1, 0], model_h, v_h, "uniform model: rift-perpendicular horizontal"),
        (axes[1, 1], model_w, v_w, "uniform model: vertical"),
        (axes[1, 2], fields["h"] - model_h, v_h, "horizontal: observed − model"),
    ]
    for ax, F, v, title in panels:
        if v is None:
            im = ax.imshow(F, cmap="viridis", origin="lower", vmin=0,
                           vmax=np.nanpercentile(F, 98), extent=extent)
            plt.colorbar(im, ax=ax, shrink=0.8, label="σ units")
        else:
            im = ax.imshow(F * 100, cmap="RdBu_r", origin="lower",
                           vmin=-100 * v, vmax=100 * v, extent=extent)
            plt.colorbar(im, ax=ax, shrink=0.8, label="cm")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("UTM 37N easting (km)", fontsize=8)
        ax.set_ylabel("UTM 37N northing (km)", fontsize=8)
    nisar_note = ("30-m reprocessed NISAR fields" if suffix == "_gslc"
                  else "GUNW NISAR fields")
    fig.suptitle("Joint NISAR + Sentinel-1 co-event decomposition on the "
                 f"1.28-km lattice\n(four LOS geometries, {nisar_note}, "
                 "decay-normalized to 18 July; rift-perpendicular horizontal "
                 "+ vertical)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = RESULTS / f"decomposition_joint{suffix}.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"figure: {out}")


def plot_3d(fields: dict, suffix: str = "") -> None:
    centres_e, centres_n = lattice_axes()
    extent = [centres_e[0] / 1e3, centres_e[-1] / 1e3,
              centres_n[0] / 1e3, centres_n[-1] / 1e3]
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    panels = [
        (axes[0, 0], fields["e"], "east (cm)"),
        (axes[0, 1], fields["n"], "north (cm) — weakly observed, see sigmas"),
        (axes[1, 0], fields["u"], "up (cm)"),
        (axes[1, 1], 100 * fields["sigma_n"], "formal sigma north (cm)"),
    ]
    for ax, F, title in panels:
        if title.startswith("formal"):
            im = ax.imshow(F, cmap="viridis", origin="lower", vmin=0,
                           vmax=np.nanpercentile(F, 98), extent=extent)
            plt.colorbar(im, ax=ax, shrink=0.8, label="cm")
        else:
            v = np.nanpercentile(np.abs(F), 99)
            im = ax.imshow(F * 100, cmap="RdBu_r", origin="lower",
                           vmin=-100 * v, vmax=100 * v, extent=extent)
            plt.colorbar(im, ax=ax, shrink=0.8, label="cm")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("UTM 37N easting (km)", fontsize=8)
        ax.set_ylabel("UTM 37N northing (km)", fontsize=8)
    fig.suptitle("Free three-component decomposition (cells with ≥3 looks; "
                 "no direction constraint)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = RESULTS / f"decomposition_joint_3d{suffix}.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"figure: {out}")


def main(nisar_source: str = "gunw") -> None:
    RESULTS.mkdir(exist_ok=True)
    suffix = "" if nisar_source == "gunw" else "_gslc"
    fit_path = RESULTS / "okada_fit.json"
    tracks_path = RESULTS / "s1_tracks.json"
    if not fit_path.exists() or not tracks_path.exists():
        raise SystemExit("run mh26.invert and mh26.s1_tracks first")
    model = np.asarray(json.loads(fit_path.read_text())["values"][:7])
    decay = json.loads(tracks_path.read_text())["post_event_decay"]

    # Rift-perpendicular azimuth from the fitted strike (pointing ENE).
    perp_az = (model[2] - 90.0) % 360.0

    cells = load_cells(model, decay, nisar_source=nisar_source)
    fields = solve_cells(cells, perp_az_deg=perp_az)
    fields_zero_north = solve_cells(cells, perp_az_deg=None)
    model_h, model_w = model_fields(model, perp_az)

    summary = summarize(fields, f"rift-perpendicular (az {perp_az:.1f}°)")
    summary_zn = summarize(fields_zero_north, "zero-north (quasi-East)")

    # Free three-component solve (>=3 looks); no direction constraint.
    fields_3d = solve_cells_3d(cells)
    az = np.radians(perp_az)
    h3d = np.sin(az) * fields_3d["e"] + np.cos(az) * fields_3d["n"]
    east3 = probe_mean(h3d, PROBES["EAST"])
    west3 = probe_mean(h3d, PROBES["WEST"])
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    near18 = (np.hypot(CE - CENTER_E, CN - CENTER_N) < 18_000) \
        & np.isfinite(fields_3d["n"])
    summary_3d = {
        "n_cells": int(np.isfinite(fields_3d["e"]).sum()),
        "n_cells_3geo": int((fields_3d["n_geometries"] == 3).sum()),
        "n_cells_4geo": int((fields_3d["n_geometries"] == 4).sum()),
        "flank_riftperp_probe_cm": {"EAST": east3, "WEST": west3},
        "opening_divergence_cm": (east3 - west3)
        if east3 is not None and west3 is not None else None,
        "axis_u_probe_cm": probe_mean(fields_3d["u"], AXIS_PROBE),
        "axis_probe_geometries": probe_geometry_counts(
            fields_3d["n_geometries"], AXIS_PROBE),
        "north_p01_p99_cm_within_18km": [
            float(v) for v in 100 * np.nanpercentile(fields_3d["n"][near18], [1, 99])
        ] if near18.sum() >= 50 else None,
        "median_sigma_cm": {
            "east": float(100 * np.nanmedian(fields_3d["sigma_e"])),
            "north": float(100 * np.nanmedian(fields_3d["sigma_n"])),
            "up": float(100 * np.nanmedian(fields_3d["sigma_u"])),
        },
        "closure": closure_regions(fields_3d),
        "closure_probes": closure_probes(fields_3d),
    }

    strike = float(model[2])
    frame = dike_frame(fields_3d, strike)
    summary_frame = summarize_dike_frame(
        frame, model_dike_frame(model), strike, fields_3d,
        wright_plane_se_strike(strike))
    for side in ("EAST", "WEST"):
        rotated = summary_frame["flank_probe_cm"]["perp"][side]
        projected = summary_3d["flank_riftperp_probe_cm"][side]
        if rotated is None or projected is None \
                or not np.isclose(rotated, projected, atol=1e-9):
            raise ValueError(f"dike-frame U⊥ {side} probe {rotated} differs "
                             f"from the projected rift-perpendicular probe "
                             f"{projected}")

    # The decay normalization is a model choice; report the unnormalized
    # variant so its effect is visible in the results themselves.
    cells_nd = load_cells(model, decay, apply_decay=False,
                          nisar_source=nisar_source)
    fields_nd = solve_cells(cells_nd, perp_az_deg=perp_az)
    summary_nd = summarize(fields_nd, "no decay normalization")

    payload = {
        "nisar_source": nisar_source,
        "datasets": [{"name": c["name"], "decay_factor": c["decay_factor"],
                      "sigma_cm": 100 * c["sigma"],
                      "n_cells": int(np.isfinite(c["los"]).sum())}
                     for c in cells],
        "no_normalization_variant": {
            "opening_divergence_cm": summary_nd["opening_divergence_cm"],
            "axis_w_probe_cm": summary_nd["axis_w_probe_cm"],
        },
        "perp_azimuth_deg": float(perp_az),
        "rift_perpendicular": summary,
        "zero_north_variant": summary_zn,
        "free_3d": summary_3d,
        "dike_frame": summary_frame,
        "model_probes_cm": {
            "flank_h": {"EAST": probe_mean(model_h, PROBES["EAST"]),
                        "WEST": probe_mean(model_h, PROBES["WEST"])},
            "axis_w": probe_mean(model_w, AXIS_PROBE),
        },
        "notes": [
            "per-cell weighted LS on >=2 of 4 geometries; two-look cells are "
            "exact solves and carry no redundancy",
            "sigma = far-field robust scatter per dataset",
            "decay normalization to co-event through 2026-07-18 (mh26.s1_tracks "
            "summary); a model choice, not a neutral epoch conversion",
            "closure requires >=3 geometries; values in sigma units",
            "north motion is constrained to the rift-perpendicular azimuth "
            "(from the fitted strike), not zero; the free three-component "
            "solve is reported separately",
            "dike_frame is an exact rotation of the free solve into "
            "(along-strike, perpendicular-left, up) with the east–north "
            "cross-covariance carried into the rotated sigmas",
        ],
    }
    (RESULTS / f"decomposition_joint{suffix}.json").write_text(
        json.dumps(payload, indent=2, default=float) + "\n")

    div = summary["opening_divergence_cm"]
    div_model = (payload["model_probes_cm"]["flank_h"]["EAST"] or 0) - \
                (payload["model_probes_cm"]["flank_h"]["WEST"] or 0)
    axis_geo = summary["axis_probe_geometries"]
    axis_label = ""
    if axis_geo.get("n_cells", 0) and axis_geo.get("max_geometries", 0) <= 2:
        axis_label = " — two-look cells only (both NISAR; no Sentinel-1 " \
                     "coverage survives at the axis)"
    close = summary["closure"]
    lines = [
        "# Joint NISAR + Sentinel-1 co-event decomposition"
        + (" (30-m reprocessed NISAR inputs)" if suffix else ""), "",
        f"- four geometries ({', '.join(c['name'] for c in cells)}), "
        f"decay-normalized to 18 July "
        f"(factors {', '.join(f'{c['decay_factor']:.3f}' for c in cells)})"
        + ("; NISAR from the 30-m reprocessed unwrapped fields" if suffix
           else "; NISAR from the GUNW product fields"),
        "",
        "## Constrained decomposition (rift-perpendicular horizontal + vertical)",
        "",
        f"- two components per cell; the horizontal direction is fixed to the "
        f"rift-perpendicular azimuth {perp_az:.1f}° taken from the fitted "
        "strike, so this variant is not independent of the source-model "
        "geometry",
        f"- cells solved: {summary['n_cells']} (two-look cells are exact "
        f"solves; exactly three geometries: {summary['n_cells_3geo']}; "
        f"four geometries: {summary['n_cells_4geo']}; total with at least "
        f"three: {summary['n_cells_ge3']})",
        f"- flank horizontal probes (2-km disks): EAST "
        f"{fmt(summary['flank_h_probe_cm']['EAST'])}, WEST "
        f"{fmt(summary['flank_h_probe_cm']['WEST'])} ⇒ "
        f"**opening divergence {fmt(div)}** "
        f"(uniform model: {div_model:+.1f} cm)",
        f"- axial vertical probe: {fmt(summary['axis_w_probe_cm'])} "
        f"(model {payload['model_probes_cm']['axis_w']:+.1f} cm)"
        + axis_label,
        f"- median formal sigmas: horizontal "
        f"{summary['median_sigma_h_cm']:.1f} cm, vertical "
        f"{summary['median_sigma_w_cm']:.1f} cm",
        "- closure (weighted residual, σ units; redundant cells only): "
        + "; ".join(
            f"{name} median {stats['median_sigma_units']:.2f} "
            f"({100 * stats['fraction_above_1sigma']:.0f}% > 1σ, "
            f"n={stats['n_cells']})"
            for name, stats in close.items()),
        "- probe-disk closure (the disks that define the reported values; "
        "redundant cells only): "
        + "; ".join(
            (f"{name} median {st['median_sigma_units']:.2f} σ "
             f"(n={st['n_redundant_cells']})"
             if "median_sigma_units" in st
             else f"{name} no redundant cells")
            for name, st in summary["closure_probes"].items()),
        "",
        f"- zero-north variant (comparability with earlier products): "
        f"divergence {fmt(summary_zn['opening_divergence_cm'])}, "
        f"axis {fmt(summary_zn['axis_w_probe_cm'])}",
        f"- no-decay-normalization variant: divergence "
        f"{fmt(summary_nd['opening_divergence_cm'])}, axis "
        f"{fmt(summary_nd['axis_w_probe_cm'])} — the normalization is a "
        "model choice, and this is its size",
        "",
        "## Free three-component decomposition (east, north, up; ≥3 looks)",
        "",
        f"- cells solved: {summary_3d['n_cells']} "
        f"(exactly three geometries: {summary_3d['n_cells_3geo']}; "
        f"four geometries: {summary_3d['n_cells_4geo']}); no direction constraint",
        f"- median formal sigmas: east "
        f"{summary_3d['median_sigma_cm']['east']:.1f} cm, north "
        f"{summary_3d['median_sigma_cm']['north']:.1f} cm, up "
        f"{summary_3d['median_sigma_cm']['up']:.1f} cm — north is weakly "
        "observed by near-polar geometries",
        f"- rift-perpendicular flank probes projected from (E, N): EAST "
        f"{fmt(summary_3d['flank_riftperp_probe_cm']['EAST'])}, WEST "
        f"{fmt(summary_3d['flank_riftperp_probe_cm']['WEST'])} ⇒ divergence "
        f"{fmt(summary_3d['opening_divergence_cm'])}",
        f"- axial vertical probe: {fmt(summary_3d['axis_u_probe_cm'])} "
        f"(three-look cells at best near the axis; "
        f"{summary_3d['axis_probe_geometries'].get('n_cells', 0)} probe cells)",
        "- probe-disk closure: "
        + "; ".join(
            (f"{name} median {st['median_sigma_units']:.2f} σ "
             f"(n={st['n_redundant_cells']})"
             if "median_sigma_units" in st
             else f"{name} no redundant cells")
            for name, st in summary_3d["closure_probes"].items()),
        "",
        "## Dike-frame components (rotated free three-component solve)",
        "",
        f"- frame from the fitted strike {strike:.1f}°: U∥ along strike "
        f"(positive toward azimuth {(strike + 180) % 360:.1f}°, NNW — the "
        "near-north direction), U⊥ perpendicular-left (positive toward "
        f"azimuth {(strike - 90) % 360:.1f}°, ENE), up; an exact rotation of "
        "the free solve, so coverage and closure are unchanged",
        f"- flank U⊥ probes: EAST "
        f"{fmt(summary_frame['flank_probe_cm']['perp']['EAST'])}, WEST "
        f"{fmt(summary_frame['flank_probe_cm']['perp']['WEST'])} ⇒ divergence "
        f"{fmt(summary_frame['opening_divergence_cm'])} (identical to the "
        "projected values above)",
        f"- flank U∥ probes: EAST "
        f"{fmt(summary_frame['flank_probe_cm']['par']['EAST'])}, WEST "
        f"{fmt(summary_frame['flank_probe_cm']['par']['WEST'])} (uniform "
        f"model: EAST {fmt(summary_frame['model_flank_probe_cm']['par']['EAST'])}, "
        f"WEST {fmt(summary_frame['model_flank_probe_cm']['par']['WEST'])})",
        f"- median formal sigmas: U⊥ "
        f"{summary_frame['median_sigma_cm']['perp']:.1f} cm, U∥ "
        f"{summary_frame['median_sigma_cm']['par']:.1f} cm, up "
        f"{summary_frame['median_sigma_cm']['u']:.1f} cm — U∥ is nearly "
        "north–south, the direction the near-polar looks observe weakly, so "
        "it inherits most of the north uncertainty",
        (f"- U∥ within 18 km: p01/p99 "
         f"{summary_frame['par_p01_p99_cm_within_18km'][0]:+.1f}/"
         f"{summary_frame['par_p01_p99_cm_within_18km'][1]:+.1f} cm; "
         f"{100 * summary_frame['par_fraction_above_2sigma_within_18km']:.0f}% "
         "of cells exceed 2 σ∥"
         if summary_frame["par_p01_p99_cm_within_18km"] else
         "- U∥ within 18 km: too few cells"),
        (f"- strike sensitivity: with the southernmost 2005 segment strike "
         f"{summary_frame['strike_sensitivity']['alternative_strike_deg']:.1f}° "
         "the flank probes become U⊥ EAST "
         f"{fmt(summary_frame['strike_sensitivity']['flank_probe_cm']['perp']['EAST'])}, "
         f"WEST {fmt(summary_frame['strike_sensitivity']['flank_probe_cm']['perp']['WEST'])}; "
         "U∥ EAST "
         f"{fmt(summary_frame['strike_sensitivity']['flank_probe_cm']['par']['EAST'])}, "
         f"WEST {fmt(summary_frame['strike_sensitivity']['flank_probe_cm']['par']['WEST'])}"
         if "strike_sensitivity" in summary_frame else
         "- strike sensitivity: no saved Wright et al. (2006) contract"),
        "",
        "## Scope",
        "",
        "- The EAST-flank redundant cells close poorly (see the probe-disk",
        "  lines): the looks disagree systematically on the east-lobe",
        "  amplitude, consistent with the Sentinel-1-vs-NISAR amplitude",
        f"  ratios ({_crosscheck_ratio_range()}) in the S1 cross-checks. The joint horizontal",
        "  value there is a cross-sensor compromise, not a closed",
        "  measurement, and the divergence numbers inherit that systematic.",
        "- The rift-perpendicular constraint is exact only if horizontal motion",
        "  is opening-normal; along-strike motion maps into the closure, not",
        "  into h. The free three-component solve has no such constraint but",
        "  carries several-cm north sigmas.",
        "- Near-axis cells lose the C-band geometries first; the axial vertical",
        "  probe of the constrained variant rests on the two NISAR looks alone",
        "  and is not a joint-sensor result.",
        "- Closure is strongly position-dependent: it is smallest in the far",
        "  field and grows in the deforming area, where window differences and",
        "  steep gradients dominate; see the regional closure line above.",
        "- Formal sigmas propagate the far-field dataset scatter only.",
        "- The dike frame depends on the chosen strike; the along-strike",
        "  component is the weakly observed one, and a few degrees of strike",
        "  moves a few percent of U⊥ into U∥ (see the sensitivity line).",
        "- The decay normalization assumes the post-event transient shares the",
        "  co-event dike pattern; it is a model choice (see LIMITATIONS.md).",
    ]
    (RESULTS / f"decomposition_joint{suffix}.md").write_text(
        "\n".join(lines) + "\n")
    print("\n".join(lines))

    plot(fields, model_h, model_w, perp_az, suffix=suffix)
    plot_3d(fields_3d, suffix=suffix)


if __name__ == "__main__":
    import sys

    main(nisar_source="gslc" if "--gslc" in sys.argv else "gunw")
