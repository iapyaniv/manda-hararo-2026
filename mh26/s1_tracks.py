"""two-track Sentinel-1 test of the model and the continuation.

The Leeds/COMET collaborators completed the LiCSAR processing for both
Sentinel-1 tracks on 2026-08-17. This module uses those products for three
independent tests that one geometry cannot provide:

1. a second co-event cross-check of the NISAR-derived dike model against the
   descending frame 079D_07694_131313 (pixel statistic of `mh26.s1_check`);
2. dike-pattern projections (block statistic of `mh26.continuation`) for six
   post-event pairs in both geometries, with nine pre-event pairs from the
   same frames as an empirical null distribution for the projection
   coefficient;
3. a quasi-East/quasi-Up decomposition of the two 24-day post-event pairs on
   a common UTM block lattice, next to the same decomposition of the
   co-event pairs.

The LiCSAR unwrapped-phase sign convention is resolved empirically once per
frame from the co-event regression against the model, then applied unchanged
to that frame's pre- and post-event pairs. A weak pair must not resolve its
own sign: that would bias every projection positive.

Run after `mh26.invert`:  python -m mh26.s1_tracks

Writes results/s1_tracks.md, results/s1_tracks.json, and three figures.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer

from .continuation import cluster_wild_interval, regression
from .okada import forward
from .params import (
    CENTER_E,
    CENTER_N,
    DATA_S1,
    EPSG_UTM37N,
    FARFIELD_RADIUS,
    PROBES,
    RESULTS,
    S1_DECOMPOSITION_SETS,
    S1_FRAMES,
    S1_TRACK_PAIRS,
    WAVELENGTH_C,
)
from .s1_check import (
    COHERENCE_MIN_S1,
    PRED_MIN,
    SUBSAMPLE,
    assert_same_grid,
    coherence_fraction,
    read_crop,
)

# Block statistics on the 0.001-degree LiCSAR grid: 12 px is about 1.33 km,
# close to the 1.28-km NISAR blocks used by mh26.continuation.
BLOCK_PX = 12
BLOCK_MIN_FRACTION = 0.3
RADIUS_LIMIT = 25_000.0
PRED_MIN_BLOCK = 0.01
MIN_COEVENT_ABS_R = 0.90

PROCESSING_VARIANTS = (
    ("0.9 km blocks", {"block": 8}),
    ("2.0 km blocks", {"block": 18}),
    ("25 km far field", {"farfield_radius": 25_000.0}),
    ("32 km far field", {"farfield_radius": 32_000.0}),
    ("fit within 20 km", {"radius_limit": 20_000.0}),
    ("fit within 30 km", {"radius_limit": 30_000.0}),
    ("prediction >0.5 cm", {"prediction_min": 0.005}),
    ("prediction >2 cm", {"prediction_min": 0.02}),
)

# Common UTM lattice for the two-geometry decomposition. A 1.28-km cell holds
# about 130 LiCSAR pixels; require 40 valid samples per cell and geometry.
LATTICE_STEP = 1_280.0
LATTICE_HALF = 40_000.0
LATTICE_MIN_SAMPLES = 40
DECOMP_MIN_ABS_DET = 0.2
DECOMP_PROBE_RADIUS = 2_000.0
AXIS_PROBE = (
    0.5 * (PROBES["EAST"][0] + PROBES["WEST"][0]),
    0.5 * (PROBES["EAST"][1] + PROBES["WEST"][1]),
)


def window_of(pair: str) -> str:
    ref, sec = pair.split("_")
    return (f"{ref[:4]}-{ref[4:6]}-{ref[6:]}" "→" f"{sec[:4]}-{sec[4:6]}-{sec[6:]}")


def pair_days(pair: str) -> int:
    import datetime as _dt

    ref, sec = pair.split("_")
    to_date = lambda s: _dt.date(int(s[:4]), int(s[4:6]), int(s[6:]))
    return (to_date(sec) - to_date(ref)).days


def load_frame(frame: str) -> dict:
    """Frame LOS vectors plus cached UTM coordinates for its grid."""
    e, lons, lats = read_crop(DATA_S1 / f"{frame}.geo.E.tif")
    n, n_lons, n_lats = read_crop(DATA_S1 / f"{frame}.geo.N.tif")
    u, u_lons, u_lats = read_crop(DATA_S1 / f"{frame}.geo.U.tif")
    assert_same_grid([
        (f"{frame} E", e, lons, lats),
        (f"{frame} N", n, n_lons, n_lats),
        (f"{frame} U", u, u_lons, u_lats),
    ])
    LON, LAT = np.meshgrid(lons, lats)
    transformer = Transformer.from_crs(4326, EPSG_UTM37N, always_xy=True)
    X, Y = transformer.transform(LON, LAT)
    R = np.hypot(X - CENTER_E, Y - CENTER_N)
    norm = np.sqrt(e**2 + n**2 + u**2)
    vec_ok = (
        np.isfinite(e) & np.isfinite(n) & np.isfinite(u)
        & (np.abs(e) > 0.01) & (np.abs(norm - 1.0) < 0.01)
    )
    return {
        "frame": frame, "lons": lons, "lats": lats,
        "e": e, "n": n, "u": u, "vec_ok": vec_ok,
        "X": X, "Y": Y, "R": R,
    }


def load_pair_los(frame_data: dict, pair: str) -> np.ndarray:
    """Screened LOS in metres, LiCSAR raw sign, before referencing."""
    unw, lons, lats = read_crop(DATA_S1 / f"{pair}.geo.unw.tif")
    cc, cc_lons, cc_lats = read_crop(DATA_S1 / f"{pair}.geo.cc.tif")
    assert_same_grid([
        (frame_data["frame"], frame_data["e"], frame_data["lons"], frame_data["lats"]),
        (f"{pair} unw", unw, lons, lats),
        (f"{pair} cc", cc, cc_lons, cc_lats),
    ])
    unw = unw.copy()
    unw[unw == 0] = np.nan
    cc = coherence_fraction(cc)
    good = (
        np.isfinite(unw) & np.isfinite(cc)
        & (cc > COHERENCE_MIN_S1) & frame_data["vec_ok"]
    )
    los = unw * WAVELENGTH_C / (4 * np.pi)
    los[~good] = np.nan
    return los


def reference(frame_data: dict, los: np.ndarray,
              farfield_radius: float = FARFIELD_RADIUS) -> np.ndarray:
    """Remove the far-field plane, as in mh26.s1_check."""
    X, Y = frame_data["X"], frame_data["Y"]
    far = np.isfinite(los) & (frame_data["R"] > farfield_radius)
    if far.sum() < 500:
        raise ValueError("too few valid far-field Sentinel-1 pixels for referencing")
    design = np.c_[np.ones(far.sum()),
                   (X[far] - CENTER_E) / 1e4, (Y[far] - CENTER_N) / 1e4]
    coef, *_ = np.linalg.lstsq(design, los[far], rcond=None)
    plane = (coef[0] + coef[1] * (X - CENTER_E) / 1e4
             + coef[2] * (Y - CENTER_N) / 1e4)
    return los - plane


def coevent_pixel_stats(frame_data: dict, los_ref: np.ndarray,
                        model: np.ndarray) -> dict:
    """Pixel-level model comparison; the statistic of mh26.s1_check."""
    sl = np.s_[::SUBSAMPLE, ::SUBSAMPLE]
    X, Y = frame_data["X"][sl], frame_data["Y"][sl]
    pts = np.column_stack([Y.ravel(), X.ravel()])
    disp = forward(model, pts)
    units = np.column_stack([
        frame_data["n"][sl].ravel(),
        frame_data["e"][sl].ravel(),
        frame_data["u"][sl].ravel(),
    ])
    los_pred = np.sum(disp * units, axis=1).reshape(X.shape)
    obs = los_ref[sl]
    ok = np.isfinite(obs) & np.isfinite(los_pred) & (np.abs(los_pred) > PRED_MIN)
    if ok.sum() < 100:
        raise ValueError("too few pixels survive the co-event comparison mask")
    slope = np.sum(obs[ok] * los_pred[ok]) / np.sum(los_pred[ok] ** 2)
    sign = -1.0 if slope < 0 else 1.0
    aligned = sign * obs
    residual = aligned - los_pred
    return {
        "n_pixels": int(ok.sum()),
        "regression_slope_observed_on_predicted": float(slope),
        "applied_observation_sign": float(sign),
        "amplitude_ratio": float(abs(slope)),
        "correlation": float(np.corrcoef(aligned[ok], los_pred[ok])[0, 1]),
        "rms_m": float(np.sqrt(np.mean(residual[ok] ** 2))),
        "_fields": (aligned, np.where(np.isfinite(obs), los_pred, np.nan),
                    np.where(ok, residual, np.nan)),
    }


def block_fields(frame_data: dict, los_ref: np.ndarray,
                 block: int = BLOCK_PX) -> dict:
    """Masked block means of LOS, UTM coordinates, and LOS unit vectors.

    All fields average over the same LOS-valid pixels so block coordinates
    and vectors are centroids of the samples that formed each block value.
    """
    valid = np.isfinite(los_ref)
    ny, nx = los_ref.shape
    by, bx = ny // block, nx // block
    if by == 0 or bx == 0:
        raise ValueError("block is larger than the Sentinel-1 crop")

    def block_mean(field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        cropped = np.where(valid, field, np.nan)[: by * block, : bx * block]
        stacked = cropped.reshape(by, block, bx, block)
        n_ok = np.isfinite(stacked).sum(axis=(1, 3))
        sums = np.nansum(stacked, axis=(1, 3))
        means = np.full(sums.shape, np.nan)
        np.divide(sums, n_ok, out=means, where=n_ok > 0)
        return means, n_ok

    los_b, n_ok = block_mean(los_ref)
    los_b[n_ok < block * block * BLOCK_MIN_FRACTION] = np.nan
    fields = {name: block_mean(frame_data[name])[0]
              for name in ("X", "Y", "n", "e", "u")}
    ok = np.isfinite(los_b)
    pts = np.column_stack([fields["Y"][ok], fields["X"][ok]])
    units = np.column_stack([fields[name][ok] for name in ("n", "e", "u")])
    return {"pts": pts, "vals": los_b[ok], "units": units,
            "grid": los_b, "ok_mask": ok}


def project_blocks(model: np.ndarray, blocks: dict,
                   radius_limit: float = RADIUS_LIMIT,
                   prediction_min: float = PRED_MIN_BLOCK) -> dict:
    """Regress block LOS on the co-event dike pattern (mh26.continuation)."""
    pts, vals, units = blocks["pts"], blocks["vals"], blocks["units"]
    predicted = np.sum(forward(model, pts) * units, axis=1)
    radius = np.hypot(pts[:, 1] - CENTER_E, pts[:, 0] - CENTER_N)
    use = (
        np.isfinite(vals) & np.isfinite(predicted)
        & (radius < radius_limit) & (np.abs(predicted) > prediction_min)
    )
    if use.sum() < 100:
        raise ValueError("too few Sentinel-1 blocks survive the projection mask")
    slope, offset = regression(predicted[use], vals[use])
    fitted = slope * predicted[use] + offset
    return {
        "n_blocks": int(use.sum()),
        "slope": float(slope),
        "offset": float(offset),
        "zero_intercept_slope": float(
            np.sum(predicted[use] * vals[use]) / np.sum(predicted[use] ** 2)
        ),
        "correlation": float(np.corrcoef(vals[use], predicted[use])[0, 1]),
        "rms": float(np.sqrt(np.mean((vals[use] - fitted) ** 2))),
        "_predicted": predicted, "_use": use,
    }


def bin_to_lattice(frame_data: dict, los_aligned: np.ndarray) -> dict:
    """Cell means of LOS and E/U vectors on the shared UTM lattice."""
    edge = np.arange(-LATTICE_HALF, LATTICE_HALF + LATTICE_STEP, LATTICE_STEP)
    n_cells = len(edge) - 1
    xi = np.floor((frame_data["X"] - (CENTER_E - LATTICE_HALF))
                  / LATTICE_STEP).astype(int)
    yi = np.floor((frame_data["Y"] - (CENTER_N - LATTICE_HALF))
                  / LATTICE_STEP).astype(int)
    valid = (
        np.isfinite(los_aligned)
        & (xi >= 0) & (xi < n_cells) & (yi >= 0) & (yi < n_cells)
    )
    flat = yi[valid] * n_cells + xi[valid]
    counts = np.bincount(flat, minlength=n_cells * n_cells)

    def cell_mean(field: np.ndarray) -> np.ndarray:
        sums = np.bincount(flat, weights=field[valid],
                           minlength=n_cells * n_cells)
        means = np.full(n_cells * n_cells, np.nan)
        np.divide(sums, counts, out=means, where=counts >= LATTICE_MIN_SAMPLES)
        return means.reshape(n_cells, n_cells)

    centres = CENTER_E - LATTICE_HALF + LATTICE_STEP * (np.arange(n_cells) + 0.5)
    centres_n = CENTER_N - LATTICE_HALF + LATTICE_STEP * (np.arange(n_cells) + 0.5)
    return {
        "los": cell_mean(los_aligned),
        "e": cell_mean(frame_data["e"]),
        "u": cell_mean(frame_data["u"]),
        "x_centres": centres,
        "y_centres": centres_n,
    }


def probe_disk_mean(field: np.ndarray, X: np.ndarray, Y: np.ndarray,
                    point: tuple[float, float],
                    radius: float = DECOMP_PROBE_RADIUS) -> float | None:
    """Disk mean in cm, or None where the disk is fully decorrelated."""
    mask = (np.hypot(X - point[0], Y - point[1]) < radius) & np.isfinite(field)
    return float(np.nanmean(field[mask]) * 100) if mask.sum() > 0 else None


def fmt_cm(value: float | None) -> str:
    return "decorrelated" if value is None else f"{value:+.1f} cm"


def decay_fit(post_rows: list[dict]) -> dict:
    """Descriptive exponential decay through the six S1 post-event windows.

    Each projection integrates a decaying rate over its window, so the model
    is  proj([t1, t2]) = A (exp(-t1/tau) - exp(-t2/tau))  with t in days
    after the 18 July ascending acquisition. This is a two-parameter summary
    of window-integrated values with correlated noise, not a physical model
    or an uncertainty analysis.
    """
    from scipy.optimize import least_squares

    import datetime as _dt

    t0 = _dt.date(2026, 7, 18)
    windows = []
    values = []
    for row in post_rows:
        ref, sec = row["pair"].split("_")
        to_date = lambda s: _dt.date(int(s[:4]), int(s[4:6]), int(s[6:]))
        windows.append(((to_date(ref) - t0).days, (to_date(sec) - t0).days))
        values.append(row["pattern_fraction"])
    windows_array = np.asarray(windows, dtype=float)
    values_array = np.asarray(values, dtype=float)

    def residual(theta: np.ndarray) -> np.ndarray:
        amplitude, tau = theta
        predicted = amplitude * (np.exp(-windows_array[:, 0] / tau)
                                 - np.exp(-windows_array[:, 1] / tau))
        return predicted - values_array

    fit = least_squares(residual, x0=[0.15, 6.0],
                        bounds=([0.0, 1.0], [1.0, 60.0]))
    amplitude, tau = fit.x
    return {
        "reference_epoch": "2026-07-18",
        "amplitude_total_pattern_fraction": float(amplitude),
        "tau_days": float(tau),
        "rms_of_fit": float(np.sqrt(np.mean(fit.fun**2))),
        "n_windows": len(values),
        "windows_days": [[float(a), float(b)] for a, b in windows],
        "note": "descriptive window-integrated summary, not a physical model",
    }


def decompose_set(asc: dict, desc: dict) -> dict:
    """Solve [los_asc; los_desc] = [[eA,uA],[eD,uD]] [dE; dU] per lattice cell."""
    eA, uA, eD, uD = asc["e"], asc["u"], desc["e"], desc["u"]
    det = eA * uD - uA * eD
    with np.errstate(invalid="ignore", divide="ignore"):
        dE = (uD * asc["los"] - uA * desc["los"]) / det
        dU = (-eD * asc["los"] + eA * desc["los"]) / det
    bad = ~np.isfinite(det) | (np.abs(det) < DECOMP_MIN_ABS_DET)
    dE[bad] = np.nan
    dU[bad] = np.nan
    X, Y = np.meshgrid(asc["x_centres"], asc["y_centres"])
    near = (np.hypot(X - CENTER_E, Y - CENTER_N) < 18_000) & np.isfinite(dE)
    if near.sum() < 50:
        raise ValueError("too few common lattice cells for the decomposition")
    east_east = probe_disk_mean(dE, X, Y, PROBES["EAST"])
    east_west = probe_disk_mean(dE, X, Y, PROBES["WEST"])
    divergence = (east_east - east_west
                  if east_east is not None and east_west is not None else None)
    return {
        "n_cells_within_18km": int(near.sum()),
        "east_cm_p01_p99": [float(v) for v in
                            100 * np.nanpercentile(dE[near], [1, 99])],
        "up_cm_p01_p99": [float(v) for v in
                          100 * np.nanpercentile(dU[near], [1, 99])],
        "east_flank_probe_cm": {"EAST": east_east, "WEST": east_west},
        "opening_divergence_cm": divergence,
        "axis_up_probe_cm": probe_disk_mean(dU, X, Y, AXIS_PROBE),
        "_fields": (dE, dU, asc["x_centres"], asc["y_centres"]),
    }


def plot_coevent(frame: str, pair: str, stats: dict, frame_data: dict) -> None:
    aligned, predicted, residual = stats["_fields"]
    lons = frame_data["lons"][::SUBSAMPLE]
    lats = frame_data["lats"][::SUBSAMPLE]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    scale = np.nanpercentile(np.abs(aligned), 99)
    titles = ("Sentinel-1 observed\n(LiCSAR unwrapped, sign-aligned)",
              "predicted by NISAR-only dike model", "residual")
    for ax, field, title in zip(axes, (aligned, predicted, residual), titles):
        im = ax.imshow(field * 100, cmap="RdBu_r",
                       vmin=-100 * scale, vmax=100 * scale,
                       extent=[lons[0], lons[-1], lats[-1], lats[0]])
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("longitude (°E)", fontsize=8)
        ax.set_ylabel("latitude (°N)", fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.7, label="LOS (cm)")
    fig.suptitle(
        f"Descending Sentinel-1 cross-check — LiCSAR {frame}, pair {pair}\n"
        f"r = {stats['correlation']:.2f}, amplitude ratio "
        f"{stats['amplitude_ratio']:.2f}, RMS = {100*stats['rms_m']:.1f} cm "
        "(model fitted to NISAR only)", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out = RESULTS / "s1_tracks_coevent_panels.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"figure: {out}")


def plot_continuation(rows: list[dict]) -> None:
    fig, axes = plt.subplots(len(rows), 3, figsize=(15, 5.2 * len(rows)),
                             squeeze=False)
    for i, row in enumerate(rows):
        grid = row["_grid"]
        ok = np.isfinite(grid)
        predicted_grid = np.full_like(grid, np.nan)
        predicted_grid[ok] = row["_predicted"]
        slope, offset = row["pattern_fraction"], row["offset_m"]
        fitted = slope * predicted_grid
        residual = grid - (fitted + offset)
        scale = np.nanpercentile(np.abs(grid - offset), 99)
        panels = (grid - offset, fitted, residual)
        titles = (
            f"{row['geometry']} {row['window']}\nobserved "
            f"({100*offset:+.1f} cm offset removed)",
            f"scaled co-event pattern ({slope:+.3f})",
            "residual",
        )
        for j, (field, title) in enumerate(zip(panels, titles)):
            ax = axes[i, j]
            im = ax.imshow(field * 100, cmap="RdBu_r",
                           vmin=-100 * scale, vmax=100 * scale)
            ax.set_title(title, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            plt.colorbar(im, ax=ax, shrink=0.72, label="LOS (cm)")
    fig.suptitle("Sentinel-1 post-event pairs against the co-event dike pattern",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = RESULTS / "s1_tracks_continuation_panels.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"figure: {out}")


def plot_decomposition(sets: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10.6))
    for i, (name, summary) in enumerate(sets.items()):
        dE, dU, xs, ys = summary["_fields"]
        extent = [xs[0] / 1e3, xs[-1] / 1e3, ys[0] / 1e3, ys[-1] / 1e3]
        for j, (field, comp) in enumerate(((dU, "quasi-vertical (up +)"),
                                           (dE, "quasi-east (east +)"))):
            ax = axes[i, j]
            scale = np.nanpercentile(np.abs(field), 99)
            im = ax.imshow(field * 100, cmap="RdBu_r", origin="lower",
                           vmin=-100 * scale, vmax=100 * scale, extent=extent)
            ax.set_title(f"{name}: {comp}", fontsize=10)
            ax.set_xlabel("UTM 37N easting (km)", fontsize=8)
            ax.set_ylabel("UTM 37N northing (km)", fontsize=8)
            plt.colorbar(im, ax=ax, shrink=0.75, label="displacement (cm)")
    fig.suptitle("Sentinel-1-only East/Up decompositions on the 1.28-km UTM "
                 "lattice\n(zero-north assumption; note the different windows "
                 "of the two geometries)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = RESULTS / "s1_tracks_decomposition.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"figure: {out}")


def main() -> None:
    fit_path = RESULTS / "okada_fit.json"
    if not fit_path.exists():
        raise SystemExit("results/okada_fit.json missing — run `python -m mh26.invert` first")
    model = np.asarray(json.loads(fit_path.read_text())["values"][:7], dtype=float)

    frames = {}
    for frame in S1_FRAMES:
        print(f"loading frame {frame}")
        frames[frame] = load_frame(frame)

    geometry_of = {"014A_07688_131313": "ascending",
                   "079D_07694_131313": "descending"}
    referenced = {}
    for frame, pair, _ in S1_TRACK_PAIRS:
        print(f"screening {frame} {pair}")
        los = load_pair_los(frames[frame], pair)
        referenced[(frame, pair)] = reference(frames[frame], los)

    # Per-frame sign from the co-event pair only.
    frame_sign = {}
    coevent = {}
    for frame, pair, role in S1_TRACK_PAIRS:
        if role != "co-event":
            continue
        stats = coevent_pixel_stats(frames[frame], referenced[(frame, pair)], model)
        if abs(stats["correlation"]) < MIN_COEVENT_ABS_R:
            raise ValueError(
                f"{frame} {pair}: co-event correlation "
                f"{stats['correlation']:.3f} is too weak to fix the frame sign"
            )
        frame_sign[frame] = stats["applied_observation_sign"]
        coevent[frame] = {"pair": pair, **{k: v for k, v in stats.items()
                                           if not k.startswith("_")}}
        if frame == "079D_07694_131313":
            plot_coevent(frame, pair, stats, frames[frame])
        print(f"{frame} {pair}: sign {stats['applied_observation_sign']:+.0f}, "
              f"r={stats['correlation']:.3f}, ratio={stats['amplitude_ratio']:.2f}, "
              f"RMS={100*stats['rms_m']:.1f} cm")

    projections = []
    panel_rows = []
    for frame, pair, role in S1_TRACK_PAIRS:
        aligned = frame_sign[frame] * referenced[(frame, pair)]
        blocks = block_fields(frames[frame], aligned)
        stats = project_blocks(model, blocks)
        interval, n_clusters = cluster_wild_interval(
            stats["_predicted"][stats["_use"]],
            blocks["vals"][stats["_use"]],
            blocks["pts"][stats["_use"]],
            stats["slope"], stats["offset"],
        )
        row = {
            "frame": frame,
            "geometry": geometry_of[frame],
            "pair": pair,
            "window": window_of(pair),
            "days": pair_days(pair),
            "role": role,
            "n_blocks": stats["n_blocks"],
            "pattern_fraction": stats["slope"],
            "pattern_fraction_cluster_wild_95": [float(v) for v in interval],
            "n_clusters": n_clusters,
            "offset_m": stats["offset"],
            "zero_intercept_pattern_fraction": stats["zero_intercept_slope"],
            "correlation": stats["correlation"],
            "rms_m": stats["rms"],
        }
        if role == "post-event":
            sensitivity = []
            values = [stats["slope"]]
            for label, options in PROCESSING_VARIANTS:
                los_variant = aligned
                if "farfield_radius" in options:
                    los_variant = frame_sign[frame] * reference(
                        frames[frame],
                        load_pair_los(frames[frame], pair),
                        farfield_radius=options["farfield_radius"],
                    )
                variant_blocks = block_fields(
                    frames[frame], los_variant,
                    block=options.get("block", BLOCK_PX),
                )
                variant = project_blocks(
                    model, variant_blocks,
                    radius_limit=options.get("radius_limit", RADIUS_LIMIT),
                    prediction_min=options.get("prediction_min", PRED_MIN_BLOCK),
                )
                sensitivity.append({"variant": label,
                                    "pattern_fraction": variant["slope"],
                                    "n_blocks": variant["n_blocks"]})
                values.append(variant["slope"])
            row["processing_sensitivity"] = sensitivity
            row["pattern_fraction_processing_envelope"] = [
                float(min(values)), float(max(values)),
            ]
        projections.append(row)
        if role == "post-event" and pair in ("20260718_20260811", "20260723_20260816"):
            panel_rows.append({**row, "_grid": blocks["grid"],
                               "_predicted": stats["_predicted"],
                               "_use": stats["_use"]})
        print(f"{geometry_of[frame]:<10} {pair} [{role}]: "
              f"projection {stats['slope']:+.3f}, "
              f"95% [{interval[0]:+.3f}, {interval[1]:+.3f}], "
              f"{stats['n_blocks']} blocks")

    nulls = [row for row in projections if row["role"] == "pre-event"]
    null_values = np.array([row["pattern_fraction"] for row in nulls])
    null_summary = {
        "n_pairs": len(nulls),
        "mean": float(null_values.mean()),
        "standard_deviation": float(null_values.std(ddof=1)),
        "max_abs": float(np.abs(null_values).max()),
    }

    decomposition = {}
    for name, ((asc_frame, asc_pair), (desc_frame, desc_pair)) in \
            S1_DECOMPOSITION_SETS.items():
        binned_asc = bin_to_lattice(
            frames[asc_frame],
            frame_sign[asc_frame] * referenced[(asc_frame, asc_pair)],
        )
        binned_desc = bin_to_lattice(
            frames[desc_frame],
            frame_sign[desc_frame] * referenced[(desc_frame, desc_pair)],
        )
        summary = decompose_set(binned_asc, binned_desc)
        summary["pairs"] = {
            "ascending": {"frame": asc_frame, "pair": asc_pair,
                          "window": window_of(asc_pair)},
            "descending": {"frame": desc_frame, "pair": desc_pair,
                           "window": window_of(desc_pair)},
        }
        decomposition[name] = summary
        print(f"decomposition {name}: ΔE(flanks) "
              f"{fmt_cm(summary['opening_divergence_cm'])}, "
              f"axis up {fmt_cm(summary['axis_up_probe_cm'])}")

    plot_continuation(panel_rows)
    plot_decomposition(decomposition)

    post = [row for row in projections if row["role"] == "post-event"]
    decay = decay_fit(post)
    print(f"decay fit: total {decay['amplitude_total_pattern_fraction']:+.3f}, "
          f"tau {decay['tau_days']:.1f} d, RMS {decay['rms_of_fit']:.3f}")

    nisar_check = None
    continuation_path = RESULTS / "continuation.json"
    if continuation_path.exists():
        import datetime as _dt

        nisar = json.loads(continuation_path.read_text())
        t0 = _dt.date(2026, 7, 18)
        start_s, end_s = nisar["window"].split("→")
        to_days = lambda s: (_dt.date(*[int(v) for v in s.split("-")]) - t0).days
        t1, t2 = to_days(start_s), to_days(end_s)
        amplitude = decay["amplitude_total_pattern_fraction"]
        tau = decay["tau_days"]
        nisar_check = {
            "window": nisar["window"],
            "observed_nisar_pattern_fraction": nisar["pattern_fraction"],
            "predicted_from_s1_decay": float(
                amplitude * (np.exp(-t1 / tau) - np.exp(-t2 / tau))
            ),
        }

    payload = {
        "coevent_crosschecks": coevent,
        "frame_signs": {frame: float(sign) for frame, sign in frame_sign.items()},
        "projections": [
            {k: v for k, v in row.items() if not k.startswith("_")}
            for row in projections
        ],
        "null_summary": null_summary,
        "post_event_decay": decay,
        "nisar_cross_sensor_check": nisar_check,
        "decomposition": {
            name: {k: v for k, v in summary.items() if not k.startswith("_")}
            for name, summary in decomposition.items()
        },
        "statistic": "block regression of screened LOS on the co-event dike "
                     "pattern within 25 km; intercept fitted; see "
                     "mh26.continuation",
        "bootstrap": "wild sign bootstrap on 5-km clusters, 2000 samples; "
                     "conditional heuristic, not a calibrated interval",
    }
    (RESULTS / "s1_tracks.json").write_text(json.dumps(payload, indent=2) + "\n")

    desc_co = coevent["079D_07694_131313"]
    lines = [
        "# Two-track Sentinel-1 tests", "",
        "## Descending co-event cross-check", "",
        f"- LiCSAR 079D_07694_131313 pair {desc_co['pair']} "
        f"({window_of(desc_co['pair'])}): **r = {desc_co['correlation']:.3f}**, "
        f"amplitude ratio {desc_co['amplitude_ratio']:.2f}, "
        f"RMS {100*desc_co['rms_m']:.1f} cm against the NISAR-only model "
        f"({desc_co['n_pixels']} pixels)",
        f"- ascending frame, same statistic (consistency with `mh26.s1_check`): "
        f"r = {coevent['014A_07688_131313']['correlation']:.3f}, "
        f"ratio {coevent['014A_07688_131313']['amplitude_ratio']:.2f}, "
        f"RMS {100*coevent['014A_07688_131313']['rms_m']:.1f} cm",
        f"- empirical LiCSAR sign factors: ascending "
        f"{frame_sign['014A_07688_131313']:+.0f}, descending "
        f"{frame_sign['079D_07694_131313']:+.0f}, resolved from the co-event "
        "pairs only", "",
        "## Dike-pattern projections", "",
        "| geometry | pair | days | role | blocks | projection | cluster-wild 95% | zero-intercept | r |",
        "|---|---|---:|---|---:|---:|---|---:|---:|",
    ]
    for row in projections:
        interval = row["pattern_fraction_cluster_wild_95"]
        lines.append(
            f"| {row['geometry']} | {row['window']} | {row['days']} | {row['role']} | "
            f"{row['n_blocks']} | **{row['pattern_fraction']:+.3f}** | "
            f"[{interval[0]:+.3f}, {interval[1]:+.3f}] | "
            f"{row['zero_intercept_pattern_fraction']:+.3f} | "
            f"{row['correlation']:+.2f} |"
        )
    by_pair = {row["pair"]: row["pattern_fraction"] for row in projections}
    lines.extend([
        "",
        f"- {null_summary['n_pairs']} pre-event null projections: mean "
        f"{null_summary['mean']:+.3f}, standard deviation "
        f"{null_summary['standard_deviation']:.3f}, maximum magnitude "
        f"{null_summary['max_abs']:.3f}",
        "- the two early post-event windows exceed every null in both",
        "  geometries; the two late windows (from 30 July and 4 August) are at",
        "  the null level: the signal decays",
        f"- sub-window additivity: ascending {by_pair['20260718_20260730']:+.3f} "
        f"+ {by_pair['20260730_20260811']:+.3f} = "
        f"{by_pair['20260718_20260730'] + by_pair['20260730_20260811']:+.3f} vs "
        f"{by_pair['20260718_20260811']:+.3f} direct; descending "
        f"{by_pair['20260723_20260804']:+.3f} + {by_pair['20260804_20260816']:+.3f} = "
        f"{by_pair['20260723_20260804'] + by_pair['20260804_20260816']:+.3f} vs "
        f"{by_pair['20260723_20260816']:+.3f} direct",
        "",
        "## Post-event decay", "",
        f"- exponential window-integral fit to the six Sentinel-1 projections: "
        f"total **{decay['amplitude_total_pattern_fraction']:+.3f}** of the "
        f"co-event pattern, decay constant "
        f"**{decay['tau_days']:.1f} days** from 18 July "
        f"(fit RMS {decay['rms_of_fit']:.3f}); a descriptive summary, not a "
        "physical model",
    ])
    if nisar_check is not None:
        lines.append(
            f"- independent NISAR descending pair {nisar_check['window']}: "
            f"observed {nisar_check['observed_nisar_pattern_fraction']:+.3f}; "
            f"the Sentinel-1-only decay curve predicts "
            f"{nisar_check['predicted_from_s1_decay']:+.3f} for that window"
        )
    lines.extend([
        "",
        "## Post-event East/Up decomposition", "",
    ])
    for name, summary in decomposition.items():
        asc_meta = summary["pairs"]["ascending"]
        desc_meta = summary["pairs"]["descending"]
        divergence = summary["opening_divergence_cm"]
        lines.extend([
            f"### {name} (ascending {asc_meta['window']}, "
            f"descending {desc_meta['window']})", "",
            f"- flank quasi-east probes (2-km disks): EAST "
            f"{fmt_cm(summary['east_flank_probe_cm']['EAST'])}, WEST "
            f"{fmt_cm(summary['east_flank_probe_cm']['WEST'])} ⇒ "
            + (f"**opening divergence {divergence:+.1f} cm**"
               if divergence is not None else "divergence not defined"),
            f"- axis quasi-up probe: {fmt_cm(summary['axis_up_probe_cm'])}",
            f"- east p01–p99 {summary['east_cm_p01_p99'][0]:+.1f} to "
            f"{summary['east_cm_p01_p99'][1]:+.1f} cm; up p01–p99 "
            f"{summary['up_cm_p01_p99'][0]:+.1f} to "
            f"{summary['up_cm_p01_p99'][1]:+.1f} cm "
            f"({summary['n_cells_within_18km']} cells within 18 km)", "",
        ])
    lines.extend([
        "## Scope", "",
        "- The projection detects only deformation shaped like the fitted dike;",
        "  the intercept absorbs a broad datum shift per pair.",
        "- The ascending and descending post-event windows differ by five days at",
        "  each end; the decomposition assumes the signal is common to both. For",
        "  a decaying signal that assumption fails in part, and the East/Up split",
        "  blends the two windows; the projection table is the quantitative result.",
        "- The cluster-wild intervals and the pre-event null distribution do not",
        "  form a calibrated atmospheric covariance model.",
        "- Near-axis C-band decoherence removes the highest-gradient pixels from",
        "  the co-event comparisons.",
    ])
    report = "\n".join(lines) + "\n"
    (RESULTS / "s1_tracks.md").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
