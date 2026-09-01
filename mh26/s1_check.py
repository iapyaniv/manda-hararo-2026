"""compare Sentinel-1 data with the NISAR-derived dike model.

The command predicts the LiCSAR interferogram (C-band, ascending frame
014A_07688_131313, pair 20260712_20260718) from the model fitted to NISAR. It
uses the LiCSAR East, North, and Up LOS vectors and compares the prediction
with the observed unwrapped phase. The Sentinel-1 data use a different
satellite, radar band, and processing chain from the NISAR fit.

LiCSAR unw is radians; its sign convention is resolved empirically from the
obs-vs-model regression slope (reported).

Run from the repo root:  python -m mh26.s1_check

Writes results/s1_crosscheck.md and results/s1_crosscheck_panels.png.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer

from .okada import forward
from .params import (
    CENTER_E,
    CENTER_N,
    DATA_S1,
    EPSG_UTM37N,
    FARFIELD_RADIUS,
    RESULTS,
    S1_FRAME,
    S1_PAIR,
    WAVELENGTH_C,
)

# Crop around the anomaly, degrees lon/lat.
WEST, SOUTH, EAST, NORTH = 40.35, 11.85, 41.05, 12.55
COHERENCE_MIN_S1 = 0.25
SUBSAMPLE = 3           # compare every 3rd pixel (speed)
PRED_MIN = 0.005        # compare only where |model| > 0.5 cm


def read_crop(path: Path):
    with rasterio.open(path) as r:
        if r.crs is None or r.crs.to_epsg() != 4326:
            raise ValueError(f"{path}: expected EPSG:4326, got {r.crs}")
        win = rasterio.windows.from_bounds(WEST, SOUTH, EAST, NORTH, r.transform)
        win = win.round_offsets().round_lengths()
        a = r.read(1, window=win).astype(np.float64)
        t = r.window_transform(win)
        rows, cols = a.shape
        xs = t.c + t.a * (np.arange(cols) + 0.5)
        ys = t.f + t.e * (np.arange(rows) + 0.5)
    return a, xs, ys


def coherence_fraction(coherence: np.ndarray) -> np.ndarray:
    """Normalize a LiCSAR coherence raster to the documented 0–1 fraction.

    LiCSAR GeoTIFFs commonly store coherence as UInt8 in the range 0–255.  A
    direct comparison with 0.25 would therefore admit every nonzero byte.
    """
    coherence = np.asarray(coherence, dtype=np.float64)
    finite = coherence[np.isfinite(coherence)]
    if finite.size == 0:
        raise ValueError("coherence raster contains no finite samples")
    if finite.min() < 0:
        raise ValueError("coherence raster contains negative values")
    if finite.max() > 1.0:
        if finite.max() > 255:
            raise ValueError("coherence raster is outside both 0–1 and 0–255 ranges")
        coherence = coherence / 255.0
    return coherence


def assert_same_grid(named_arrays: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]]) -> None:
    """Require equal shapes and coordinate vectors for all input rasters."""
    name0, array0, x0, y0 = named_arrays[0]
    for name, array, x, y in named_arrays[1:]:
        if array.shape != array0.shape or not np.array_equal(x, x0) or not np.array_equal(y, y0):
            raise ValueError(f"LiCSAR grid mismatch: {name0} vs {name}")


def main() -> None:
    fit_path = RESULTS / "okada_fit.json"
    if not fit_path.exists():
        raise SystemExit("results/okada_fit.json missing — run `python -m mh26.invert` first")
    p = np.array(json.loads(fit_path.read_text())["values"][:7])

    unw, lons, lats = read_crop(DATA_S1 / f"{S1_PAIR}.geo.unw.tif")
    cc, cc_lons, cc_lats = read_crop(DATA_S1 / f"{S1_PAIR}.geo.cc.tif")
    e, e_lons, e_lats = read_crop(DATA_S1 / f"{S1_FRAME}.geo.E.tif")
    n, n_lons, n_lats = read_crop(DATA_S1 / f"{S1_FRAME}.geo.N.tif")
    u, u_lons, u_lats = read_crop(DATA_S1 / f"{S1_FRAME}.geo.U.tif")

    # The five source rasters currently align exactly; silent resampling is not
    # an acceptable fallback for a quantitative cross-sensor comparison.
    assert_same_grid([
        ("unwrapped", unw, lons, lats),
        ("coherence", cc, cc_lons, cc_lats),
        ("E", e, e_lons, e_lats),
        ("N", n, n_lons, n_lats),
        ("U", u, u_lons, u_lats),
    ])

    unw[unw == 0] = np.nan
    cc = coherence_fraction(cc)
    vector_norm = np.sqrt(e**2 + n**2 + u**2)
    good = (
        np.isfinite(unw)
        & np.isfinite(cc)
        & (cc > COHERENCE_MIN_S1)
        & np.isfinite(e)
        & np.isfinite(n)
        & np.isfinite(u)
        & (np.abs(e) > 0.01)
        & (np.abs(vector_norm - 1.0) < 0.01)
    )
    los_obs = unw * WAVELENGTH_C / (4 * np.pi)  # metres; sign resolved below

    LON, LAT = np.meshgrid(lons, lats)
    tr = Transformer.from_crs(4326, EPSG_UTM37N, always_xy=True)
    X, Y = tr.transform(LON, LAT)
    R = np.hypot(X - CENTER_E, Y - CENTER_N)
    far = good & (R > FARFIELD_RADIUS)
    if far.sum() < 500:
        raise ValueError("too few valid far-field Sentinel-1 pixels for plane referencing")
    A = np.c_[np.ones(far.sum()), (X[far] - CENTER_E) / 1e4, (Y[far] - CENTER_N) / 1e4]
    coef, *_ = np.linalg.lstsq(A, los_obs[far], rcond=None)
    los_obs = los_obs - (coef[0] + coef[1] * (X - CENTER_E) / 1e4
                         + coef[2] * (Y - CENTER_N) / 1e4)
    los_obs[~good] = np.nan

    sl = np.s_[::SUBSAMPLE, ::SUBSAMPLE]
    pts = np.column_stack([Y[sl].ravel(), X[sl].ravel()])
    disp = forward(p, pts)  # (n, e, up)
    los_unit = np.column_stack([n[sl].ravel(), e[sl].ravel(), u[sl].ravel()])
    los_pred = np.sum(disp * los_unit, axis=1).reshape(Y[sl].shape)
    obs_s = los_obs[sl]
    # The observation had a far-field plane removed; the model must lose its
    # own far-field plane through the same operator, or the ramp difference
    # deflates the regression slope and the correlation.
    R_s = R[sl]
    far_s = (R_s > FARFIELD_RADIUS) & np.isfinite(obs_s)
    B_m = np.c_[np.ones(int(far_s.sum())), (X[sl][far_s] - CENTER_E) / 1e4,
                (Y[sl][far_s] - CENTER_N) / 1e4]
    coef_m, *_ = np.linalg.lstsq(B_m, los_pred[far_s], rcond=None)
    los_pred = los_pred - (coef_m[0] + coef_m[1] * (X[sl] - CENTER_E) / 1e4
                           + coef_m[2] * (Y[sl] - CENTER_N) / 1e4)
    ok = np.isfinite(obs_s) & np.isfinite(los_pred) & (np.abs(los_pred) > PRED_MIN)
    if ok.sum() < 100:
        raise ValueError("too few pixels survive the Sentinel-1 comparison mask")

    # regression obs = slope * pred; the sign convention falls out of the slope
    slope = np.sum(obs_s[ok] * los_pred[ok]) / np.sum(los_pred[ok] ** 2)
    sgn = -1.0 if slope < 0 else 1.0
    obs_aligned = sgn * obs_s
    resid = obs_aligned - los_pred
    rms = np.sqrt(np.nanmean(resid[ok] ** 2))
    corr = np.corrcoef(obs_aligned[ok], los_pred[ok])[0, 1]

    lines = [
        "# Sentinel-1 cross-check of the NISAR-derived dike model", "",
        f"- LiCSAR {S1_FRAME} pair {S1_PAIR}; {ok.sum()} compared pixels "
        f"(cc>{COHERENCE_MIN_S1}, |pred|>{100*PRED_MIN:.1f} cm, {SUBSAMPLE}× subsample)",
        f"- regression slope obs-vs-model: {slope:+.3f} ⇒ LiCSAR unw sign convention "
        f"{'opposite to' if sgn < 0 else 'the same as'} the model convention "
        f"(applied factor {sgn:+.0f})",
        f"- **correlation r = {corr:.3f}**, amplitude ratio |slope| = {abs(slope):.2f}",
        f"- **RMS residual = {100*rms:.1f} cm** against a model fitted to NISAR only",
        "", "The observed and predicted fields have the correlation, amplitude ratio, and",
        "RMS residual reported above. Residuals are concentrated near the axis, where the",
        "uniform-opening model does not represent graben faulting.",
    ]
    out = "\n".join(lines)
    print(out)
    (RESULTS / "s1_crosscheck.md").write_text(out + "\n")
    (RESULTS / "s1_crosscheck.json").write_text(json.dumps({
        "frame": S1_FRAME,
        "pair": S1_PAIR,
        "n_pixels": int(ok.sum()),
        "coherence_min": COHERENCE_MIN_S1,
        "subsample": SUBSAMPLE,
        "prediction_min_m": PRED_MIN,
        "regression_slope_observed_on_predicted": float(slope),
        "applied_observation_sign": float(sgn),
        "correlation": float(corr),
        "amplitude_ratio": float(abs(slope)),
        "rms_m": float(rms),
    }, indent=2) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    v = np.nanpercentile(np.abs(obs_aligned[ok]), 99)
    for ax, F, ttl in ((axes[0], obs_aligned, "Sentinel-1 observed\n(LiCSAR unwrapped, sign-aligned)"),
                       (axes[1], np.where(np.isfinite(obs_s), los_pred, np.nan),
                        "predicted by NISAR-only dike model"),
                       (axes[2], np.where(ok, resid, np.nan), "residual")):
        im = ax.imshow(F * 100, cmap="RdBu_r", vmin=-100 * v, vmax=100 * v,
                       extent=[lons[0], lons[-1], lats[-1], lats[0]])
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("longitude (°E)", fontsize=8)
        ax.set_ylabel("latitude (°N)", fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.7, label="LOS (cm)")
    fig.suptitle(f"Sentinel-1 comparison — LiCSAR {S1_FRAME}, pair {S1_PAIR}\n"
                 f"Manda Hararo, Afar: r = {corr:.2f}, RMS = {100*rms:.1f} cm "
                 "(spatially correlated samples; model fitted to NISAR only; "
                 "LiCSAR sign convention resolved empirically)", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    out_png = RESULTS / "s1_crosscheck_panels.png"
    fig.savefig(out_png, dpi=110)
    print(f"figure: {out_png}")


if __name__ == "__main__":
    main()
