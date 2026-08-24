"""Okada dike inversion, joint ascending + descending NISAR.

Single rectangular tensile dislocation (vertical dike), uniform opening,
solved by multistart least squares on block-averaged, ionosphere-corrected,
far-field-referenced and valid-sample-masked LOS from both event geometries,
with a per-geometry offset nuisance parameter.

Run from the repo root:  python -m mh26.invert
(Run `python -m mh26.okada` once first to confirm sign conventions.)

Writes results/okada_fit.md, results/okada_fit.json (read by mh26.s1_check),
and results/okada_fit_panels.png.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from .gunw import block_average, corrected_los, load_pair, los_unit_neu, reference_farfield
from .okada import forward
from .params import (
    BLOCK_PIXELS,
    CENTER_E,
    CENTER_N,
    FARFIELD_RADIUS,
    GRID_POSTING,
    INVERSION_PAIRS,
    LOS_HEIGHT_AE,
    RESULTS,
    SHEAR_MODULUS,
)

PARAM_NAMES = ["east_c", "north_c", "strike", "length", "top", "bottom", "opening",
               "offs_asc", "offs_desc"]

# Bounds and multistart grid for [e_c, n_c, strike, length, top, bottom,
# opening, offs_asc, offs_desc]; metres/degrees.
BOUNDS_LO = [660000, 1330000, 120, 5000, 0, 1500, 0.0, -0.2, -0.2]
BOUNDS_HI = [705000, 1365000, 185, 60000, 4000, 12000, 8.0, 0.2, 0.2]
X_SCALE = [1e3, 1e3, 5, 1e3, 200, 500, 0.2, 0.01, 0.01]
X_ORIGIN = [683000, 1350000, 150, 18000, 800, 7000, 1.7, 0.0, 0.0]
STARTS = [[684000, 1348000, strike, length, 500, 6000, 2.0, 0.0, 0.0]
          for strike in (150, 165, 135) for length in (20000, 35000)]
FINITE_DIFF_STEP = 1e-3


def load_observations(block: int = BLOCK_PIXELS,
                      farfield_radius: float = FARFIELD_RADIUS,
                      los_height: float = LOS_HEIGHT_AE,
                      max_radius: float | None = None) -> dict:
    """Load the two inversion geometries with configurable sensitivity choices."""
    data = {}
    for tag in INVERSION_PAIRS:
        d = load_pair(tag)
        if d is None:
            raise SystemExit(f"{tag}.npz missing from data/nisar — see data/README.md")
        los, _ = corrected_los(d)
        ref = reference_farfield(los, d["x"], d["y"], r_in=farfield_radius)
        pts, vals, grid = block_average(ref, d["x"], d["y"], block=block)
        if max_radius is not None:
            radius = np.hypot(pts[:, 1] - CENTER_E, pts[:, 0] - CENTER_N)
            use = radius < max_radius
            pts, vals = pts[use], vals[use]
        units = los_unit_neu(d, pts, height=los_height)
        data[tag] = dict(pts=pts, vals=vals, units=units, grid=grid)
        print(f"{tag}: {len(vals)} block obs, mean local LOS unit (n,e,up)="
              f"{units.mean(axis=0).round(3)}")
    return data


def residual_vector(theta: np.ndarray, data: dict,
                    equal_geometry_weight: bool = False) -> np.ndarray:
    """Model-minus-observed objective vector for one physical parameter vector."""
    p, offsets = theta[:7], theta[7:]
    out = []
    for i, tag in enumerate(INVERSION_PAIRS):
        dd = data[tag]
        displacement = forward(p, dd["pts"])
        los_pred = np.sum(displacement * dd["units"], axis=1)
        residual = los_pred + offsets[i] - dd["vals"]
        if equal_geometry_weight:
            residual = residual / np.sqrt(len(residual))
        out.append(residual)
    return np.concatenate(out)


def fit(data: dict, starts: list[list[float]] | list[np.ndarray] | None = None,
        equal_geometry_weight: bool = False):
    """Fit the dike using numerically scaled variables and report diagnostics.

    The former implementation passed metre-scale and unit-scale parameters
    directly to finite differencing.  It could terminate on ``xtol`` with a
    large gradient and a start-dependent model.  Optimising dimensionless
    variables with explicit finite-difference steps makes all six default
    starts converge to the same basin.
    """
    scale = np.asarray(X_SCALE, dtype=float)
    origin = np.asarray(X_ORIGIN, dtype=float)
    lower = (np.asarray(BOUNDS_LO, dtype=float) - origin) / scale
    upper = (np.asarray(BOUNDS_HI, dtype=float) - origin) / scale
    starts = STARTS if starts is None else starts

    results = []
    solutions = []
    for s0 in starts:
        q0 = (np.asarray(s0, dtype=float) - origin) / scale
        try:
            result = least_squares(
                lambda q: residual_vector(
                    origin + scale * q, data,
                    equal_geometry_weight=equal_geometry_weight,
                ),
                q0,
                bounds=(lower, upper),
                x_scale="jac",
                diff_step=FINITE_DIFF_STEP,
                max_nfev=2000,
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            )
        except ValueError as error:
            print("start failed:", error)
            continue
        results.append(result)
        solutions.append(origin + scale * result.x)
    if not results:
        raise SystemExit("all inversion starts failed")

    best_index = int(np.argmin([result.cost for result in results]))
    best = results[best_index]
    theta = solutions[best_index]
    if not best.success or best.optimality > 0.1:
        raise RuntimeError(
            "inversion did not converge cleanly: "
            f"success={best.success}, optimality={best.optimality:.3g}, "
            f"message={best.message}"
        )

    # Always report physical (unweighted) residuals and RMS.
    res = residual_vector(theta, data, equal_geometry_weight=False)
    n_obs = len(res)
    rms = float(np.sqrt(np.mean(res**2)))
    # White-noise-only diagnostic covariance.  The Jacobian is with respect to
    # scaled variables, so transform the resulting sigmas back to physical units.
    jacobian = best.jac
    try:
        covariance_q = np.linalg.inv(jacobian.T @ jacobian)
        covariance_q *= 2 * best.cost / max(n_obs - len(theta), 1)
        sig = scale * np.sqrt(np.diag(covariance_q))
    except np.linalg.LinAlgError:
        sig = np.full(len(theta), np.nan)

    solution_array = np.asarray(solutions)
    diagnostics = {
        "starts_attempted": len(starts),
        "starts_completed": len(results),
        "best_cost": float(best.cost),
        "status": int(best.status),
        "message": best.message,
        "nfev": int(best.nfev),
        "optimality": float(best.optimality),
        "max_start_solution_spread": [
            float(value) for value in np.ptp(solution_array, axis=0)
        ],
        "equal_geometry_weight": equal_geometry_weight,
    }
    return theta, sig, rms, n_obs, diagnostics


def write_report(p: np.ndarray, sig: np.ndarray, rms: float, n_obs: int,
                 diagnostics: dict) -> None:
    lines = ["# Okada dike fit — Manda Hararo July 2026", "",
             f"- observations: {n_obs} (block-avg {BLOCK_PIXELS * GRID_POSTING:.0f} m), "
             f"RMS misfit **{100*rms:.1f} cm**",
             "- screening is entirely a priori: coherence, connected component, and the",
             "  product valid-sample mask; no blocks were rejected using model residuals",
             f"- optimizer: {diagnostics['starts_completed']}/{diagnostics['starts_attempted']} "
             f"starts completed; best optimality {diagnostics['optimality']:.3g}",
             "", "| param | value | formal 1σ* | unit |", "|---|---|---|---|"]
    for name, v, s in zip(PARAM_NAMES, p, sig):
        unit = "deg" if name == "strike" else "m"
        prec = 4 if name.startswith("offs_") else (2 if name == "opening" else 1)
        lines.append(f"| {name} | {v:.{prec}f} | ±{s:.{prec}f} | {unit} |")
    lines.extend([
        "",
        "*The formal errors are a numerical white-noise diagnostic only. InSAR noise is",
        "spatially correlated and the uniform-opening half-space is an approximation, so",
        "these values are **not** parameter confidence intervals. Processing-choice",
        "sensitivity is reported separately in `model_sensitivity.md`.",
    ])
    length, width, opening = p[3], p[5] - p[4], p[6]
    volume = length * width * opening
    moment = SHEAR_MODULUS * volume  # tensile geodetic moment scale
    lines.append(f"\n- dike area {length/1e3:.1f} km × {width/1e3:.1f} km, "
                 f"opening {opening:.2f} m ⇒ **geometric volume ≈ {volume/1e9:.2f} km³**")
    lines.append(f"- geodetic moment scale μV ≈ {moment:.2e} N·m "
                 f"(Mw≈{(2/3)*(np.log10(moment)-9.1):.1f} only as a μV comparison; "
                 "this is not a measured seismic moment)")
    md = "\n".join(lines)
    print(md)
    (RESULTS / "okada_fit.md").write_text(md + "\n")
    (RESULTS / "okada_fit.json").write_text(json.dumps({
        "names": PARAM_NAMES,
        "values": [float(v) for v in p],
        "formal_white_noise_sigma": [float(s) for s in sig],
        "rms_m": rms,
        "n_obs": n_obs,
        "volume_km3": volume / 1e9,
        "geodetic_moment_Nm": moment,
        "screening": "coherence + connected component + GUNW valid-sample mask",
        "unit_vectors": "cube-interpolated @500m",
        "optimizer": diagnostics,
    }, indent=2) + "\n")


def plot_panels(data: dict, p: np.ndarray) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.6))
    row_names = {"172_008_A": "NISAR ascending 07-08→07-20",
                 "165_082_D": "NISAR descending 07-07→07-19"}
    for i, tag in enumerate(INVERSION_PAIRS):
        dd = data[tag]
        xs, ys, m = dd["grid"]
        displacement = forward(p[:7], dd["pts"])
        los_pred = np.sum(displacement * dd["units"], axis=1) + p[7 + i]
        pred_grid = np.full_like(m, np.nan)
        pred_grid[np.isfinite(m)] = los_pred
        resid = m - pred_grid
        v = np.nanpercentile(np.abs(m), 99)
        for j, (F, title) in enumerate(((m, f"{row_names.get(tag, tag)}\nobserved"),
                                        (pred_grid, "dike model"),
                                        (resid, "residual"))):
            ax = axes[i, j]
            im = ax.imshow(F * 100, cmap="RdBu_r", vmin=-100 * v, vmax=100 * v,
                           extent=[xs[0] / 1e3, xs[-1] / 1e3, ys[-1] / 1e3, ys[0] / 1e3])
            ax.set_title(title, fontsize=10)
            ax.set_xlabel("UTM 37N easting (km)", fontsize=8)
            if j == 0:
                ax.set_ylabel("UTM 37N northing (km)", fontsize=8)
            plt.colorbar(im, ax=ax, shrink=0.7, label="LOS (cm)")
    length, width, opening = p[3], p[5] - p[4], p[6]
    fig.suptitle(f"Okada dike fit, Manda Hararo — strike {p[2]:.0f}°, "
                 f"{length/1e3:.1f} km × {p[4]/1e3:.1f}–{p[5]/1e3:.1f} km depth, "
                 f"opening {opening:.2f} m ({length*width*opening/1e9:.2f} km³)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = RESULTS / "okada_fit_panels.png"
    fig.savefig(out, dpi=110)
    print(f"figure: {out}")


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    data = load_observations()
    p, sig, rms, n_obs, diagnostics = fit(data)
    write_report(p, sig, rms, n_obs, diagnostics)
    plot_panels(data, p)


if __name__ == "__main__":
    main()
