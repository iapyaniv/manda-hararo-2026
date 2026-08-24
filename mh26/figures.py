"""publication-style fringe comparison of data and model.

Rewrapped-LOS panels in the presentation of Xu et al. (2026,
doi:10.1016/j.geog.2026.06.002): each colour cycle is 10 cm of LOS
displacement toward/away from the satellite, so fringe counts read directly
as displacement, and observed and modelled fields can be compared fringe by
fringe. Rows: the two reprocessed NISAR geometries (30 m, decay-normalized
as in the inversion); columns: observed, model prediction, continuous
residual; bottom: a rift-perpendicular profile through the maximum-opening
section.

The model is the distributed fit to the GUNW inputs with the idealized
fault geometry (results/distributed_fit.json) — the primary reported model
— evaluated on the reprocessed 30-m grids, so the figure is a
cross-processing comparison, not a self-fit. It is a qualitative
fringe-level check; the residual panels are clipped at ±10 cm and keep
near-axis residuals beyond that.

Run after `mh26.distributed` and `mh26.gslc --unwrap`:
    python -m mh26.figures

Writes results/fringes_model_comparison.png.
"""

from __future__ import annotations

import datetime as dt
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyrocko.modelling import okada_ext

from .distributed import decay_factor, geometry_from_json, geometry_sha256
from .gslc import DATA_GSLC, UNWRAP_COHERENCE_MIN
from .gunw import load_pair, los_unit_neu, pair_window
from .params import (
    CENTER_E,
    CENTER_N,
    FARFIELD_RADIUS,
    INVERSION_PAIRS,
    LAME_LAMBDA,
    RESULTS,
    SHEAR_MODULUS,
)

CYCLE_M = 0.10          # one colour cycle of rewrapped LOS
SUBSAMPLE = 5           # 30-m pixels -> 150-m model evaluation
PROFILE_S_M = 1_000.0   # along-strike position of the profile
PROFILE_HALF_M = 20_000.0
PROFILE_BIN_M = 500.0
LABEL = {"172_008_A": "NISAR ascending 07-08→07-20",
         "165_082_D": "NISAR descending 07-07→07-19"}


def rewrap(los: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(2j * np.pi * los / CYCLE_M))


def main() -> None:
    fit = json.loads((RESULTS / "distributed_fit.json").read_text())
    decay = json.loads((RESULTS / "s1_tracks.json").read_text())["post_event_decay"]

    geometry = geometry_from_json(fit["geometry"])
    if fit.get("geometry_sha256") != geometry_sha256(geometry):
        raise ValueError("distributed_fit.json geometry hash does not match")
    magnitudes = np.asarray(fit["dislocation_magnitudes_m"], dtype=float)
    if magnitudes.shape != geometry["kinds"].shape:
        raise ValueError("saved dislocations do not match the saved geometry")
    disl = geometry["disl"] * magnitudes[:, None]

    fig = plt.figure(figsize=(16, 13))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.55])
    profiles = {}

    for row, tag in enumerate(INVERSION_PAIRS):
        f = np.load(DATA_GSLC / f"{tag}_unw.npz")
        los = np.where(f["coherence"] > UNWRAP_COHERENCE_MIN, f["los"], np.nan)
        d = load_pair(tag)
        window = pair_window(d)
        end = dt.date(*[int(v) for v in window.split("→")[1].split("-")])
        los = los / decay_factor(end, decay)

        sl = np.s_[::SUBSAMPLE, ::SUBSAMPLE]
        x, y = f["x"][::SUBSAMPLE], f["y"][::SUBSAMPLE]
        obs = los[sl]
        X, Y = np.meshgrid(x, y)
        pts = np.column_stack([Y.ravel(), X.ravel()])
        rec = np.zeros((len(pts), 3))
        rec[:, 0], rec[:, 1] = pts[:, 0], pts[:, 1]
        res = okada_ext.okada(geometry["patches"], disl, rec, LAME_LAMBDA,
                              SHEAR_MODULUS, nthreads=0, rotate_sdn=0,
                              stack_sources=1)
        disp = np.column_stack([res[:, 0], res[:, 1], -res[:, 2]])
        units = los_unit_neu(d, pts)
        pred = np.sum(disp * units, axis=1).reshape(obs.shape)
        far = (np.hypot(X - CENTER_E, Y - CENTER_N) > FARFIELD_RADIUS) \
            & np.isfinite(obs)
        pred = pred + np.nanmedian(obs[far] - pred[far])
        resid = obs - pred

        extent = [x[0] / 1e3, x[-1] / 1e3, y[-1] / 1e3, y[0] / 1e3]
        for col, (field, title, wrapped) in enumerate((
                (obs, f"{LABEL[tag]}\nobserved (rewrapped, "
                      f"{100 * CYCLE_M:.0f} cm/cycle)", True),
                (np.where(np.isfinite(obs), pred, np.nan),
                 "distributed model (fit to GUNW inputs)", True),
                (resid, "residual (clipped ±10 cm)", False))):
            ax = fig.add_subplot(gs[row, col])
            if wrapped:
                im = ax.imshow(rewrap(field), cmap="twilight_shifted",
                               vmin=-np.pi, vmax=np.pi, extent=extent)
                cbar = plt.colorbar(im, ax=ax, shrink=0.75)
                cbar.set_ticks([-np.pi, 0, np.pi])
                cbar.set_ticklabels(["0", f"{50 * CYCLE_M:.0f}",
                                     f"{100 * CYCLE_M:.0f} cm"])
            else:
                im = ax.imshow(field * 100, cmap="RdBu_r", vmin=-10, vmax=10,
                               extent=extent)
                plt.colorbar(im, ax=ax, shrink=0.75, label="cm")
            ax.set_title(title, fontsize=10)
            ax.set_xlabel("UTM 37N easting (km)", fontsize=8)
            if col == 0:
                ax.set_ylabel("UTM 37N northing (km)", fontsize=8)

        # Rift-perpendicular profile through the maximum-opening section
        # (signed strike-frame coordinates, perp positive ENE as in
        # mh26.faults.strike_frame).
        e_c, n_c = geometry["trace_centre"]
        s_r = np.radians(geometry["strike"])
        signed_perp = ((Y - n_c) * np.cos(s_r - np.pi / 2)
                       + (X - e_c) * np.sin(s_r - np.pi / 2))
        along_signed = ((Y - n_c) * np.cos(s_r)
                        + (X - e_c) * np.sin(s_r))
        band = np.abs(along_signed - PROFILE_S_M) < 1_000.0
        edges = np.arange(-PROFILE_HALF_M, PROFILE_HALF_M + 1, PROFILE_BIN_M)
        centres = 0.5 * (edges[:-1] + edges[1:])
        prof_obs, prof_mod = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            cell = band & (signed_perp >= lo) & (signed_perp < hi)
            prof_obs.append(np.nanmedian(obs[cell]) if np.isfinite(obs[cell]).sum() > 20
                            else np.nan)
            prof_mod.append(np.nanmedian(pred[cell]) if cell.sum() > 20 else np.nan)
        profiles[tag] = (centres, np.array(prof_obs), np.array(prof_mod))

    ax = fig.add_subplot(gs[2, :])
    colours = {"172_008_A": "tab:red", "165_082_D": "tab:blue"}
    for tag, (centres, obs_p, mod_p) in profiles.items():
        ax.plot(centres / 1e3, 100 * obs_p, color=colours[tag], lw=1.5,
                label=f"{LABEL[tag]} observed")
        ax.plot(centres / 1e3, 100 * mod_p, color=colours[tag], lw=1.2,
                ls="--", label="model")
    ax.axvline(0, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel("distance from the dike trace, ENE positive (km)")
    ax.set_ylabel("LOS (cm)")
    ax.set_title(f"rift-perpendicular profile at +{PROFILE_S_M / 1e3:.0f} km "
                 "along strike (500-m bins)", fontsize=10)
    ax.legend(fontsize=8, ncols=2)
    ax.grid(alpha=0.3)

    fig.suptitle("Reprocessed 30-m NISAR interferograms against the "
                 "distributed dike + graben-fault model fitted to the GUNW "
                 "inputs (idealized fault geometry)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = RESULTS / "fringes_model_comparison.png"
    fig.savefig(out, dpi=110)
    print(f"figure: {out}")


if __name__ == "__main__":
    main()
