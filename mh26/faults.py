"""phase-gradient ridges flanking the rift axis.

Phase gradients of the co-event interferograms enhance localized high-strain
features (Xu et al., 2026, Geodesy and Geodynamics,
doi:10.1016/j.geog.2026.06.002, delineate surface ruptures from
phase-gradient maps; this module automates a narrower task and is not a
reimplementation of their mapping). A shallow fault break appears as a
narrow phase-gradient ridge, but a smooth steep gradient band can also
produce one, so the picker output is a candidate scarp line, not a mapped
fault:

1. wrapped-phase gradients from the 30-m reprocessed interferograms (the
   wrapped gradient needs no unwrapping and is immune to ambiguity errors);
2. per 500-m along-strike bin and per side of the axis, the strongest
   coherent gradient ridge between 0.3 and 4 km from the dike trace — by
   construction at most one ridge per side and bin;
3. ridges detected independently in both NISAR geometries and merged;
4. a LOS-step diagnostic across each ridge from the unwrapped fields: a
   surface-breaking normal fault with the axis side down gives a localized
   negative step (axis side minus outer side, local trend removed) with the
   same sign in both geometries; a smooth gradient band gives a step near
   zero.

The ridge lines feed the ridge-anchored inversion (mh26.distributed
--mapped) and the stress-change receivers (mh26.stress); their
interpretation as graben-bounding fault scarps rests on the cross-geometry
agreement and the step diagnostic, and remains unconfirmed against an
independent scarp map.

Run after `mh26.gslc --ifg --unwrap`:  python -m mh26.faults

Writes results/fault_traces.json and results/fault_traces.png.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .params import INVERSION_PAIRS, RESULTS

COHERENCE_MIN_GRADIENT = 0.3
S_HALF = 14_000.0        # m, along-strike extent searched
S_BIN = 500.0            # m
P_MIN, P_MAX = 300.0, 4_000.0   # m, ridge search band each side of the axis
P_BIN = 100.0            # m
MIN_SAMPLES = 40         # pixels per (s, p) cell
QUALITY_MIN = 1.8        # ridge peak over band median
QUALITY_SENSITIVITY = (1.5, QUALITY_MIN, 2.1)
COHERENCE_SENSITIVITY = (0.25, COHERENCE_MIN_GRADIENT, 0.35)
SMOOTH_BINS = 5          # rolling-median window along strike
STEP_INNER_M = 100.0     # m, exclusion half-width around the ridge
STEP_BAND_M = 400.0      # m, outer edge of the step-measurement bands
STEP_TREND_M = 1_200.0   # m, outer edge of the local trend fit


def wrapped_gradient(phase: np.ndarray, spacing: float) -> np.ndarray:
    """|grad(phi)| from wrapped one-sided differences, rad/m."""
    gx = np.full(phase.shape, np.nan, dtype=np.float32)
    gy = np.full(phase.shape, np.nan, dtype=np.float32)
    dx = np.angle(np.exp(1j * np.diff(phase, axis=1)))
    dy = np.angle(np.exp(1j * np.diff(phase, axis=0)))
    gx[:, :-1] = dx / spacing
    gy[:-1, :] = dy / spacing
    return np.hypot(gx, gy)


def strike_frame(X: np.ndarray, Y: np.ndarray, model: np.ndarray):
    """(along, perp) metres relative to the dike trace; perp positive ENE.

    The perpendicular sign must match mh26.distributed._strike_vectors
    (unit vector at azimuth strike − 90°, ENE for this dike): a flipped
    sign silently swaps the east/west trace labels everywhere downstream.
    """
    e_c, n_c, strike = model[0], model[1], model[2]
    s = np.radians(strike)
    d_n, d_e = Y - n_c, X - e_c
    along = d_n * np.cos(s) + d_e * np.sin(s)
    perp = d_n * np.cos(s - np.pi / 2) + d_e * np.sin(s - np.pi / 2)
    return along, perp


def extract_side(gmag: np.ndarray, along: np.ndarray, perp: np.ndarray,
                 sign: float, quality_min: float = QUALITY_MIN
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per along-strike bin: ridge offset, quality, and sample count."""
    s_edges = np.arange(-S_HALF, S_HALF + 1, S_BIN)
    p_edges = np.arange(P_MIN, P_MAX + 1, P_BIN)
    n_s, n_p = len(s_edges) - 1, len(p_edges) - 1
    band = np.isfinite(gmag) & (sign * perp >= P_MIN) & (sign * perp < P_MAX) \
        & (np.abs(along) < S_HALF)
    si = np.floor((along[band] + S_HALF) / S_BIN).astype(int)
    pi = np.floor((sign * perp[band] - P_MIN) / P_BIN).astype(int)
    flat = si * n_p + pi
    sums = np.bincount(flat, weights=gmag[band], minlength=n_s * n_p)
    counts = np.bincount(flat, minlength=n_s * n_p)
    with np.errstate(invalid="ignore"):
        mean = (sums / counts).reshape(n_s, n_p)
    mean[counts.reshape(n_s, n_p) < MIN_SAMPLES] = np.nan

    offsets = np.full(n_s, np.nan)
    quality = np.full(n_s, np.nan)
    p_centres = 0.5 * (p_edges[:-1] + p_edges[1:])
    for i in range(n_s):
        row = mean[i]
        if np.isfinite(row).sum() < 10:
            continue
        k = int(np.nanargmax(row))
        background = np.nanmedian(row)
        if background <= 0:
            continue
        q = row[k] / background
        if q >= quality_min:
            offsets[i] = sign * p_centres[k]
            quality[i] = q
    return offsets, quality, 0.5 * (s_edges[:-1] + s_edges[1:])


def los_step(los: np.ndarray, along: np.ndarray, perp: np.ndarray,
             s_centres: np.ndarray, offsets: np.ndarray,
             sign: float) -> np.ndarray:
    """Per-bin LOS step across the ridge, axis side minus outer side (m).

    The local linear LOS trend (fitted beyond the step bands, out to
    STEP_TREND_M) is removed first, so a smooth gradient contributes ~0 and
    only a localized discontinuity survives. Negative = axis side moved away
    from the satellite (down, for these steep look vectors).
    """
    steps = np.full(len(s_centres), np.nan)
    for i, (s_c, p_r) in enumerate(zip(s_centres, offsets)):
        if not np.isfinite(p_r):
            continue
        u = (perp - p_r) * (-sign)          # positive toward the axis
        cand = (np.abs(along - s_c) < S_BIN / 2) & np.isfinite(los) \
            & (np.abs(u) < STEP_TREND_M)
        if cand.sum() < 60:
            continue
        uu, vv = u[cand], los[cand]
        trend = np.abs(uu) > STEP_BAND_M
        if trend.sum() < 30:
            continue
        M = np.column_stack([np.ones(int(trend.sum())), uu[trend]])
        coef, *_ = np.linalg.lstsq(M, vv[trend], rcond=None)
        detr = vv - (coef[0] + coef[1] * uu)
        near_axis = (uu > STEP_INNER_M) & (uu <= STEP_BAND_M)
        near_out = (uu < -STEP_INNER_M) & (uu >= -STEP_BAND_M)
        if near_axis.sum() < 15 or near_out.sum() < 15:
            continue
        steps[i] = float(np.median(detr[near_axis])
                         - np.median(detr[near_out]))
    return steps


def rolling_median(values: np.ndarray, window: int) -> np.ndarray:
    out = np.full_like(values, np.nan)
    half = window // 2
    for i in range(len(values)):
        seg = values[max(0, i - half): i + half + 1]
        if np.isfinite(seg).sum() >= max(2, window // 2):
            out[i] = np.nanmedian(seg)
    return out


def column_nanmean(values: np.ndarray) -> np.ndarray:
    """Column means with NaN for empty columns and no empty-slice warning."""
    counts = np.isfinite(values).sum(axis=0)
    sums = np.nansum(values, axis=0)
    means = np.full(sums.shape, np.nan, dtype=float)
    np.divide(sums, counts, out=means, where=counts > 0)
    return means


def cross_geometry_alignment_null(first: np.ndarray, second: np.ndarray) -> dict:
    """Circular-shift null for along-strike alignment of two ridge traces."""
    common = np.isfinite(first) & np.isfinite(second)
    if common.sum() < 5:
        return {"observed_median_spread_m": None, "n_null_shifts": 0,
                "fraction_null_as_good_or_better": None}
    observed = float(np.median(np.abs(first[common] - second[common])))
    null = []
    for shift in range(1, len(second)):
        shifted = np.roll(second, shift)
        valid = np.isfinite(first) & np.isfinite(shifted)
        if valid.sum() >= 5:
            null.append(float(np.median(np.abs(first[valid] - shifted[valid]))))
    return {
        "method": "circularly shift the second geometry along strike; preserve "
                  "each trace's values and smoothness",
        "observed_median_spread_m": observed,
        "n_common_bins": int(common.sum()),
        "n_null_shifts": len(null),
        "null_median_spread_m_p05_p50_p95": [
            float(value) for value in np.percentile(null, [5, 50, 95])
        ] if null else None,
        "fraction_null_as_good_or_better": (
            float(np.mean(np.asarray(null) <= observed)) if null else None
        ),
    }


def main() -> None:
    fit_path = RESULTS / "okada_fit.json"
    model = np.asarray(json.loads(fit_path.read_text())["values"][:7])
    from .gslc import DATA_GSLC, UNWRAP_COHERENCE_MIN

    per_geometry = {}
    grids = {}
    for tag in INVERSION_PAIRS:
        f = np.load(DATA_GSLC / f"{tag}_ifg.npz")
        phase = np.where(f["coherence"] > COHERENCE_MIN_GRADIENT,
                         f["phase"], np.nan)
        gmag = wrapped_gradient(phase, 30.0)
        X, Y = np.meshgrid(f["x"], f["y"])
        along, perp = strike_frame(X, Y, model)
        sides = {}
        for name, sign in (("east", +1.0), ("west", -1.0)):
            offsets, quality, s_centres = extract_side(gmag, along, perp, sign)
            sides[name] = {"offset_m": offsets, "quality": quality}
            print(f"{tag} {name}: {int(np.isfinite(offsets).sum())}/"
                  f"{len(offsets)} bins with a ridge, median |offset| "
                  f"{np.nanmedian(np.abs(offsets)):.0f} m")
        per_geometry[tag] = sides
        grids[tag] = (gmag, f["x"], f["y"])

    # Merge the two geometries: mean where both detect, single value where
    # one does; then a rolling median along strike.
    s_centres = np.arange(-S_HALF + S_BIN / 2, S_HALF, S_BIN)
    merged = {}
    for name in ("east", "west"):
        stack = np.vstack([per_geometry[tag][name]["offset_m"]
                           for tag in INVERSION_PAIRS])
        mean = column_nanmean(stack)
        both = np.isfinite(stack).all(axis=0)
        spread = np.abs(stack[0] - stack[1])
        merged[name] = {
            "offset_m": rolling_median(mean, SMOOTH_BINS),
            "n_bins_detected": int(np.isfinite(mean).sum()),
            "n_bins_both_geometries": int(both.sum()),
            "median_cross_geometry_spread_m": float(np.nanmedian(spread[both]))
            if both.any() else None,
        }
        print(f"merged {name}: {merged[name]['n_bins_detected']} bins, "
              f"{merged[name]['n_bins_both_geometries']} in both geometries, "
              f"median cross-geometry spread "
              f"{merged[name]['median_cross_geometry_spread_m']} m")

    threshold_sensitivity = []
    configurations = [
        (f"quality={quality:g}", COHERENCE_MIN_GRADIENT, quality)
        for quality in QUALITY_SENSITIVITY if quality != QUALITY_MIN
    ] + [
        (f"coherence={coherence:g}", coherence, QUALITY_MIN)
        for coherence in COHERENCE_SENSITIVITY
        if coherence != COHERENCE_MIN_GRADIENT
    ]
    for label, coherence_min, quality_min in configurations:
        picked = {side: [] for side in ("east", "west")}
        for tag in INVERSION_PAIRS:
            field = np.load(DATA_GSLC / f"{tag}_ifg.npz")
            phase = np.where(field["coherence"] > coherence_min,
                             field["phase"], np.nan)
            gradient = wrapped_gradient(phase, 30.0)
            X, Y = np.meshgrid(field["x"], field["y"])
            along, perp = strike_frame(X, Y, model)
            for side, sign in (("east", +1.0), ("west", -1.0)):
                offsets, _, _ = extract_side(
                    gradient, along, perp, sign, quality_min=quality_min
                )
                picked[side].append(offsets)
        entry = {"variant": label, "coherence_min": coherence_min,
                 "quality_min": quality_min, "sides": {}}
        for side in ("east", "west"):
            stack = np.vstack(picked[side])
            alternative = rolling_median(column_nanmean(stack), SMOOTH_BINS)
            baseline = np.asarray(merged[side]["offset_m"], dtype=float)
            common = np.isfinite(alternative) & np.isfinite(baseline)
            entry["sides"][side] = {
                "n_merged_bins": int(np.isfinite(alternative).sum()),
                "n_bins_compared_with_baseline": int(common.sum()),
                "median_abs_change_from_baseline_m": (
                    float(np.median(np.abs(alternative[common] - baseline[common])))
                    if common.any() else None
                ),
            }
        threshold_sensitivity.append(entry)

    alignment_null = {
        side: cross_geometry_alignment_null(
            np.asarray(per_geometry[INVERSION_PAIRS[0]][side]["offset_m"]),
            np.asarray(per_geometry[INVERSION_PAIRS[1]][side]["offset_m"]),
        )
        for side in ("east", "west")
    }

    # LOS-step diagnostic across the merged ridges, per geometry, from the
    # unwrapped fields.
    step_stats = {}
    for tag in INVERSION_PAIRS:
        fu = np.load(DATA_GSLC / f"{tag}_unw.npz")
        ulos = np.where(fu["coherence"] > UNWRAP_COHERENCE_MIN,
                        fu["los"], np.nan)
        Xu, Yu = np.meshgrid(fu["x"], fu["y"])
        alo, per = strike_frame(Xu, Yu, model)
        step_stats[tag] = {}
        for name, sign in (("east", +1.0), ("west", -1.0)):
            offs = np.asarray(merged[name]["offset_m"], dtype=float)
            steps = los_step(ulos, alo, per, s_centres, offs, sign)
            okf = np.isfinite(steps)
            step_stats[tag][name] = {
                "n_bins": int(okf.sum()),
                "median_step_cm": (float(100 * np.median(steps[okf]))
                                   if okf.any() else None),
                "fraction_axis_side_down": (float(np.mean(steps[okf] < 0))
                                            if okf.any() else None),
            }
            s = step_stats[tag][name]
            print(f"{tag} {name} LOS step: median "
                  f"{s['median_step_cm']} cm across {s['n_bins']} bins, "
                  f"axis side down in {s['fraction_axis_side_down']}")

    trace_support = {}
    for side in ("east", "west"):
        values = [step_stats[tag][side] for tag in INVERSION_PAIRS]
        consistent = all(value["median_step_cm"] is not None
                         and value["median_step_cm"] < 0
                         and value["fraction_axis_side_down"] >= 0.5
                         for value in values)
        trace_support[side] = {
            "consistent_axis_side_down_step": consistent,
            "interpretation": (
                "ridge and LOS-step evidence are consistent with a scarp"
                if consistent else
                "ridge only; the LOS-step evidence is not consistent across "
                "geometries"
            ),
        }

    def _null_safe(values: np.ndarray) -> list:
        return [float(v) if np.isfinite(v) else None for v in values]

    payload = {
        "method": "wrapped-phase-gradient ridge per 500-m along-strike bin, "
                  "0.3-4 km each side of the fitted dike trace, detected per "
                  "NISAR geometry and merged. The picker follows the "
                  "phase-gradient idea of Xu et al. 2026 "
                  "(doi:10.1016/j.geog.2026.06.002) but automates a narrower "
                  "task (per-bin ridge localization with a peak/median "
                  "threshold); it selects at most one ridge per side and "
                  "does not by itself establish a fault break — see "
                  "los_step_across_ridges",
        "strike_deg": float(model[2]),
        "trace_centre_utm": [float(model[0]), float(model[1])],
        "s_centres_m": s_centres.tolist(),
        "east": {k: (_null_safe(v) if isinstance(v, np.ndarray) else v)
                 for k, v in merged["east"].items()},
        "west": {k: (_null_safe(v) if isinstance(v, np.ndarray) else v)
                 for k, v in merged["west"].items()},
        "per_geometry_quality_median": {
            tag: {name: float(np.nanmedian(per_geometry[tag][name]["quality"]))
                  for name in ("east", "west")}
            for tag in INVERSION_PAIRS},
        "los_step_across_ridges": {
            "method": "median LOS difference between 0.1-0.4 km bands on the "
                      "axis and outer sides of each ridge bin, local linear "
                      "trend (0.4-1.2 km) removed; negative = axis side "
                      "moved away from the satellite (down)",
            "per_geometry": step_stats,
        },
        "trace_support": trace_support,
        "threshold_sensitivity": threshold_sensitivity,
        "cross_geometry_alignment_null": alignment_null,
        "parameters": {"s_bin_m": S_BIN, "p_band_m": [P_MIN, P_MAX],
                       "quality_min": QUALITY_MIN,
                       "coherence_min": COHERENCE_MIN_GRADIENT,
                       "quality_sensitivity": list(QUALITY_SENSITIVITY),
                       "coherence_sensitivity": list(COHERENCE_SENSITIVITY),
                       "step_bands_m": [STEP_INNER_M, STEP_BAND_M,
                                        STEP_TREND_M]},
    }
    (RESULTS / "fault_traces.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n")

    # Figure: gradient maps with merged traces, plus one cross-axis profile.
    e_c, n_c, strike = model[0], model[1], model[2]
    sr = np.radians(strike)
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.5))
    for ax, tag in zip(axes[:2], INVERSION_PAIRS):
        gmag, x, y = grids[tag]
        im = ax.imshow(gmag * 1e3, cmap="gray_r", vmin=0, vmax=8,
                       extent=[x[0] / 1e3, x[-1] / 1e3,
                               y[-1] / 1e3, y[0] / 1e3])
        for name, colour in (("east", "tab:red"), ("west", "tab:blue")):
            p = np.asarray(merged[name]["offset_m"], dtype=float)
            ok = np.isfinite(p)
            # Inverse of strike_frame: along = (cos s, sin s) and
            # perp = (sin s, -cos s) in (north, east). The opposite perp sign
            # here once plotted the east ridge on the west side.
            e_line = (e_c + s_centres * np.sin(sr) - p * np.cos(sr))[ok]
            n_line = (n_c + s_centres * np.cos(sr) + p * np.sin(sr))[ok]
            ax.plot(e_line / 1e3, n_line / 1e3, ".", color=colour, ms=2,
                    label=f"{name} ridge")
        ax.set_xlim(660, 705)
        ax.set_ylim(1330, 1370)
        ax.set_title(f"{tag}: wrapped-phase gradient + merged ridges",
                     fontsize=10)
        ax.set_xlabel("UTM 37N easting (km)")
        ax.set_ylabel("UTM 37N northing (km)")
        ax.legend(fontsize=8, loc="lower left")
        plt.colorbar(im, ax=ax, shrink=0.75, label="|∇φ| (mrad/m)")

    # Cross-axis profile of gradient at the centre bins, ascending geometry.
    gmag, x, y = grids[INVERSION_PAIRS[0]]
    X, Y = np.meshgrid(x, y)
    along, perp = strike_frame(X, Y, model)
    sel = np.abs(along) < 2_000
    p_edges = np.arange(-4_000, 4_001, 100.0)
    pi = np.digitize(perp[sel], p_edges)
    prof = np.array([np.nanmedian(gmag[sel][pi == k])
                     if (pi == k).sum() > 30 else np.nan
                     for k in range(1, len(p_edges))])
    centres = 0.5 * (p_edges[:-1] + p_edges[1:])
    axes[2].plot(centres / 1e3, prof * 1e3, "k-")
    for name, colour in (("east", "tab:red"), ("west", "tab:blue")):
        p = np.asarray(merged[name]["offset_m"], dtype=float)
        mid = np.nanmedian(p[np.abs(s_centres) < 2_000])
        if np.isfinite(mid):
            axes[2].axvline(mid / 1e3, color=colour, ls="--",
                            label=f"{name} ridge ({mid:+.0f} m)")
    axes[2].axvline(0, color="0.5", ls=":", label="dike trace")
    axes[2].set_xlabel("distance from dike trace, ENE positive (km)")
    axes[2].set_ylabel("median |∇φ| (mrad/m)")
    axes[2].set_title("cross-axis gradient profile (|along| < 2 km, asc)",
                      fontsize=10)
    axes[2].legend(fontsize=8)
    fig.suptitle("Phase-gradient ridges flanking the rift axis "
                 "(30-m reprocessed interferograms)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = RESULTS / "fault_traces.png"
    fig.savefig(out, dpi=110)
    print(f"figure: {out}")


if __name__ == "__main__":
    main()
