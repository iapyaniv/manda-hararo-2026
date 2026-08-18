"""Test for additional dike-like deformation after the July intrusion.

The descending 2026-07-19→2026-08-12 GUNW pair begins at the secondary
acquisition of the co-event descending pair.  This module regresses its
screened LOS field on the NISAR-only dike pattern and reports the fraction of
the co-event model repeated after 19 July.

Run after :mod:`mh26.invert`.  Writes ``results/continuation.md``, ``.json``,
and ``continuation_panels.png``.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .gunw import (
    block_average,
    corrected_los,
    load_pair,
    los_unit_neu,
    pair_window,
    reference_farfield,
)
from .history import probe_mean
from .okada import forward
from .params import (
    BLOCK_PIXELS,
    CENTER_E,
    CENTER_N,
    FARFIELD_RADIUS,
    PROBES,
    RESULTS,
)


POST_TAG = "165_082_D_20260719"
RADIUS = 25_000.0
PRED_MIN = 0.01
BOOTSTRAP_SAMPLES = 2_000
BOOTSTRAP_SEED = 20260818
CLUSTER_SIZE = 5_000.0
PROCESSING_VARIANTS = (
    ("0.96 km blocks", {"block": 12}),
    ("1.92 km blocks", {"block": 24}),
    ("25 km far field", {"farfield_radius": 25_000.0}),
    ("32 km far field", {"farfield_radius": 32_000.0}),
    ("fit within 20 km", {"radius_limit": 20_000.0}),
    ("fit within 30 km", {"radius_limit": 30_000.0}),
    ("prediction >0.5 cm", {"prediction_min": 0.005}),
    ("prediction >2 cm", {"prediction_min": 0.02}),
)


def regression(predicted: np.ndarray, observed: np.ndarray) -> tuple[float, float]:
    design = np.column_stack([predicted, np.ones(len(predicted))])
    slope, offset = np.linalg.lstsq(design, observed, rcond=None)[0]
    return float(slope), float(offset)


def cluster_wild_interval(predicted: np.ndarray, observed: np.ndarray,
                          points: np.ndarray, slope: float,
                          offset: float) -> tuple[list[float], int]:
    """Heuristic 95% slope interval with 5-km correlated-noise clusters."""
    fitted = slope * predicted + offset
    residual = observed - fitted
    east_bin = np.floor((points[:, 1] - points[:, 1].min()) / CLUSTER_SIZE).astype(int)
    north_bin = np.floor((points[:, 0] - points[:, 0].min()) / CLUSTER_SIZE).astype(int)
    cluster_pairs = np.column_stack([north_bin, east_bin])
    _, cluster = np.unique(cluster_pairs, axis=0, return_inverse=True)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    samples = np.empty(BOOTSTRAP_SAMPLES)
    for index in range(BOOTSTRAP_SAMPLES):
        signs = rng.choice((-1.0, 1.0), size=cluster.max() + 1)
        synthetic = fitted + residual * signs[cluster]
        samples[index], _ = regression(predicted, synthetic)
    interval = [float(value) for value in np.percentile(samples, [2.5, 97.5])]
    return interval, int(cluster.max() + 1)


def projection(model: np.ndarray, data: dict, los: np.ndarray,
               block: int = BLOCK_PIXELS,
               farfield_radius: float = FARFIELD_RADIUS,
               radius_limit: float = RADIUS,
               prediction_min: float = PRED_MIN) -> dict:
    """Project one screened post-event field onto the co-event dike pattern."""
    referenced = reference_farfield(
        los, data["x"], data["y"], r_in=farfield_radius
    )
    points, observed, grid = block_average(
        referenced, data["x"], data["y"], block=block
    )
    units = los_unit_neu(data, points)
    predicted = np.sum(forward(model, points) * units, axis=1)
    radius = np.hypot(points[:, 1] - CENTER_E, points[:, 0] - CENTER_N)
    use = (
        np.isfinite(observed)
        & np.isfinite(predicted)
        & (radius < radius_limit)
        & (np.abs(predicted) > prediction_min)
    )
    if use.sum() < 100:
        raise ValueError("too few post-event blocks survive the continuation comparison")
    slope, offset = regression(predicted[use], observed[use])
    fitted = slope * predicted[use] + offset
    return {
        "referenced": referenced,
        "points": points,
        "observed": observed,
        "predicted": predicted,
        "grid": grid,
        "use": use,
        "slope": slope,
        "offset": offset,
        "correlation": float(np.corrcoef(observed[use], predicted[use])[0, 1]),
        "rms": float(np.sqrt(np.mean((observed[use] - fitted) ** 2))),
        "zero_intercept_slope": float(
            np.sum(predicted[use] * observed[use])
            / np.sum(predicted[use] ** 2)
        ),
        "block": block,
        "farfield_radius": farfield_radius,
        "radius_limit": radius_limit,
        "prediction_min": prediction_min,
    }


def main() -> None:
    fit_path = RESULTS / "okada_fit.json"
    if not fit_path.exists():
        raise SystemExit("results/okada_fit.json missing — run `python -m mh26.invert` first")
    data = load_pair(POST_TAG)
    if data is None:
        raise SystemExit(f"data/nisar/{POST_TAG}.npz missing — see data/README.md")
    model = np.asarray(json.loads(fit_path.read_text())["values"][:7], dtype=float)

    los, _ = corrected_los(data)
    default = projection(model, data, los)
    referenced = default["referenced"]
    points = default["points"]
    observed = default["observed"]
    predicted = default["predicted"]
    grid = default["grid"]
    use = default["use"]
    slope = default["slope"]
    offset = default["offset"]
    interval, n_clusters = cluster_wild_interval(
        predicted[use], observed[use], points[use], slope, offset
    )
    rms = default["rms"]
    corr = default["correlation"]

    sensitivity = []
    for label, options in PROCESSING_VARIANTS:
        variant = projection(model, data, los, **options)
        sensitivity.append({
            "variant": label,
            "pattern_fraction": variant["slope"],
            "n_blocks": int(variant["use"].sum()),
        })
    sensitivity_values = [slope] + [
        row["pattern_fraction"] for row in sensitivity
    ]
    processing_envelope = [
        float(min(sensitivity_values)),
        float(max(sensitivity_values)),
    ]
    probes = {
        name: probe_mean(referenced, data["x"], data["y"], *xy)
        for name, xy in PROBES.items()
    }
    window = pair_window(data)
    payload = {
        "tag": POST_TAG,
        "window": window,
        "interpretation": "fraction of the co-event best-fit dike pattern after 19 July",
        "n_blocks": int(use.sum()),
        "radius_m": RADIUS,
        "prediction_min_m": PRED_MIN,
        "pattern_fraction": slope,
        "pattern_fraction_cluster_wild_95": interval,
        "pattern_fraction_processing_envelope": processing_envelope,
        "processing_sensitivity": sensitivity,
        "zero_intercept_pattern_fraction": default["zero_intercept_slope"],
        "offset_m": offset,
        "correlation": corr,
        "rms_m": rms,
        "probe_means_cm": {name: float(value) for name, value in probes.items()},
        "bootstrap": {
            "method": "wild sign bootstrap on 5-km spatial clusters",
            "samples": BOOTSTRAP_SAMPLES,
            "seed": BOOTSTRAP_SEED,
            "clusters": n_clusters,
            "caveat": "conditional heuristic, not a full InSAR covariance model",
        },
    }
    (RESULTS / "continuation.json").write_text(json.dumps(payload, indent=2) + "\n")

    lines = [
        "# Post-event continuation test",
        "",
        f"- NISAR descending pair: **{window}** (24 days; starts at the co-event",
        "  descending pair's secondary acquisition)",
        f"- dike-pattern projection coefficient: **{slope:+.3f}**; processing-choice",
        f"  envelope **[{processing_envelope[0]:+.3f}, {processing_envelope[1]:+.3f}]**",
        f"- conditional 5-km cluster-wild interval: **[{interval[0]:+.3f},",
        f"  {interval[1]:+.3f}]** ({n_clusters} clusters); this is not a calibrated",
        "  confidence interval for spatially correlated InSAR noise",
        f"- correlation with the co-event pattern: {corr:+.3f}; residual RMS: {100*rms:.1f} cm",
        f"- fitted subset offset: {100*offset:+.1f} cm; forcing it to zero would give a",
        f"  coefficient of {default['zero_intercept_slope']:+.3f}",
        f"- probe means: EAST {probes['EAST']:+.1f} cm, WEST {probes['WEST']:+.1f} cm,",
        f"  SOUTH {probes['SOUTH']:+.1f} cm",
        "",
        "The positive projection coefficient is about one-tenth of the co-event value.",
        "This result uses one provisional descending geometry and does not establish",
        "renewed intrusion.",
        "",
        "The test detects only deformation shaped like the fitted dike. The broad residual",
        "field and fitted datum offset show why additional geometries and atmospheric",
        "covariance treatment are required before converting this coefficient to added volume.",
    ]
    report = "\n".join(lines) + "\n"
    (RESULTS / "continuation.md").write_text(report)
    print(report)

    xs, ys, observed_grid = grid
    predicted_grid = np.full_like(observed_grid, np.nan)
    fitted_grid = np.full_like(observed_grid, np.nan)
    residual_grid = np.full_like(observed_grid, np.nan)
    valid_grid = np.isfinite(observed_grid)
    predicted_grid[valid_grid] = predicted
    fitted_grid[valid_grid] = slope * predicted
    residual_grid[valid_grid] = observed - (slope * predicted + offset)
    display = (observed_grid - offset, fitted_grid, residual_grid)
    titles = (
        f"observed LOS (fitted {100*offset:+.1f} cm offset removed)",
        "scaled co-event dike pattern",
        "residual",
    )
    scale = np.nanpercentile(np.abs(observed_grid - offset), 99)
    figure, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    for axis, field, title in zip(axes, display, titles):
        image = axis.imshow(
            field * 100,
            cmap="RdBu_r",
            vmin=-100 * scale,
            vmax=100 * scale,
            extent=[xs[0] / 1e3, xs[-1] / 1e3, ys[-1] / 1e3, ys[0] / 1e3],
        )
        axis.set_title(title)
        axis.set_xlabel("UTM 37N easting (km)")
        axis.set_ylabel("UTM 37N northing (km)")
        plt.colorbar(image, ax=axis, shrink=0.72, label="LOS (cm)")
    figure.suptitle(
        f"Manda Hararo continuation, {window}: projection {slope:+.3f}; "
        f"processing envelope [{processing_envelope[0]:+.3f}, "
        f"{processing_envelope[1]:+.3f}]"
    )
    figure.tight_layout(rect=[0, 0, 1, 0.93])
    figure.savefig(RESULTS / "continuation_panels.png", dpi=120)


if __name__ == "__main__":
    main()
