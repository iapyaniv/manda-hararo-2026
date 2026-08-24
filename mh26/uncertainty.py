"""correlated-noise Monte Carlo for the distributed model.

The sensitivity envelopes bracket processing and geometry choices; they say
nothing about data noise. This module propagates spatially correlated noise
through the fixed-geometry inversion:

1. per dataset, a nugget + exponential semivariogram
   gamma(h) = n + c (1 - exp(-h/L)) is fitted to the empirical far-field
   semivariogram of the block observations (L searched to 30 km; the fit
   parameters and binned semivariogram are stored with the result). The
   observations are already far-field-plane referenced. The fit estimates a
   stationary kernel from those values but does not reproduce the reference
   constraints by itself;
2. the best-model prediction and every stationary covariance realization are
   passed through the block-space far-field plane projection. The synthetic
   fields therefore have the same datum and ramp constraints as the block
   observations;
3. each realization is re-inverted with the same geometry, bounds, and
   smoothing weight.

The noiseless refit differs from the committed best fit for two combined
reasons: the smoothing penalty acts again on a smooth input, and the plane
referencing removes part of the model's own far-field signal. The stored
noiseless refit therefore measures the COMBINED estimator bias (smoothing +
referencing), and raw percentiles of the refits are NOT intervals around
the best fit. The reported 95% intervals are pivotal (basic bootstrap): the
refit errors about the noiseless refit, re-centred on the best-fit value.
They describe conditional data-noise dispersion under this simulator, with
geometry and lambda fixed; they are not a full physical-parameter
confidence interval. Positivity truncates the noise response.

Target fit: results/distributed_fit.json (idealized fault geometry, GUNW
inputs) — the primary reported model. The mapped-ridge and GSLC variants are
bracketed separately by the sensitivity envelopes.

Run after `mh26.distributed`:  python -m mh26.uncertainty [--n 1000]

Writes results/distributed_uncertainty.{json,md}.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np

from .distributed import (
    assemble,
    dataset_rms,
    geometry_from_json,
    geometry_sha256,
    load_datasets,
    moment_summary,
    solve,
)
from .params import CENTER_E, CENTER_N, FARFIELD_RADIUS, RESULTS

RNG_SEED = 20260819
VARIOGRAM_MAX_BLOCKS = 1500
VARIOGRAM_BIN_M = 2_000.0
VARIOGRAM_MAX_LAG_M = 40_000.0
RANGE_GRID_M = np.arange(500.0, 30_001.0, 250.0)
QUANTITY_KEYS = ("dike_volume_km3", "max_opening_m", "fault_east_max_slip_m",
                 "fault_west_max_slip_m", "fault_to_dike_potency_ratio",
                 "opening_centroid_depth_km")


def fit_covariance(pts: np.ndarray, vals: np.ndarray,
                   rng: np.random.Generator) -> dict:
    """Nugget + exponential covariance from the far-field semivariogram."""
    radius = np.hypot(pts[:, 1] - CENTER_E, pts[:, 0] - CENTER_N)
    far = radius > FARFIELD_RADIUS
    p = pts[far]
    v = vals[far] - np.mean(vals[far])
    if len(v) > VARIOGRAM_MAX_BLOCKS:
        pick = rng.choice(len(v), VARIOGRAM_MAX_BLOCKS, replace=False)
        p, v = p[pick], v[pick]
    d = np.hypot(p[:, 0, None] - p[None, :, 0], p[:, 1, None] - p[None, :, 1])
    g = 0.5 * (v[:, None] - v[None, :]) ** 2
    iu = np.triu_indices(len(v), k=1)
    d, g = d[iu], g[iu]
    use = d < VARIOGRAM_MAX_LAG_M
    bins = np.floor(d[use] / VARIOGRAM_BIN_M).astype(int)
    n_bins = int(VARIOGRAM_MAX_LAG_M / VARIOGRAM_BIN_M)
    counts = np.bincount(bins, minlength=n_bins)
    gamma = np.bincount(bins, weights=g[use], minlength=n_bins) / np.maximum(counts, 1)
    centres = (np.arange(n_bins) + 0.5) * VARIOGRAM_BIN_M
    ok = counts > 200
    best = None
    for length in RANGE_GRID_M:
        basis = 1.0 - np.exp(-centres[ok] / length)
        # gamma ~ nugget + c * basis, both non-negative.
        M = np.column_stack([np.ones(ok.sum()), basis])
        coef, *_ = np.linalg.lstsq(M, gamma[ok], rcond=None)
        nugget, c = max(coef[0], 0.0), max(coef[1], 0.0)
        if coef[0] < 0:
            c = max(float(np.sum(basis * gamma[ok]) / np.sum(basis ** 2)), 0.0)
        elif coef[1] < 0:
            nugget, c = float(np.mean(gamma[ok])), 0.0
        err = float(np.mean((nugget + c * basis - gamma[ok]) ** 2))
        if best is None or err < best["sse"]:
            best = {"sse": err, "nugget_m2": nugget, "partial_sill_m2": c,
                    "range_m": float(length)}
    best["range_at_search_limit"] = bool(best["range_m"] >= RANGE_GRID_M[-1] - 1)
    best["n_far_blocks"] = int(far.sum())
    best["semivariogram"] = {
        "lag_centres_m": centres[ok].tolist(),
        "gamma_m2": gamma[ok].tolist(),
        "pair_counts": counts[ok].tolist(),
    }
    return best


def remove_farfield_plane(pts: np.ndarray, vals: np.ndarray) -> np.ndarray:
    """Apply the block-space form of the far-field datum/ramp projection."""
    radius = np.hypot(pts[:, 1] - CENTER_E, pts[:, 0] - CENTER_N)
    far = radius > FARFIELD_RADIUS
    M = np.column_stack([np.ones(far.sum()), pts[far, 1] - CENTER_E,
                         pts[far, 0] - CENTER_N])
    coef, *_ = np.linalg.lstsq(M, vals[far], rcond=None)
    full = np.column_stack([np.ones(len(vals)), pts[:, 1] - CENTER_E,
                            pts[:, 0] - CENTER_N])
    return vals - full @ coef


def referenced_noise_draw(pts: np.ndarray, factor: np.ndarray,
                          rng: np.random.Generator) -> np.ndarray:
    """Draw stationary noise and impose the far-field plane constraints."""
    raw = factor @ rng.standard_normal(factor.shape[0])
    return remove_farfield_plane(pts, raw)


def quantities_from(summary: dict) -> dict:
    potency = (summary.get("fault_east_potency_m3", 0.0)
               + summary.get("fault_west_potency_m3", 0.0))
    return {
        "dike_volume_km3": summary["dike_volume_km3"],
        "max_opening_m": summary["max_opening_m"],
        "fault_east_max_slip_m": summary.get("fault_east_max_slip_m", 0.0),
        "fault_west_max_slip_m": summary.get("fault_west_max_slip_m", 0.0),
        "fault_to_dike_potency_ratio":
            potency / (summary["dike_volume_km3"] * 1e9),
        "opening_centroid_depth_km": summary["opening_centroid_depth_km"],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=1000)
    args = ap.parse_args()

    target_path = RESULTS / "distributed_fit.json"
    target_bytes = target_path.read_bytes()
    base = json.loads(target_bytes)
    target_sha256 = hashlib.sha256(target_bytes).hexdigest()
    lam = float(base["lambda_chosen"])
    model = np.asarray(json.loads((RESULTS / "okada_fit.json").read_text())["values"][:7])
    decay = json.loads((RESULTS / "s1_tracks.json").read_text())["post_event_decay"]

    geometry = geometry_from_json(base["geometry"])
    if base.get("geometry_sha256") != geometry_sha256(geometry):
        raise ValueError("distributed_fit.json geometry hash does not match")
    datasets = load_datasets(model, decay, nisar_source="gunw")
    n_sets = len(datasets)
    A, b, keep = assemble(datasets, geometry, with_faults=True)
    x_best = solve(A, b, keep, geometry, lam, n_sets)
    dataset_rms(datasets, geometry, keep, x_best)  # caches predictions
    q_best = quantities_from(moment_summary(geometry, keep, x_best))

    rng = np.random.default_rng(RNG_SEED)
    factors, covariances = [], {}
    for ds in datasets:
        cov_fit = fit_covariance(ds["pts"], ds["vals"], rng)
        covariances[ds["name"]] = cov_fit
        d = np.hypot(ds["pts"][:, 0, None] - ds["pts"][None, :, 0],
                     ds["pts"][:, 1, None] - ds["pts"][None, :, 1])
        cov = cov_fit["partial_sill_m2"] * np.exp(-d / cov_fit["range_m"])
        cov[np.diag_indices_from(cov)] += cov_fit["nugget_m2"] + 1e-10
        factors.append(np.linalg.cholesky(cov))
        print(f"{ds['name']}: sigma "
              f"{100 * np.sqrt(cov_fit['nugget_m2'] + cov_fit['partial_sill_m2']):.1f} cm "
              f"(nugget {100 * np.sqrt(cov_fit['nugget_m2']):.1f} cm), range "
              f"{cov_fit['range_m'] / 1e3:.1f} km"
              + (" [at search limit]" if cov_fit["range_at_search_limit"] else ""))

    predictions = [ds["_greens"][:, keep] @ x_best[:keep.sum()]
                   + x_best[keep.sum() + i]
                   for i, ds in enumerate(datasets)]
    # Apply the deterministic reference projection to the model. Each random
    # draw receives this projection separately below.
    predictions_ref = [remove_farfield_plane(ds["pts"], pred)
                       for ds, pred in zip(datasets, predictions)]

    def solve_simulated(noises: list[np.ndarray] | None) -> dict:
        rows = []
        for i, (ds, pred) in enumerate(zip(datasets, predictions_ref)):
            sim = pred if noises is None else pred + noises[i]
            rows.append(sim / ds["sigma"])
        x = solve(A, np.concatenate(rows), keep, geometry, lam, n_sets)
        return quantities_from(moment_summary(geometry, keep, x))

    # Noiseless refit: the combined estimator bias (smoothing + referencing).
    q_ref = solve_simulated(None)
    print("noiseless refit: " + ", ".join(f"{k}={v:.4g}"
                                          for k, v in q_ref.items()))

    draws = {k: [] for k in QUANTITY_KEYS}

    def write_payload(n_done: int) -> None:
        percentiles = [2.5, 16, 50, 84, 97.5]
        raw = {k: {f"p{p:g}": float(np.percentile(v, p)) for p in percentiles}
               for k, v in draws.items()}
        pivotal = {}
        for k, v in draws.items():
            err = np.asarray(v) - q_ref[k]
            pivotal[k] = {
                "lo95": float(q_best[k] - np.percentile(err, 97.5)),
                "hi95": float(q_best[k] - np.percentile(err, 2.5)),
                "lo68": float(q_best[k] - np.percentile(err, 84)),
                "hi68": float(q_best[k] - np.percentile(err, 16)),
            }
        payload = {
            "n_realizations": n_done,
            "target_fit": "distributed_fit.json (idealized fault geometry, "
                          "GUNW inputs)",
            "target_fit_sha256": target_sha256,
            "geometry_sha256": base["geometry_sha256"],
            "lambda": lam,
            "seed": RNG_SEED,
            "noise_model": "per-dataset nugget + exponential stationary "
                           "kernel fitted to the far-field semivariogram of "
                           "the referenced observations; every stationary "
                           "draw is passed through the block-space far-field "
                           "datum/ramp projection before it is added to the "
                           "projected best-model prediction; fixed saved "
                           "geometry and lambda",
            "reference_projection": {
                "basis": ["constant", "easting", "northing"],
                "fit_region": f"radius > {FARFIELD_RADIUS / 1e3:g} km",
                "applied_to": ["best-model prediction", "every noise draw"],
                "space": "block observations",
            },
            "covariances": covariances,
            "best_fit": q_best,
            "noiseless_refit": q_ref,
            "noiseless_refit_meaning": "combined estimator bias: the "
                                       "smoothing penalty acting again plus "
                                       "the far-field referencing removing "
                                       "part of the model's own signal",
            "interval_method": "pivotal (basic bootstrap): refit errors "
                               "about the noiseless refit, re-centred on the "
                               "best fit; conditional data-noise dispersion "
                               "under this simulator, not a full "
                               "physical-parameter confidence interval; raw "
                               "refit percentiles are the estimator "
                               "dispersion, not an interval around the best "
                               "fit",
            "intervals_95": pivotal,
            "refit_percentiles_raw": raw,
            "draws": {k: [float(x) for x in v] for k, v in draws.items()},
        }
        (RESULTS / "distributed_uncertainty.json").write_text(
            json.dumps(payload, indent=2, allow_nan=False) + "\n")
        return payload

    for it in range(args.n):
        noises = [referenced_noise_draw(ds["pts"], factor, rng)
                  for ds, factor in zip(datasets, factors)]
        q = solve_simulated(noises)
        for k in QUANTITY_KEYS:
            draws[k].append(q[k])
        if (it + 1) % 50 == 0:
            print(f"  {it + 1}/{args.n} realizations", flush=True)
        if (it + 1) % 100 == 0:
            write_payload(it + 1)

    payload = write_payload(args.n)

    lines = [
        "# Correlated-noise Monte Carlo (distributed model)", "",
        f"- target: the idealized-geometry GUNW fit (distributed_fit.json); "
        f"{args.n} realizations at fixed geometry and λ = {lam:.2f}",
        "- noise: per-dataset nugget + exponential stationary kernel fitted "
        "to the far-field semivariogram of the referenced observations; the "
        "best-model prediction and every covariance draw receive the "
        "block-space far-field datum/ramp projection",
        "- intervals: pivotal (basic bootstrap) about the noiseless refit, "
        "re-centred on the best fit. The noiseless refit differs from the "
        "best fit through the combined estimator bias — the smoothing "
        "penalty acting again plus the referencing removing part of the "
        "model's own far-field signal — so raw refit percentiles are "
        "reported separately as estimator dispersion. The intervals are "
        "conditional data-noise dispersion under this simulator, not a "
        "full physical-parameter confidence interval",
        "",
        "| quantity | best fit | noiseless refit | 95% interval |",
        "|---|---|---|---|",
    ]
    for k in QUANTITY_KEYS:
        iv = payload["intervals_95"][k]
        lines.append(f"| {k} | {payload['best_fit'][k]:.3f} | "
                     f"{payload['noiseless_refit'][k]:.3f} | "
                     f"{iv['lo95']:.3f}–{iv['hi95']:.3f} |")
    lines.extend([
        "",
        "| quantity | raw refit 2.5% | 50% | 97.5% |",
        "|---|---|---|---|",
    ])
    for k in QUANTITY_KEYS:
        row = payload["refit_percentiles_raw"][k]
        lines.append(f"| {k} | {row['p2.5']:.3f} | {row['p50']:.3f} | "
                     f"{row['p97.5']:.3f} |")
    limits = [name for name, c in payload["covariances"].items()
              if c["range_at_search_limit"]]
    lines.extend([
        "",
        "The intervals describe data noise only, for fixed geometry and",
        "regularization; positivity truncates the noise response near zero",
        "slip, and the pivotal re-centring assumes the refit error about the",
        "noiseless refit approximates the estimator error about the truth.",
        "Geometry and processing choices are bracketed separately by the",
        "sensitivity envelopes in the distributed-fit reports."
        + ("" if not limits else
           f" Range search hit its 30 km limit for: {', '.join(limits)}."),
    ])
    (RESULTS / "distributed_uncertainty.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
