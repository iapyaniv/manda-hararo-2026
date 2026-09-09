"""Depth-extent sensitivity of the trace-anchored graben-fault term.

The distributed fits leave the east graben fault free from the surface down
to where each column meets the dike plane (mh26.distributed). Smoothing and
positivity spread slip down that surface, so the depth reach of the fitted
slip — and with it the fault potency Σ area·slip and the scalar moment
μ · Σ area·slip — is partly regularization, not data. This module refits the
mapped-trace model with the fault truncated at a sweep of depth caps and
records what each truncation costs in misfit and removes in potency, so a
statement about the fault's moment can cite the data-required minimum
instead of the regularized value.

Two facts frame the numbers:

* Okada surface displacements depend on the Lamé constants only through
  alpha = (λ+μ)/(λ+2μ), so uniformly rescaling (λ, μ) changes no prediction
  and rescales every quoted moment linearly; the data do not constrain the
  assumed μ = 30 GPa. A depth-varying μ is outside the half-space model.
* The fault term acts almost entirely inside the graben band (within 6 km of
  the trace, 20 km along strike); the all-block RMS dilutes it, so misfit is
  scored on both.

Run after the mapped distributed fits (mh26.distributed --mapped [--gslc]):

    python -m mh26.fault_depth           # GUNW-chain variant
    python -m mh26.fault_depth --gslc    # 30-m reprocessed chain

Writes results/fault_depth_mapped[_gslc].{json,md,png}.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np

from .distributed import (
    assemble,
    build_geometry,
    geometry_sha256,
    load_datasets,
    load_fault_trace_source,
    moment_summary,
    solve,
)
from .gslc import NEAR_AXIS_ALONG_M, _axis_distance
from .params import RESULTS, SHEAR_MODULUS

import matplotlib.pyplot as plt

# Truncation depths for the fault columns, metres. The uncapped geometry ends
# each column at the dike plane or 4 km (mh26.distributed), so the deepest cap
# below equals a pure re-solve check only where trace offsets allow 4 km.
DEPTH_CAPS_M = (500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0)

# The graben band doubles the near-axis half-width (3 km) because the mapped
# trace columns sit up to 3.3 km off the axis and their hanging-wall lobes
# extend beyond the trace; the along-strike window matches the near-axis band.
GRABEN_PERP_M = 6_000.0
GRABEN_ALONG_M = NEAR_AXIS_ALONG_M

IMPROVEMENT_TARGET = 0.9


def band_masks(datasets: list[dict], model: np.ndarray) -> dict:
    masks = {}
    for ds in datasets:
        perp, along = _axis_distance(ds["pts"][:, 1], ds["pts"][:, 0], model)
        masks[ds["tag"]] = (perp < GRABEN_PERP_M) & (along < GRABEN_ALONG_M)
    return masks


def rms_all_and_band(datasets: list[dict], keep: np.ndarray, x: np.ndarray,
                     masks: dict) -> tuple[float, float]:
    """Physical RMS over all blocks and over the graben band, metres.

    Uses the Green's functions cached by the immediately preceding
    assemble() call, so it must run before the next geometry is assembled.
    """
    n_par = keep.sum()
    res, res_band = [], []
    for i, ds in enumerate(datasets):
        pred = ds["_greens"][:, keep] @ x[:n_par] + x[n_par + i]
        r = pred - ds["vals"]
        res.append(r)
        res_band.append(r[masks[ds["tag"]]])
    cat = np.concatenate(res)
    cat_band = np.concatenate(res_band)
    return (float(np.sqrt(np.mean(cat ** 2))),
            float(np.sqrt(np.mean(cat_band ** 2))))


def fault_depth_profile(geometry: dict, x: np.ndarray,
                        side: str = "east") -> tuple[np.ndarray, np.ndarray]:
    """Cumulative fault potency above depth z, piecewise linear in z.

    Potency density is uniform within a patch (slip is constant on it), so
    the cumulative curve is exact for the fitted solution.
    """
    kinds = geometry["kinds"]
    sel = kinds == f"fault_{side}"
    patches = geometry["patches"][sel]
    areas = geometry["areas"][sel]
    slip = x[:len(kinds)][sel]
    dip_r = np.radians(geometry["fault_dip"])
    z_top = patches[:, 2] + patches[:, 7] * np.sin(dip_r)
    z_bot = patches[:, 2] + patches[:, 8] * np.sin(dip_r)
    grid = np.linspace(0.0, float(z_bot.max()), 81)
    cum = np.zeros_like(grid)
    for a, s, za, zb in zip(areas, slip, z_top, z_bot):
        cum += a * s * np.clip((grid - za) / (zb - za), 0.0, 1.0)
    return grid, cum


def run(nisar_source: str = "gunw") -> None:
    suffix = "" if nisar_source == "gunw" else "_gslc"
    ref_path = RESULTS / f"distributed_fit_mapped{suffix}.json"
    if not ref_path.exists():
        raise SystemExit(f"{ref_path.name} missing — run the mapped "
                         "distributed fit first")
    ref_raw = ref_path.read_bytes()
    ref = json.loads(ref_raw)
    lam = float(ref["lambda_chosen"])
    model = np.asarray(
        json.loads((RESULTS / "okada_fit.json").read_text())["values"][:7])
    decay = json.loads(
        (RESULTS / "s1_tracks.json").read_text())["post_event_decay"]

    trace_offsets, trace_prov = load_fault_trace_source(None)
    support = trace_prov.get("trace_support") or {}
    sides = tuple(side for side in ("east", "west")
                  if support.get(side, {}).get("consistent_axis_side_down_step"))
    if list(sides) != ref["geometry"]["enabled_fault_sides"]:
        raise SystemExit("trace support rule and the published mapped fit "
                         "disagree on modeled fault sides — rerun the mapped "
                         "fit before the sweep")

    datasets = load_datasets(model, decay, nisar_source=nisar_source)
    masks = band_masks(datasets, model)
    n_band = int(sum(int(m.sum()) for m in masks.values()))

    # Baseline: the published mapped geometry, re-solved, and required to
    # reproduce the stored fit so the sweep cannot drift from the published
    # numbers.
    geometry = build_geometry(model, fault_trace_offsets=trace_offsets,
                              fault_trace_source="mapped",
                              enabled_fault_sides=sides)
    if geometry_sha256(geometry) != ref["geometry_sha256"]:
        raise SystemExit(f"rebuilt mapped geometry differs from "
                         f"{ref_path.name} — rerun the mapped fit, then the "
                         "sweep")
    A, b, keep = assemble(datasets, geometry, with_faults=True)
    x = solve(A, b, keep, geometry, lam, len(datasets))
    base_rms_all, base_rms_band = rms_all_and_band(datasets, keep, x, masks)
    base = moment_summary(geometry, keep, x)
    stored_rms = ref["rms_m"]["distributed_with_faults"]["all"]
    stored_potency = ref["dike"][f"fault_{sides[0]}_potency_m3"]
    if (abs(base_rms_all - stored_rms) > 1e-6
            or abs(base["dike_volume_km3"]
                   - ref["dike"]["dike_volume_km3"]) > 1e-6
            or abs(base[f"fault_{sides[0]}_potency_m3"] - stored_potency)
            > 1e-6 * stored_potency + 1.0):
        raise SystemExit(f"baseline re-solve does not reproduce "
                         f"{ref_path.name} — rerun the mapped fit, then the "
                         "sweep")
    profile_z, profile_cum = fault_depth_profile(geometry, x, sides[0])
    base_potency = float(base[f"fault_{sides[0]}_potency_m3"])

    # Opening-only floor at the same lambda: what removing the fault
    # entirely costs (the geometry, and so the cached Green's functions,
    # are unchanged).
    A0, b0, keep0 = assemble(datasets, geometry, with_faults=False)
    x0 = solve(A0, b0, keep0, geometry, lam, len(datasets))
    open_rms_all, open_rms_band = rms_all_and_band(datasets, keep0, x0, masks)
    open_summary = moment_summary(geometry, keep0, x0)
    band_gain = open_rms_band - base_rms_band
    if band_gain <= 0:
        raise SystemExit("the fault term does not improve the graben band — "
                         "the sweep's retained-improvement metric is undefined")

    rows = []
    for cap in DEPTH_CAPS_M:
        g = build_geometry(model, fault_trace_offsets=trace_offsets,
                           fault_trace_source="mapped",
                           enabled_fault_sides=sides, fault_depth_cap=cap)
        A_c, b_c, keep_c = assemble(datasets, g, with_faults=True)
        x_c = solve(A_c, b_c, keep_c, g, lam, len(datasets))
        rms_all_c, rms_band_c = rms_all_and_band(datasets, keep_c, x_c, masks)
        s_c = moment_summary(g, keep_c, x_c)
        potency = float(s_c.get(f"fault_{sides[0]}_potency_m3", 0.0))
        rows.append({
            "depth_cap_km": cap / 1e3,
            "fault_patches": int(g["fault_retained_patches"][sides[0]]),
            "fault_max_slip_m": float(
                s_c.get(f"fault_{sides[0]}_max_slip_m", 0.0)),
            "fault_potency_m3": potency,
            "fault_moment_Nm": SHEAR_MODULUS * potency,
            "dike_volume_km3": float(s_c["dike_volume_km3"]),
            "rms_all_cm": 100 * rms_all_c,
            "rms_band_cm": 100 * rms_band_c,
            "band_improvement_retained": (open_rms_band - rms_band_c)
            / band_gain,
            "geometry_sha256": geometry_sha256(g),
        })
        print(f"cap {cap / 1e3:.1f} km: band rms {100 * rms_band_c:.2f} cm, "
              f"potency {potency:.3e} m3, retained "
              f"{rows[-1]['band_improvement_retained']:.2f}")

    reaching = [row for row in rows
                if row["band_improvement_retained"] >= IMPROVEMENT_TARGET]
    headline = reaching[0] if reaching else None

    payload = {
        "nisar_source": nisar_source,
        "lambda_chosen": lam,
        "reference_fit": ref_path.name,
        "reference_fit_sha256": hashlib.sha256(ref_raw).hexdigest(),
        "trace_provenance": {"path": trace_prov["path"],
                             "sha256": trace_prov["sha256"]},
        "modeled_fault_sides": list(sides),
        "shear_modulus_Pa": SHEAR_MODULUS,
        "graben_band": {"perp_m": GRABEN_PERP_M, "along_m": GRABEN_ALONG_M,
                        "n_blocks": n_band},
        "opening_only": {
            "rms_all_cm": 100 * open_rms_all,
            "rms_band_cm": 100 * open_rms_band,
            "dike_volume_km3": float(open_summary["dike_volume_km3"]),
        },
        "baseline": {
            "fault_patches": int(geometry["fault_retained_patches"][sides[0]]),
            "rms_all_cm": 100 * base_rms_all,
            "rms_band_cm": 100 * base_rms_band,
            "dike_volume_km3": float(base["dike_volume_km3"]),
            "fault_max_slip_m": float(base[f"fault_{sides[0]}_max_slip_m"]),
            "fault_potency_m3": base_potency,
            "fault_moment_Nm": SHEAR_MODULUS * base_potency,
            "geometry_sha256": ref["geometry_sha256"],
        },
        "rows": rows,
        "depth_profile": {
            "z_km": [float(v) / 1e3 for v in profile_z],
            "cumulative_potency_fraction": [float(v) / base_potency
                                            for v in profile_cum],
        },
        "potency_fraction_above": {
            f"{z_km:.0f}_km": float(np.interp(z_km * 1e3, profile_z,
                                              profile_cum) / base_potency)
            for z_km in (1.0, 2.0, 3.0)
        },
        "headline": {
            "improvement_target": IMPROVEMENT_TARGET,
            "shallowest_cap_km": (headline["depth_cap_km"]
                                  if headline else None),
            "cap_fault_potency_m3": (headline["fault_potency_m3"]
                                     if headline else None),
            "cap_fault_moment_Nm": (headline["fault_moment_Nm"]
                                    if headline else None),
        },
        "notes": [
            "fault moment = shear_modulus x potency (scalar moment for "
            "shear slip); the dike's muV comparison scale is reported by "
            "the distributed fits",
            "Okada surface displacements depend on the Lame constants only "
            "through alpha=(lambda+mu)/(lambda+2mu); uniformly rescaling "
            "(lambda, mu) changes no prediction and rescales the moments "
            "linearly, so the data do not constrain mu",
            "slip magnitudes remain regularized (positivity + smoothing at "
            "the reused lambda); the sweep varies depth extent only, on the "
            "fixed mapped trace and dip",
        ],
    }
    out = RESULTS / f"fault_depth_mapped{suffix}.json"
    out.write_text(json.dumps(payload, indent=2) + "\n")
    write_md(payload, suffix)
    figure(payload, suffix)
    print(f"wrote {out}")


def write_md(payload: dict, suffix: str) -> None:
    chain = ("GUNW product chain" if payload["nisar_source"] == "gunw"
             else "30-m reprocessed (GSLC) chain")
    band = payload["graben_band"]
    base = payload["baseline"]
    opening = payload["opening_only"]
    lines = [
        f"# Graben-fault depth-extent sweep (mapped trace; {chain})",
        "",
        f"- reference fit: {payload['reference_fit']} (λ = "
        f"{payload['lambda_chosen']:.2f} reused; SHA-256 "
        f"{payload['reference_fit_sha256'][:12]}…); the uncapped baseline "
        "must reproduce it exactly before the sweep runs",
        f"- graben band: within {band['perp_m'] / 1e3:.0f} km of the trace, "
        f"{band['along_m'] / 1e3:.0f} km along strike "
        f"({band['n_blocks']} blocks); the fault term acts here, and the "
        "all-block RMS dilutes it",
        f"- fault moment = μ · Σ area·slip at μ = "
        f"{payload['shear_modulus_Pa'] / 1e9:.0f} GPa. Okada surface "
        "displacements depend on (λ, μ) only through α = (λ+μ)/(λ+2μ), so "
        "a uniform rescaling of both changes no prediction and rescales "
        "the moment linearly: the data do not constrain μ. A "
        "depth-varying μ is outside the half-space model.",
        "",
        "| fault depth extent | patches | max slip | potency | μ·potency "
        "| dike volume | RMS all | RMS band | band improvement retained |",
        "|---|---|---|---|---|---|---|---|---|",
        (f"| none (opening only) | 0 | — | — | — | "
         f"{opening['dike_volume_km3']:.3f} km³ | "
         f"{opening['rms_all_cm']:.2f} cm | {opening['rms_band_cm']:.2f} cm "
         "| 0.00 |"),
    ]
    for row in payload["rows"]:
        lines.append(
            f"| ≤ {row['depth_cap_km']:.1f} km | {row['fault_patches']} | "
            f"{row['fault_max_slip_m']:.2f} m | "
            f"{row['fault_potency_m3'] / 1e6:.1f}×10⁶ m³ | "
            f"{row['fault_moment_Nm'] / 1e17:.1f}×10¹⁷ N·m | "
            f"{row['dike_volume_km3']:.3f} km³ | {row['rms_all_cm']:.2f} cm "
            f"| {row['rms_band_cm']:.2f} cm | "
            f"{row['band_improvement_retained']:.2f} |")
    lines.append(
        f"| trace-limited (≤ 4 km) | {base['fault_patches']} | "
        f"{base['fault_max_slip_m']:.2f} m | "
        f"{base['fault_potency_m3'] / 1e6:.1f}×10⁶ m³ | "
        f"{base['fault_moment_Nm'] / 1e17:.1f}×10¹⁷ N·m | "
        f"{base['dike_volume_km3']:.3f} km³ | {base['rms_all_cm']:.2f} cm | "
        f"{base['rms_band_cm']:.2f} cm | 1.00 |")
    above = payload["potency_fraction_above"]
    lines += [
        "",
        "## Depth reach of the uncapped solution",
        "",
        (f"- cumulative potency fraction above 1 km: {above['1_km']:.2f}; "
         f"above 2 km: {above['2_km']:.2f}; above 3 km: {above['3_km']:.2f} "
         "(piecewise-linear in depth within patches)"),
    ]
    head = payload["headline"]
    if head["shallowest_cap_km"] is not None:
        lines.append(
            f"- shallowest cap retaining ≥ {head['improvement_target']:.0%} "
            f"of the fault term's graben-band improvement: "
            f"{head['shallowest_cap_km']:.1f} km, with potency "
            f"{head['cap_fault_potency_m3'] / 1e6:.1f}×10⁶ m³ "
            f"(μ·potency {head['cap_fault_moment_Nm'] / 1e17:.1f}×10¹⁷ N·m, "
            f"vs {base['fault_moment_Nm'] / 1e17:.1f}×10¹⁷ N·m uncapped)")
    else:
        lines.append(
            f"- no cap in the sweep retains "
            f"{head['improvement_target']:.0%} of the graben-band "
            "improvement; the full trace-limited depth extent is required")
    lines += [
        "",
        "## Scope",
        "",
        "- The sweep truncates the depth extent of the fixed mapped-trace, "
        "fixed-dip fault surface and re-solves at the reused smoothing "
        "weight. It does not test other fault geometries (listric surfaces, "
        "distributed fissuring, the unsupported west trace), and the slip "
        "values remain regularized magnitudes, not resolved slip.",
        "- Truncation removes parameters, so a capped fit can only misfit "
        "more; the question is how little the graben band degrades while "
        "the potency (and so the moment) drops. The dike-volume column "
        "shows how much removed fault signal the opening absorbs.",
        "- The band's RMS level is dominated by axial signal no model term "
        "fits (graben-floor cells, where the loop-closure check also flags "
        "a one-fringe unwrap risk), so differences between rows are "
        "meaningful; the level is not.",
        "- Both input chains run the same sweep "
        "(fault_depth_mapped.md / fault_depth_mapped_gslc.md).",
    ]
    out = RESULTS / f"fault_depth_mapped{suffix}.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


def figure(payload: dict, suffix: str) -> None:
    rows = payload["rows"]
    base = payload["baseline"]
    opening = payload["opening_only"]
    caps = [row["depth_cap_km"] for row in rows]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    x_axis = [0.0] + caps + [4.0]
    band = ([opening["rms_band_cm"]] + [row["rms_band_cm"] for row in rows]
            + [base["rms_band_cm"]])
    ax1.plot(x_axis, band, "o-", color="tab:red", label="RMS, graben band")
    ax1.axhline(opening["rms_band_cm"], color="tab:red", lw=0.8, ls=":",
                alpha=0.7)
    ax1.set_xlabel("fault depth extent cap (km; 0 = no fault, "
                   "4 = trace-limited)")
    ax1.set_ylabel("RMS (cm)", color="tab:red")
    ax1.tick_params(axis="y", labelcolor="tab:red")
    ax1b = ax1.twinx()
    moment = ([0.0] + [row["fault_moment_Nm"] / 1e17 for row in rows]
              + [base["fault_moment_Nm"] / 1e17])
    ax1b.plot(x_axis, moment, "s--", color="tab:blue",
              label="fault μ·potency")
    ax1b.set_ylabel("fault μ·potency (10¹⁷ N·m at 30 GPa)",
                    color="tab:blue")
    ax1b.tick_params(axis="y", labelcolor="tab:blue")
    ax1.set_title("misfit vs fault depth extent")

    profile = payload["depth_profile"]
    ax2.plot(profile["cumulative_potency_fraction"], profile["z_km"],
             color="tab:blue")
    ax2.invert_yaxis()
    ax2.set_xlabel("cumulative potency fraction above depth")
    ax2.set_ylabel("depth (km)")
    ax2.set_title("depth reach of the uncapped slip")
    ax2.grid(alpha=0.3)

    fig.suptitle(f"Graben-fault depth-extent sweep "
                 f"({'GSLC' if suffix else 'GUNW'} chain, mapped trace)",
                 y=1.02)
    fig.tight_layout()
    out = RESULTS / f"fault_depth_mapped{suffix}.png"
    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gslc", action="store_true",
                        help="use the 30-m reprocessed NISAR chain")
    args = parser.parse_args()
    run(nisar_source="gslc" if args.gslc else "gunw")


if __name__ == "__main__":
    main()
