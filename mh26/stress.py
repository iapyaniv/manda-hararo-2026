"""Dike-induced stress change resolved on saved fault surfaces.

If the graben faults slipped because the dike loaded them, the dike's static
stress change must favour inward-dipping normal slip on those surfaces (the
triggered-shallow-faulting test of Xu et al., 2026,
doi:10.1016/j.geog.2026.06.002). This module takes the distributed dike
opening (faults excluded from the source), computes the internal stress
field from the dislocation strain derivatives, and resolves normal and
Coulomb stress change at receivers on the fault planes stored with that fit.
Every receiver is a fault patch centre or a point along a saved trace column.
Receivers closer than 400 m horizontally to the dike plane are excluded as
near-singular.

The result has depth structure and the report keeps it: statements are made
per depth row, not at one depth, and the resolved stresses are computed for
receiver dips 55/65/75 degrees.

Run after `mh26.distributed --mapped`:
    python -m mh26.stress
Use another saved distributed fit without replacing the default report:
    python -m mh26.stress --fit results/distributed_fit_traces_<hash>.json
Self-check only (strain derivatives vs finite differences):
    python -m mh26.stress --selftest

Writes results/stress_change.{json,md,png}.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyrocko.modelling import okada_ext

from .distributed import (
    _strike_vectors,
    geometry_from_json,
    geometry_sha256,
)
from .params import LAME_LAMBDA, REPO, RESULTS, SHEAR_MODULUS

FRICTION = 0.4
DIKE_STANDOFF_M = 400.0     # m horizontal; closer receivers are near-singular
PROFILE_DZ_M = 250.0
RECEIVER_DIPS = (55.0, 65.0, 75.0)
MAIN_DIP = 65.0


def stress_at(patches: np.ndarray, disl: np.ndarray,
              rec: np.ndarray) -> np.ndarray:
    """Stress tensors (n, 3, 3) in NED at the receiver points."""
    res = okada_ext.okada(patches, disl, rec, LAME_LAMBDA, SHEAR_MODULUS,
                          nthreads=0, rotate_sdn=0, stack_sources=1)
    # res[:, 3:12].reshape(3, 3)[i, j] = du_j/dx_i in (north, east, down);
    # checked against finite differences in selftest().
    grad = res[:, 3:12].reshape(-1, 3, 3)
    strain = 0.5 * (grad + np.transpose(grad, (0, 2, 1)))
    trace = np.trace(strain, axis1=1, axis2=2)
    return (LAME_LAMBDA * trace[:, None, None] * np.eye(3)[None]
            + 2 * SHEAR_MODULUS * strain)


def selftest() -> bool:
    """Check the gradient block that every stress tensor rests on.

    Columns 3:12 of the okada_ext result must be du_j/dx_i in NED at the
    receiver; compare them against central differences of the displacement
    columns 0:3 for a tensile patch and off-plane receivers.
    """
    patch = np.array([[0.0, 0.0, 3_000.0, 152.6, 90.0,
                       -900.0, 900.0, -500.0, 500.0]])
    disl = np.array([[0.0, 0.0, 1.0]])
    rec = np.array([[1_500.0, 2_500.0, 1_200.0],
                    [-2_200.0, 1_100.0, 2_600.0],
                    [800.0, -1_900.0, 3_400.0]])

    def run(points: np.ndarray) -> np.ndarray:
        return okada_ext.okada(patch, disl, points, LAME_LAMBDA,
                               SHEAR_MODULUS, nthreads=0, rotate_sdn=0,
                               stack_sources=1)

    grad = run(rec)[:, 3:12].reshape(-1, 3, 3)
    h = 2.0     # small enough for the km-scale field; FD noise grows below ~1 m
    fd = np.empty_like(grad)
    for i in range(3):
        step = np.zeros(3)
        step[i] = h
        fd[:, i, :] = (run(rec + step)[:, :3]
                       - run(rec - step)[:, :3]) / (2.0 * h)
    err = float(np.max(np.abs(fd - grad)))
    scale = float(np.max(np.abs(grad)))
    # An index-order or sign error disagrees at order 1; agreement is
    # roundoff-limited near 1e-4 relative.
    ok = 0.0 < err <= 1e-3 * scale
    print(f"selftest: gradient vs finite difference max |diff| = {err:.2e} "
          f"(max |gradient| = {scale:.2e}, limit 1e-3 relative) "
          f"-> {'OK' if ok else 'FAIL'}")
    return ok


def resolve_normal_slip(stress: np.ndarray, perp: np.ndarray, sign: float,
                        dip_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """(sigma_n tension-positive, Coulomb for normal slip) on one side's plane."""
    dip_r = np.radians(dip_deg)
    h = -sign * perp                      # horizontal down-dip: toward the axis
    down_dip = np.array([h[0] * np.cos(dip_r), h[1] * np.cos(dip_r),
                         np.sin(dip_r)])
    normal = np.array([-h[0] * np.sin(dip_r), -h[1] * np.sin(dip_r),
                       np.cos(dip_r)])
    if normal[2] > 0:                     # orient the normal upward in NED
        normal = -normal                  # = into the hanging wall
    tractions = stress @ normal
    sigma_n = tractions @ normal
    tau_slip = tractions @ down_dip
    return sigma_n, tau_slip + FRICTION * sigma_n


def side_stats(sigma_n: np.ndarray, cfs: np.ndarray) -> dict:
    if len(sigma_n) == 0:
        raise ValueError("receiver set is empty")
    return {
        "n_receivers": int(len(sigma_n)),
        "median_sigma_n_mpa": float(np.median(sigma_n)),
        "fraction_unclamped": float(np.mean(sigma_n > 0)),
        "median_cfs_mpa": float(np.median(cfs)),
        "fraction_cfs_positive": float(np.mean(cfs > 0)),
    }


def fault_patch_receivers(geometry: dict, dip_deg: float) -> dict[str, dict]:
    """Build patch-centre receivers from the geometry saved with a fit."""
    e_c, n_c = geometry["trace_centre"]
    along, perp = _strike_vectors(geometry["strike"])
    s_centres = np.asarray(geometry["fault_s_centres_m"], dtype=float)
    z_edges = np.asarray(geometry["fault_z_edges_m"], dtype=float)
    min_thickness = float(geometry["min_fault_patch_thickness_m"])
    dip_r = np.radians(dip_deg)
    enabled = set(geometry.get("enabled_fault_sides", ("east", "west")))
    out = {}
    for side, sign in (("east", +1.0), ("west", -1.0)):
        if side not in enabled:
            continue
        offsets = np.asarray(geometry["fault_column_offsets_m"][side], dtype=float)
        if offsets.shape != s_centres.shape:
            raise ValueError(f"saved {side} trace columns have the wrong shape")
        points, rows, columns, receiver_offsets = [], [], [], []
        z_cut = np.minimum(z_edges[-1], np.abs(offsets) * np.tan(dip_r))
        for row, (z0, z1) in enumerate(zip(z_edges[:-1], z_edges[1:])):
            for column, (s_c, trace_offset) in enumerate(zip(s_centres, offsets)):
                z1_eff = min(z1, z_cut[column])
                if z1_eff - z0 < min_thickness:
                    continue
                depth = 0.5 * (z0 + z1_eff)
                offset = trace_offset - sign * depth / np.tan(dip_r)
                position = (np.array([n_c, e_c]) + along * s_c + perp * offset)
                points.append([position[0], position[1], depth])
                rows.append(row)
                columns.append(column)
                receiver_offsets.append(offset)
        out[side] = {
            "points": np.asarray(points, dtype=float).reshape(-1, 3),
            "rows": np.asarray(rows, dtype=int),
            "columns": np.asarray(columns, dtype=int),
            "offsets_m": np.asarray(receiver_offsets, dtype=float),
        }
    return out


def fault_profile_receivers(geometry: dict, dip_deg: float,
                            step_m: float = PROFILE_DZ_M) -> dict[str, dict]:
    """Build regular-depth receivers from saved trace columns."""
    e_c, n_c = geometry["trace_centre"]
    along, perp = _strike_vectors(geometry["strike"])
    s_centres = np.asarray(geometry["fault_s_centres_m"], dtype=float)
    max_depth = float(np.asarray(geometry["fault_z_edges_m"])[-1])
    dip_r = np.radians(dip_deg)
    enabled = set(geometry.get("enabled_fault_sides", ("east", "west")))
    out = {}
    for side, sign in (("east", +1.0), ("west", -1.0)):
        if side not in enabled:
            continue
        offsets = np.asarray(geometry["fault_column_offsets_m"][side], dtype=float)
        points, depths, receiver_offsets = [], [], []
        for s_c, trace_offset in zip(s_centres, offsets):
            z_cut = min(max_depth, abs(trace_offset) * np.tan(dip_r))
            for depth in np.arange(step_m / 2, z_cut, step_m):
                offset = trace_offset - sign * depth / np.tan(dip_r)
                position = (np.array([n_c, e_c]) + along * s_c + perp * offset)
                points.append([position[0], position[1], depth])
                depths.append(depth)
                receiver_offsets.append(offset)
        out[side] = {
            "points": np.asarray(points, dtype=float).reshape(-1, 3),
            "depths_m": np.asarray(depths, dtype=float),
            "offsets_m": np.asarray(receiver_offsets, dtype=float),
        }
    return out


def _result_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO.resolve()))
    except ValueError:
        return str(resolved)


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(requested_fit: Path | None = None) -> None:
    if not selftest():
        raise SystemExit("selftest failed — do not trust the stress derivatives")
    default_fit = RESULTS / "distributed_fit_mapped.json"
    if not default_fit.exists():
        default_fit = RESULTS / "distributed_fit.json"
    fit_path = default_fit if requested_fit is None else requested_fit
    payload = json.loads(fit_path.read_text())
    geometry = geometry_from_json(payload["geometry"])
    saved_geometry_hash = payload.get("geometry_sha256")
    if saved_geometry_hash != geometry_sha256(geometry):
        raise ValueError(f"{fit_path}: geometry hash does not match its geometry")
    magnitudes = np.asarray(payload["dislocation_magnitudes_m"], dtype=float)
    if magnitudes.shape != geometry["kinds"].shape:
        raise ValueError(f"{fit_path}: dislocation vector does not match geometry")

    _, perp = _strike_vectors(geometry["strike"])
    dike = geometry["kinds"] == "dike"
    src_patches = geometry["patches"][dike]
    src_disl = geometry["disl"][dike] * magnitudes[dike, None]

    per_dip = {}
    for dip_v in RECEIVER_DIPS:
        per_dip[f"{dip_v:.0f}"] = {}
        for side, receivers in fault_patch_receivers(geometry, dip_v).items():
            sign = +1.0 if side == "east" else -1.0
            keep = np.abs(receivers["offsets_m"]) > DIKE_STANDOFF_M
            rec = receivers["points"][keep]
            rows = receivers["rows"][keep]
            stress = stress_at(src_patches, src_disl, rec)
            sigma_n, cfs = resolve_normal_slip(stress, perp, sign, dip_v)
            sigma_n, cfs = sigma_n / 1e6, cfs / 1e6
            entry = {
                "overall": side_stats(sigma_n, cfs),
                "n_excluded_near_dike": int((~keep).sum()),
                "rows": [],
            }
            for row in sorted(set(rows.tolist())):
                in_row = rows == row
                entry["rows"].append({
                    "row": int(row),
                    "mean_receiver_depth_km": float(np.mean(rec[in_row, 2])) / 1e3,
                    **side_stats(sigma_n[in_row], cfs[in_row]),
                })
            per_dip[f"{dip_v:.0f}"][side] = entry
            overall = entry["overall"]
            print(f"dip {dip_v:.0f}° {side}: median σn "
                  f"{overall['median_sigma_n_mpa']:+.2f} MPa "
                  f"({100 * overall['fraction_unclamped']:.0f}% unclamped), "
                  f"median ΔCFS {overall['median_cfs_mpa']:+.2f} MPa")

    profile = {}
    max_fault_depth = float(np.asarray(geometry["fault_z_edges_m"])[-1])
    z_bins = np.arange(0.0, max_fault_depth + PROFILE_DZ_M, PROFILE_DZ_M)
    for side, receivers in fault_profile_receivers(geometry, MAIN_DIP).items():
        sign = +1.0 if side == "east" else -1.0
        keep = np.abs(receivers["offsets_m"]) > DIKE_STANDOFF_M
        rec = receivers["points"][keep]
        depths = receivers["depths_m"][keep]
        stress = stress_at(src_patches, src_disl, rec)
        sigma_n, cfs = resolve_normal_slip(stress, perp, sign, MAIN_DIP)
        sigma_n, cfs = sigma_n / 1e6, cfs / 1e6
        depth_index = np.digitize(depths, z_bins) - 1
        rowsum = {"depth_km": [], "median_sigma_n_mpa": [],
                  "p16_sigma_n_mpa": [], "p84_sigma_n_mpa": [],
                  "median_cfs_mpa": []}
        for index in range(len(z_bins) - 1):
            selected = depth_index == index
            if selected.sum() < 4:
                continue
            rowsum["depth_km"].append(
                float(z_bins[index] + PROFILE_DZ_M / 2) / 1e3
            )
            rowsum["median_sigma_n_mpa"].append(float(np.median(sigma_n[selected])))
            rowsum["p16_sigma_n_mpa"].append(float(np.percentile(sigma_n[selected], 16)))
            rowsum["p84_sigma_n_mpa"].append(float(np.percentile(sigma_n[selected], 84)))
            rowsum["median_cfs_mpa"].append(float(np.median(cfs[selected])))
        profile[side] = rowsum

    def row_signs(side: str) -> list[str]:
        return [f"{row['mean_receiver_depth_km']:.1f} km: "
                f"{row['median_sigma_n_mpa']:+.1f} MPa"
                for row in per_dip[f"{MAIN_DIP:.0f}"][side]["rows"]]

    opening = magnitudes[dike].reshape(geometry["dike_shape"])
    fit_provenance = {"path": _result_path(fit_path),
                      "sha256": _file_sha256(fit_path)}
    source_label = geometry["fault_trace_source"]
    geometry_label = ("idealized" if source_label == "idealized"
                      else "trace-column")
    interpretation = (
        "Normal stress varies with receiver depth. The shallowest saved row "
        "is mainly unclamped, while most deeper saved rows are clamped in the "
        "65-degree case. This result is conditional on the saved receiver "
        "surfaces. The distributed opening starts at the surface and has no "
        "subsurface hard top, so this report makes no conclusion relative to "
        "the uniform rectangle's fitted top."
    )
    out_payload = {
        "source": f"distributed dike opening ({geometry_label} geometry); "
                  "fault slip excluded from the stress source",
        "source_fit": fit_provenance,
        "geometry_sha256": saved_geometry_hash,
        "source_parameterization": {
            "dike_depth_edges_km": [float(value) / 1e3
                                    for value in geometry["dike_z_edges_m"]],
            "opening_sum_by_depth_row_m": [float(value)
                                            for value in opening.sum(axis=1)],
            "parameterized_top_depth_km":
                float(geometry["dike_z_edges_m"][0]) / 1e3,
            "subsurface_hard_top_depth_km": None,
            "note": "The distributed plane starts at the surface; it has no "
                    "subsurface hard top at the uniform rectangle's fitted "
                    "top depth.",
        },
        "receivers": {
            "placement": "fault patch centres rebuilt only from the geometry "
                         "stored in source_fit; dip-specific columns end at "
                         "the dike plane; receivers within 400 m horizontally "
                         "of the dike plane are excluded",
            "dips_deg": list(RECEIVER_DIPS),
            "rake": "normal (hanging wall down-dip)",
            "friction": FRICTION,
            "standoff_m": DIKE_STANDOFF_M,
        },
        "per_dip": per_dip,
        "depth_profile_dip65": profile,
        "interpretation": interpretation,
    }

    custom = requested_fit is not None and fit_path.resolve() != default_fit.resolve()
    tag = fit_path.stem.removeprefix("distributed_fit").strip("_") or "idealized"
    safe_tag = "".join(char if char.isalnum() or char in "_-" else "_"
                       for char in tag)
    suffix = f"_{safe_tag}" if custom else ""
    json_path = RESULTS / f"stress_change{suffix}.json"
    md_path = RESULTS / f"stress_change{suffix}.md"
    figure_path = RESULTS / f"stress_change{suffix}.png"
    json_path.write_text(json.dumps(out_payload, indent=2, allow_nan=False) + "\n")

    d65 = per_dip[f"{MAIN_DIP:.0f}"]
    lines = [
        "# Dike-induced stress change on the saved fault surfaces", "",
        f"- source fit: `{fit_provenance['path']}` (SHA-256 "
        f"`{fit_provenance['sha256']}`)",
        f"- source: distributed opening with {geometry_label} geometry; "
        "fault slip is excluded from the stress source",
        "- receivers: patch centres from the fit's saved trace columns; "
        "dip-specific truncation at the dike plane; 400 m dike standoff",
        "",
        "## Median normal-stress change by receiver depth (dip 65°; "
        "unclamping positive)", "",
    ]
    for side in d65:
        lines.append(
            f"- {side}: " + "; ".join(row_signs(side))
            + f" (overall {d65[side]['overall']['median_sigma_n_mpa']:+.1f} MPa, "
            f"{100 * d65[side]['overall']['fraction_unclamped']:.0f}% "
            "of receivers unclamped)"
        )
    lines.extend(["", "## Receiver-dip sensitivity", ""])
    for dip_v in RECEIVER_DIPS:
        entries = per_dip[f"{dip_v:.0f}"]
        values = [f"{side} {entry['overall']['median_sigma_n_mpa']:+.1f} / "
                  f"{entry['overall']['median_cfs_mpa']:+.1f} MPa"
                  for side, entry in entries.items()]
        lines.append(f"- dip {dip_v:.0f}°: " + "; ".join(values))
    lines.extend([
        "", interpretation, "", "## Scope", "",
        "- Static elastic half-space stress change from the inverted opening",
        "  only; no fault interaction, pore pressure, or topography.",
        "- Source patches and receiver trace columns come from the saved fit",
        "  geometry. Current source-code geometry constants are not used.",
        "- Receiver dip changes are sensitivity cases. The source geometry",
        "  does not change between these cases.",
    ])
    md_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    colours = {"east": "tab:red", "west": "tab:blue"}
    for ax, key, title in ((axes[0], "median_sigma_n_mpa",
                            "normal stress change (unclamping +)"),
                           (axes[1], "median_cfs_mpa",
                            "Coulomb stress change (normal slip)")):
        for side, values in profile.items():
            ax.plot(values[key], values["depth_km"], "-",
                    color=colours[side], label=f"{side} plane")
            if key == "median_sigma_n_mpa":
                ax.fill_betweenx(values["depth_km"], values["p16_sigma_n_mpa"],
                                 values["p84_sigma_n_mpa"],
                                 color=colours[side], alpha=0.15)
        ax.axvline(0, color="0.4", ls=":", lw=0.8)
        ax.set_xlabel("MPa")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("receiver depth on the fault plane (km)")
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=9)
    fig.suptitle("Stress change from the saved distributed opening, resolved "
                 "on saved trace-column receivers (65° dip; median with "
                 "16–84% band)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(figure_path, dpi=110)
    plt.close(fig)
    print(f"figure: {figure_path}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit", type=Path, metavar="FILE",
                        help="saved distributed-fit JSON; a non-default fit "
                             "gets a separate output suffix")
    parser.add_argument("--selftest", action="store_true")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    if args.selftest:
        raise SystemExit(0 if selftest() else 1)
    main(args.fit)
