"""North/south-half along-rift (U∥) averages on the dike flanks.

Along-rift stretching moves the northern half of the deforming zone
toward the NNW and the southern half toward the SSE — a north/south-
antisymmetric pattern. The decomposition's flank probes are 2-km disks
at the co-event opening peaks (3.8 and 1.9 km south of the fitted
centre) and a full-flank average integrates across the sign change, so
neither statistic measures that contrast. This module averages the
dike-frame U∥ field over (flank × along-strike half) regions anchored
to the fitted dike and compares each region with the model's prediction
over the same cells. The elastic dike itself predicts the antisymmetric
pattern (flank material is dragged toward the tips), so the test is
observed vs model region means, not observed vs zero.

The model column is processed like the observations: the model is
projected into each geometry's LOS at the cell centres, the same
far-field reference plane is removed per dataset, and the same free
three-component solve and rotation are applied. Raw elastic
displacement differs from the referenced column by up to ~1 cm in a
region mean and is kept in the JSON for comparison.

Formal per-cell sigmas assume independent cells and understate a region
mean's error under spatially correlated noise (atmosphere and ionosphere
survive the far-field referencing at these scales). A correlated-noise
floor comes from control regions: boxes with the primary regions'
dimensions tiled over the far field (|t| ≥ 12 km), evaluated on the
observed-minus-model residual. The scatter of their means is the floor
for one region mean, and the scatter of along-strike-adjacent pair
differences is the floor for a north-minus-south differential. It is a
floor, not a complete error model: the control boxes sit almost
entirely east of the dike, and their cross-sensor closure level matches
the west-flank regions but is well below the east-flank regions, so
east-flank errors are larger than the floor by an unquantified amount.

Run after mh26.invert and mh26.s1_tracks:  python -m mh26.along_rift [--gslc]
Writes results/along_rift_halves{,_gslc}.{json,md,png}.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .decompose3d import (
    dike_frame,
    lattice_axes,
    load_cells,
    model_dike_frame,
    solve_cells_3d,
    wright_plane_se_strike,
)
from .distributed import _strike_vectors
from .gunw import remove_farfield_plane
from .okada import forward
from .params import RESULTS

# Same U∥ mask as the decomposition summary figure.
PAR_SIGMA_MAX_CM = 5.0
# Flank band: one lattice cell off the axial strip out to the opening
# peaks' far side; the variant narrows the band and extends past the tips.
FLANK_T_MIN_M = 1_280.0
FLANK_T_MAX_M = 8_000.0
VARIANT_T_MAX_M = 6_000.0
VARIANT_S_MAX_M = 12_000.0
PLACEBO_T_MIN_M = 12_000.0
PLACEBO_MIN_CELLS = 10


def frame_coords(east_c: float, north_c: float,
                 strike_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Dike-frame cell coordinates (s, t) about the fitted dike centre.

    s is along strike with the U∥ sign convention (positive toward
    strike + 180°, NNW here), so the "north" half is s > 0; t is
    perpendicular-left (ENE positive), matching U⊥.
    """
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    along, perp = _strike_vectors(strike_deg)  # (north, east) unit vectors
    de, dn = CE - east_c, CN - north_c
    s = -(along[1] * de + along[0] * dn)
    t = perp[1] * de + perp[0] * dn
    return s, t


def region_selectors(s: np.ndarray, t: np.ndarray, s_max: float,
                     t_min: float, t_max: float) -> dict[str, np.ndarray]:
    flanks = {"WEST": (t <= -t_min) & (t >= -t_max),
              "EAST": (t >= t_min) & (t <= t_max)}
    halves = {"north": (s > 0) & (s < s_max),
              "south": (s < 0) & (s > -s_max)}
    out = {f"{fl} {half}": flanks[fl] & halves[half]
           for fl in flanks for half in halves}
    out["BOTH north"] = (flanks["WEST"] | flanks["EAST"]) & halves["north"]
    out["BOTH south"] = (flanks["WEST"] | flanks["EAST"]) & halves["south"]
    return out


def referenced_model_fields(cells: list[dict], model: np.ndarray) -> dict:
    """The model solved through the observations' own operators.

    The observed fields had a far-field plane removed per dataset before
    the per-cell solve, so raw elastic displacement is not comparable to
    them. This projects the model into each geometry's LOS at the cell
    centres, removes the same far-field plane over that dataset's valid
    cells, and runs the same free solve, giving model fields with the
    observations' coverage, reference, and rotation behaviour.
    """
    centres_e, centres_n = lattice_axes()
    CE, CN = np.meshgrid(centres_e, centres_n)
    disp = forward(model, np.column_stack([CN.ravel(), CE.ravel()]))
    e = disp[:, 1].reshape(CE.shape)
    n = disp[:, 0].reshape(CE.shape)
    u = disp[:, 2].reshape(CE.shape)
    model_cells = []
    for c in cells:
        los = c["e"] * e + c["n"] * n + c["u"] * u
        vmask = np.isfinite(c["los"]) & np.isfinite(c["e"])
        referenced, _ = remove_farfield_plane(
            np.column_stack([CN[vmask], CE[vmask]]), los[vmask])
        los_ref = np.full(los.shape, np.nan)
        los_ref[vmask] = referenced
        model_cells.append({**c, "los": los_ref})
    return solve_cells_3d(model_cells)


def region_stats(sel: np.ndarray, upar_cm: np.ndarray, spar_cm: np.ndarray,
                 mpar_cm: np.ndarray,
                 mpar_raw_cm: np.ndarray | None = None) -> dict | None:
    n = int(sel.sum())
    if n == 0:
        return None
    out = {
        "n_cells": n,
        "mean_cm": float(np.mean(upar_cm[sel])),
        "median_cm": float(np.median(upar_cm[sel])),
        # Independent-cell bound; the control-box scatter is the floor
        # actually quoted for a region mean.
        "formal_sigma_mean_cm": float(
            np.sqrt(np.sum(spar_cm[sel] ** 2)) / n),
        "model_mean_cm": float(np.mean(mpar_cm[sel])),
    }
    if mpar_raw_cm is not None:
        out["model_raw_mean_cm"] = float(np.mean(mpar_raw_cm[sel]))
    return out


def max_single_cell_influence(upar_cm: np.ndarray, north: np.ndarray,
                              south: np.ndarray) -> float:
    """Largest change of the north-minus-south differential when any one
    cell is removed — a fragility measure for small, ragged regions."""
    largest = 0.0
    for sel, sign in ((north, 1.0), (south, -1.0)):
        values = upar_cm[sel]
        if len(values) < 2:
            continue
        mean = values.mean()
        shifted = (values.sum() - values) / (len(values) - 1)
        largest = max(largest, float(np.max(np.abs(shifted - mean))))
    return largest


def placebo_scatter(resid_cm: np.ndarray, valid: np.ndarray, s: np.ndarray,
                    t: np.ndarray, box_s_m: float,
                    box_t_m: float) -> dict | None:
    """Far-field box means of the residual, in the real regions' geometry.

    Boxes tile non-overlapping (s, t) strips with |t| ≥ PLACEBO_T_MIN_M.
    Pairs are along-strike-adjacent boxes in the same strip, the geometry
    of a north-minus-south differential. Where coverage puts the boxes
    (almost all on one side here) is reported, because the scatter only
    represents regions whose noise regime the boxes share.
    """
    t_extent = float(np.max(np.abs(t)))
    strips = []
    edge = PLACEBO_T_MIN_M
    while edge + box_t_m <= t_extent:
        strips += [(edge, edge + box_t_m), (-edge - box_t_m, -edge)]
        edge += box_t_m
    s_extent = float(np.max(np.abs(s)))
    n_steps = int(2 * s_extent / box_s_m)
    means, pairs = [], []
    n_east = n_west = 0
    box_cells = np.zeros(resid_cm.shape, dtype=bool)
    for t_lo, t_hi in strips:
        strip_means = []
        for k in range(n_steps):
            s_lo = -s_extent + k * box_s_m
            sel = valid & (s >= s_lo) & (s < s_lo + box_s_m) \
                & (t >= t_lo) & (t < t_hi)
            if sel.sum() >= PLACEBO_MIN_CELLS:
                strip_means.append(float(np.mean(resid_cm[sel])))
                box_cells |= sel
                n_east, n_west = n_east + (t_lo > 0), n_west + (t_lo < 0)
            else:
                strip_means.append(None)
        means += [m for m in strip_means if m is not None]
        pairs += [b - a for a, b in zip(strip_means, strip_means[1:])
                  if a is not None and b is not None]
    if len(means) < 5 or len(pairs) < 5:
        return None
    return {
        "box_dimensions_km": [box_s_m / 1e3, box_t_m / 1e3],
        "n_boxes": len(means),
        "n_boxes_east_of_dike": n_east,
        "n_boxes_west_of_dike": n_west,
        "mean_sigma_cm": float(np.std(means, ddof=1)),
        "pair_differential_sigma_cm": float(np.std(pairs, ddof=1)),
        "n_adjacent_pairs": len(pairs),
        "box_mean_extrema_cm": [float(min(means)), float(max(means))],
        "box_cells": box_cells,
    }


def differential(regions: dict, flank: str,
                 pair_sigma_cm: float | None) -> dict | None:
    north = regions.get(f"{flank} north")
    south = regions.get(f"{flank} south")
    if north is None or south is None:
        return None
    observed = north["mean_cm"] - south["mean_cm"]
    model = north["model_mean_cm"] - south["model_mean_cm"]
    out = {
        "observed_cm": observed,
        "model_cm": model,
        "formal_sigma_cm": float(np.hypot(
            north["formal_sigma_mean_cm"], south["formal_sigma_mean_cm"])),
    }
    if pair_sigma_cm:
        out["control_pair_sigma_cm"] = pair_sigma_cm
        out["detection_sigma_vs_floor"] = observed / pair_sigma_cm
        out["model_departure_sigma_vs_floor"] = \
            (observed - model) / pair_sigma_cm
    return out


def fmt(value: float) -> str:
    return f"{value:+.1f}"


def figure(upar_cm: np.ndarray, selectors: dict, model: np.ndarray,
           regions: dict, diffs: dict, suffix: str) -> None:
    centres_e, centres_n = lattice_axes()
    fig, (ax_map, ax_bar) = plt.subplots(
        1, 2, figsize=(12.5, 5.4), width_ratios=[1.05, 1])

    extent = [centres_e[0] / 1e3, centres_e[-1] / 1e3,
              centres_n[0] / 1e3, centres_n[-1] / 1e3]
    cmap = matplotlib.colormaps["RdBu_r"].copy()
    cmap.set_bad("#b3b3b3")
    im = ax_map.imshow(upar_cm, origin="lower", cmap=cmap, extent=extent,
                       vmin=-15, vmax=15)
    east_c, north_c, strike = model[0], model[1], model[2]
    along, perp = _strike_vectors(float(strike))
    half = model[3] / 2
    ax_map.plot([(east_c - along[1] * half) / 1e3,
                 (east_c + along[1] * half) / 1e3],
                [(north_c - along[0] * half) / 1e3,
                 (north_c + along[0] * half) / 1e3],
                color="k", lw=1.2)
    u_s = np.array([-along[1], -along[0]])  # (east, north), NNW positive
    u_t = np.array([perp[1], perp[0]])
    corners = {
        "WEST north": (0, half, -FLANK_T_MAX_M, -FLANK_T_MIN_M),
        "WEST south": (-half, 0, -FLANK_T_MAX_M, -FLANK_T_MIN_M),
        "EAST north": (0, half, FLANK_T_MIN_M, FLANK_T_MAX_M),
        "EAST south": (-half, 0, FLANK_T_MIN_M, FLANK_T_MAX_M),
    }
    colors = {"north": "#1b7837", "south": "#762a83"}
    for name, (s_lo, s_hi, t_lo, t_hi) in corners.items():
        ring = [(s_lo, t_lo), (s_hi, t_lo), (s_hi, t_hi), (s_lo, t_hi),
                (s_lo, t_lo)]
        xy = np.array([(east_c + sv * u_s[0] + tv * u_t[0],
                        north_c + sv * u_s[1] + tv * u_t[1])
                       for sv, tv in ring]) / 1e3
        ax_map.plot(xy[:, 0], xy[:, 1], color=colors[name.split()[1]],
                    lw=1.4)
    ax_map.set_xlim((east_c - 24_000) / 1e3, (east_c + 24_000) / 1e3)
    ax_map.set_ylim((north_c - 24_000) / 1e3, (north_c + 24_000) / 1e3)
    ax_map.set_xlabel("UTM 37N easting (km)")
    ax_map.set_ylabel("UTM 37N northing (km)")
    ax_map.set_title("U∥ (σ∥ ≤ 5 cm), fitted dike, and the half-flank "
                     "regions\n(green north halves, purple south halves)",
                     fontsize=10)
    fig.colorbar(im, ax=ax_map, shrink=0.85, label="U∥ (cm, NNW positive)")

    names = ["WEST north", "WEST south", "EAST north", "EAST south"]
    xpos = np.arange(len(names))
    observed = [regions[n]["mean_cm"] for n in names]
    modeled = [regions[n]["model_mean_cm"] for n in names]
    errors = [regions[n].get("control_floor_sigma_cm",
                             regions[n]["formal_sigma_mean_cm"])
              for n in names]
    ax_bar.bar(xpos - 0.17, observed, width=0.34, yerr=errors, capsize=3,
               color="#2166ac", label="observed")
    ax_bar.bar(xpos + 0.17, modeled, width=0.34, color="#b2182b",
               alpha=0.75, label="model (same processing)")
    ax_bar.axhline(0, color="k", lw=0.6)
    ax_bar.set_xticks(xpos)
    ax_bar.set_xticklabels([n.replace(" ", "\n") for n in names])
    ax_bar.set_ylabel("region-mean U∥ (cm, NNW positive)")
    lines = ["Region means; error bars are the far-field control floor"]
    for flank in ("WEST", "EAST"):
        d = diffs[flank]
        if "control_pair_sigma_cm" in d:
            lines.append(
                f"{flank} N−S: {d['observed_cm']:+.1f} ± "
                f"{d['control_pair_sigma_cm']:.1f} cm "
                f"(model {d['model_cm']:+.1f})")
    ax_bar.set_title("\n".join(lines), fontsize=9)
    ax_bar.legend(fontsize=9)
    fig.tight_layout()
    out = RESULTS / f"along_rift_halves{suffix}.png"
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
    east_c, north_c, strike, length = (float(model[0]), float(model[1]),
                                       float(model[2]), float(model[3]))

    cells = load_cells(model, decay, nisar_source=nisar_source)
    fields_3d = solve_cells_3d(cells)
    frame = dike_frame(fields_3d, strike)
    fields_model = referenced_model_fields(cells, model)
    model_frame = dike_frame(fields_model, strike)

    upar_cm = 100.0 * frame["par"]
    spar_cm = 100.0 * frame["sigma_par"]
    mpar_cm = 100.0 * model_frame["par"]
    mpar_raw_cm = 100.0 * model_dike_frame(model)["par"]
    valid = np.isfinite(upar_cm) & (spar_cm <= PAR_SIGMA_MAX_CM)
    s, t = frame_coords(east_c, north_c, strike)

    half_length = length / 2
    placebo = placebo_scatter(np.where(valid, upar_cm - mpar_cm, np.nan),
                              valid, s, t, half_length,
                              FLANK_T_MAX_M - FLANK_T_MIN_M)
    box_cells = placebo.pop("box_cells") if placebo else None
    mean_sigma = placebo["mean_sigma_cm"] if placebo else None
    pair_sigma = placebo["pair_differential_sigma_cm"] if placebo else None

    variants = {}
    for key, s_max, t_min, t_max in (
            ("primary", half_length, FLANK_T_MIN_M, FLANK_T_MAX_M),
            ("variant", VARIANT_S_MAX_M, FLANK_T_MIN_M, VARIANT_T_MAX_M)):
        selectors = region_selectors(s, t, s_max, t_min, t_max)
        regions = {name: region_stats(sel & valid, upar_cm, spar_cm,
                                      mpar_cm, mpar_raw_cm)
                   for name, sel in selectors.items()}
        regions = {name: r for name, r in regions.items() if r}
        if key == "primary" and mean_sigma:
            for name, r in regions.items():
                if not name.startswith("BOTH"):
                    r["control_floor_sigma_cm"] = mean_sigma
        # The control floor applies to single contiguous regions; the
        # BOTH statistic spans both flanks with unequal weights and gets
        # no floor of its own.
        diffs = {flank: differential(
                     regions, flank,
                     pair_sigma if key == "primary" and flank != "BOTH"
                     else None)
                 for flank in ("WEST", "EAST", "BOTH")}
        for flank, d in diffs.items():
            if d and flank != "BOTH":
                d["max_single_cell_influence_cm"] = max_single_cell_influence(
                    upar_cm, selectors[f"{flank} north"] & valid,
                    selectors[f"{flank} south"] & valid)
        variants[key] = {
            "s_max_m": s_max, "t_min_m": t_min, "t_max_m": t_max,
            "regions": regions,
            "north_minus_south_cm": {k: v for k, v in diffs.items() if v},
        }
        if key == "primary":
            primary_selectors, primary_regions = selectors, regions
            primary_diffs = {k: v for k, v in diffs.items() if v}

    # Cross-sensor closure (weighted residual, σ units) tells whether the
    # control boxes share the target regions' noise regime.
    closure = {}
    if box_cells is not None:
        regions_closure = {"controls": box_cells}
        regions_closure.update(
            {name: primary_selectors[name] & valid
             for name in ("WEST north", "WEST south",
                          "EAST north", "EAST south")})
        for name, sel in regions_closure.items():
            values = fields_3d["closure"][sel & np.isfinite(fields_3d["closure"])]
            closure[name] = float(np.median(values)) if len(values) else None

    # U⊥ leaks into U∥ with the same sign along a whole flank when the
    # strike is wrong, so the differential should survive the alternative
    # strike even though individual region means move by several cm. The
    # model is rotated into the same alternative frame.
    alt = None
    alt_strike = wright_plane_se_strike(strike)
    if alt_strike is not None:
        alt_frame = dike_frame(fields_3d, alt_strike)
        alt_model = dike_frame(fields_model, alt_strike)
        alt_u = 100.0 * alt_frame["par"]
        alt_s = 100.0 * alt_frame["sigma_par"]
        alt_valid = np.isfinite(alt_u) & (alt_s <= PAR_SIGMA_MAX_CM)
        alt_regions = {
            name: region_stats(sel & alt_valid, alt_u, alt_s,
                               100.0 * alt_model["par"])
            for name, sel in primary_selectors.items()}
        alt_regions = {k: v for k, v in alt_regions.items() if v}
        alt = {"strike_deg": float(alt_strike),
               "north_minus_south_cm": {
                   flank: d for flank in ("WEST", "EAST")
                   if (d := differential(alt_regions, flank, None))}}

    payload = {
        "nisar_source": nisar_source,
        "parameters": {
            "par_sigma_max_cm": PAR_SIGMA_MAX_CM,
            "dike_centre_utm": [east_c, north_c],
            "strike_deg": strike,
            "half_length_m": half_length,
            "placebo_t_min_m": PLACEBO_T_MIN_M,
            "placebo_min_cells": PLACEBO_MIN_CELLS,
        },
        "variants": variants,
        "control_noise": placebo,
        "closure_median_sigma_units": closure,
        "strike_sensitivity": alt,
        "notes": [
            "regions are anchored to the fitted dike: s along strike "
            "(NNW positive, the U∥ sign convention), t perpendicular "
            "(ENE positive); north half 0 < s < L/2, south half "
            "-L/2 < s < 0",
            "model_mean_cm is the uniform-opening fit processed like the "
            "observations (per-geometry LOS, per-dataset far-field "
            "reference plane, same solve and rotation) over the same "
            "cells; model_raw_mean_cm is the unreferenced elastic "
            "displacement, kept for comparison",
            "formal sigmas assume independent cells; the control-box "
            "scatter is a correlated-noise floor, not a complete error "
            "model — the boxes sit almost entirely east of the dike and "
            "their closure matches the west-flank regions but not the "
            "east-flank regions (closure_median_sigma_units), so "
            "east-flank errors exceed the floor by an unquantified "
            "amount, and the sigma-unit statements are relative to the "
            "floor",
            "a strike error leaks U⊥ into U∥ with one sign per flank, "
            "shifting both halves of a flank together; the "
            "north-minus-south differential is the strike-robust "
            "statistic (see strike_sensitivity)",
            "max_single_cell_influence_cm is the largest change of a "
            "differential when any single cell is removed; the EAST "
            "differential is fragile by this measure, the WEST one is "
            "not",
        ],
    }
    (RESULTS / f"along_rift_halves{suffix}.json").write_text(
        json.dumps(payload, indent=2, default=float) + "\n")

    region_order = ["WEST north", "WEST south", "EAST north", "EAST south",
                    "BOTH north", "BOTH south"]
    lines = [
        "# Along-rift (U∥) north/south-half averages on the flanks"
        + (" (30-m reprocessed NISAR inputs)" if suffix else ""), "",
        "Along-rift stretching moves the northern half of the deforming "
        "zone NNW and",
        "the southern half SSE. The flank probes are local 2-km disks at "
        "the opening",
        "peaks (3.8 and 1.9 km south of the fitted centre) and a "
        "full-flank average",
        "integrates across the sign change, so neither measures that "
        "contrast.",
        "These are half-flank averages of the dike-frame U∥ field "
        f"(σ∥ ≤ {PAR_SIGMA_MAX_CM:.0f} cm, the",
        "summary-figure mask), against the fitted model passed through "
        "the same",
        "processing (per-geometry LOS, per-dataset far-field reference "
        "plane, same",
        "solve and rotation) over the same cells. The elastic dike "
        "itself predicts",
        "the antisymmetric pattern, so the test is observed vs model, "
        "not observed",
        "vs zero.", "",
        f"Regions: |s| < L/2 = {half_length / 1e3:.1f} km from the fitted "
        f"centre, {FLANK_T_MIN_M / 1e3:.1f}–{FLANK_T_MAX_M / 1e3:.0f} km "
        "off axis.", "",
        "| region | n | observed (cm) | model (cm) |",
        "|---|---|---|---|",
    ]
    for name in region_order:
        r = primary_regions.get(name)
        if r is None:
            continue
        err = (f" ± {r['control_floor_sigma_cm']:.1f}"
               if "control_floor_sigma_cm" in r else "")
        lines.append(f"| {name} | {r['n_cells']} | {fmt(r['mean_cm'])}{err} "
                     f"(median {fmt(r['median_cm'])}) | "
                     f"{fmt(r['model_mean_cm'])} |")
    lines.append("")
    for flank in ("WEST", "EAST"):
        d = primary_diffs.get(flank)
        if d is None or "control_pair_sigma_cm" not in d:
            continue
        lines.append(
            f"- **{flank} north − south: {fmt(d['observed_cm'])} ± "
            f"{d['control_pair_sigma_cm']:.1f} cm** (model "
            f"{fmt(d['model_cm'])}); relative to the control floor: "
            f"{d['detection_sigma_vs_floor']:.1f} σ from zero, "
            f"{d['model_departure_sigma_vs_floor']:+.1f} σ from the "
            f"model; largest single-cell influence "
            f"{d['max_single_cell_influence_cm']:.1f} cm")
    both = primary_diffs.get("BOTH")
    if both:
        lines.append(
            f"- BOTH flanks north − south: {fmt(both['observed_cm'])} cm "
            f"(model {fmt(both['model_cm'])}) — no control floor is "
            "assigned: the statistic spans the two flanks with unequal "
            "weights")
    if placebo:
        lines += ["",
                  f"- control floor: {placebo['n_boxes']} far-field boxes "
                  f"of the regions' dimensions "
                  f"({placebo['box_dimensions_km'][0]:.1f} × "
                  f"{placebo['box_dimensions_km'][1]:.1f} km) on the "
                  f"observed − model residual give a region-mean scatter "
                  f"of {placebo['mean_sigma_cm']:.1f} cm and an "
                  f"adjacent-pair differential scatter of "
                  f"{placebo['pair_differential_sigma_cm']:.1f} cm "
                  f"({placebo['n_adjacent_pairs']} pairs). "
                  f"{placebo['n_boxes_east_of_dike']} of "
                  f"{placebo['n_boxes']} boxes lie east of the dike; "
                  "median cross-sensor closure is "
                  f"{closure['controls']:.2f} σ in the boxes vs "
                  f"{closure['WEST north']:.2f}/"
                  f"{closure['WEST south']:.2f} σ in the west regions "
                  f"and {closure['EAST north']:.2f}/"
                  f"{closure['EAST south']:.2f} σ in the east regions, "
                  "so the floor represents the west-flank statistics and "
                  "understates the east-flank errors; the formal "
                  "independent-cell errors are several times smaller "
                  "still and are not quoted"]
    variant = variants["variant"]
    lines.append(
        "- region sensitivity (|s| < "
        f"{variant['s_max_m'] / 1e3:.0f} km, "
        f"{variant['t_min_m'] / 1e3:.1f}–{variant['t_max_m'] / 1e3:.0f} km "
        "off axis): north − south "
        + "; ".join(
            f"{flank} {fmt(d['observed_cm'])} cm (model {fmt(d['model_cm'])})"
            for flank, d in variant["north_minus_south_cm"].items()
            if flank != "BOTH"))
    if alt:
        lines.append(
            f"- strike sensitivity ({alt['strike_deg']:.1f}° from the "
            "southernmost 2005 segment): north − south "
            + "; ".join(
                f"{flank} {fmt(d['observed_cm'])} cm"
                for flank, d in alt["north_minus_south_cm"].items())
            + " — the differential survives the strike change that moves "
            "individual probes by several cm")
    east = primary_diffs.get("EAST")
    lines += [
        "",
        "## Notes",
        "",
        "- The free three-component solve observes the near-north "
        "direction weakly (median σ∥ ≈ 4 cm per cell), and the cells are "
        "not independent; the control-box scatter is a correlated-noise "
        "floor, not a complete error model, so the σ statements above "
        "are relative to that floor. Model-parameter uncertainty and "
        "cross-sensor systematics are not propagated.",
        "- The WEST differential is the robust result: the control "
        "floor's closure level matches the west regions, and no single "
        "cell changes the differential by more than "
        f"{primary_diffs['WEST']['max_single_cell_influence_cm']:.1f} cm.",
        "- The EAST differential is fragile and does not establish the "
        "reversal: both east halves average positive, the south region "
        "is the raggedest (fewest cells, mean and median farthest "
        "apart), one cell moves the differential by "
        f"{east['max_single_cell_influence_cm']:.1f} cm, and the "
        "east-flank closure is well above the control level (the "
        "east-flank cross-sensor disagreement of "
        "results/decomposition_joint*.md).",
        "- A strike error leaks U⊥ into U∥ with one sign per flank, "
        "shifting both halves together, so the north-minus-south "
        "differential is the strike-robust statistic.",
        "- The model column shares each region's exact cells and the "
        "observations' reference treatment; the unreferenced elastic "
        "means (model_raw_mean_cm in the JSON) differ by up to about "
        "1 cm per region.",
    ]
    (RESULTS / f"along_rift_halves{suffix}.md").write_text(
        "\n".join(lines) + "\n")
    print("\n".join(lines))

    figure(np.where(valid, upar_cm, np.nan), primary_selectors, model,
           primary_regions, primary_diffs, suffix)


if __name__ == "__main__":
    import sys

    main(nisar_source="gslc" if "--gslc" in sys.argv else "gunw")
