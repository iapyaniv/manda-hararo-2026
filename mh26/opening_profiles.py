"""Depth-averaged opening per intrusion, after Hamling et al. (2009).

Fig. 13 of Hamling et al. (2009) shows the top-10-km depth-averaged
opening of each Dabbahu-Manda Hararo intrusion versus distance along
the segment, with a bottom panel stacking the per-intrusion openings
into the cumulative profile. This module renders that figure updated
through the whole sequence: the September 2005 megadike (Ian Hamling's
10-km model, decoded in mh26.hamling), the thirteen canonical
2006-2010 intrusions of his 2026-09-17 mapping, and the July 2026 dike
from this repository's fixed-plane fit, all on the along-plane s axis
of the mh26.wright2006 contract.

The depth average is over 0-10 km everywhere, so on the native grid a
profile integrates exactly to the epoch's above-10-km volume
(cross-checked against mh26.hamling on every run): average(s) x 10 km
x column width = column volume. The 2026 model plane spans 0-10 km,
so its average needs no truncation. The drawn curves interpolate
linearly between column centres on the projected s axis for display
and are not volume-preserving (the payload's display block records
every drawn area next to its native volume; worst event ~6%), so
quoted volumes always come from the native-grid tables. The delivered
inversions carry no uncertainties, so unlike the published figure the
panels have no error bars; the comparability caveats of
`results/uperp_three_epoch*.md` (different regularization and greens
functions, floored post-2005 files) apply unchanged.

All input decoding, the canonical-13 selection with its duplicate-pair
guards, and the grid registration are reused from
mh26.hamling.compute_report, so a changed delivery or a stale fit
fails here exactly as it fails there.

Report + figures:  python -m mh26.opening_profiles           (GUNW chain)
                   python -m mh26.opening_profiles --gslc    (30-m GSLC chain)
Outputs: results/opening_profiles[_gslc].{json,md,png,pdf} and the
standalone cumulative panel
results/opening_profiles_cumulative[_gslc].{png,pdf}.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .distributed import _display_source_path
from .hamling import (DATA_DIR, EPOCH_TAGS, PRIMARY_EXCLUDED, RECEIVED,
                      SOURCE_URL, UPDATED_SEPT05, _find_fit, _sha256,
                      compute_report as _three_epoch_report)
from .params import RESULTS

AVERAGE_DEPTH_KM = 10.0

# Panel labels are the intrusion months of the Hamling et al. (2009,
# 2010) sequence, not the delivery's file tags: the july09 file is the
# re-estimate of the June 2009 intrusion (its defunct earlier version
# was tagged june09). Exact onset dates for captions belong to the
# source catalogue.
EVENT_LABELS = {
    "sept05": "September 2005",
    "june06": "June 2006",
    "july06": "July 2006",
    "sept06": "September 2006",
    "dec06": "December 2006",
    "jan07": "January 2007",
    "aug07": "August 2007",
    "nov07": "November 2007",
    "apr08": "April 2008",
    "july08": "July 2008",
    "oct08": "October 2008",
    "feb09": "February 2009",
    "july09": "June 2009",
    "may10": "May 2010",
}

# One color per intrusion, identical in its own panel and in the
# cumulative stack. The first eight follow Fig. 13 of Hamling et
# al. (2009) — orange 2005, blue June 2006, red July 2006, green
# September 2006, black December 2006, yellow January 2007, cyan
# August 2007, pink November 2007 — and the 2008-2010 and 2026 events
# extend the set with hues chosen so that neighbours in the stack stay
# separable under color-vision deficiency (checked with the OKLab
# adjacent-pair validator; identity never rests on color alone because
# every event has its labelled panel).
EVENT_COLORS = {
    "sept05": "#ef8a1c",
    "june06": "#2c4d9e",
    "july06": "#cf3423",
    "sept06": "#6abf5e",
    "dec06": "#2b2b2b",
    "jan07": "#d4af1a",
    "aug07": "#43b7d8",
    "nov07": "#c9439c",
    "apr08": "#9a5b16",
    "july08": "#6f52c1",
    "oct08": "#0f9180",
    "feb09": "#97a11f",
    "july09": "#6fa8e8",
    "may10": "#345995",
    "epoch_2026": "#a50f15",
}
LABEL_2026 = "July 2026"


DISPLAY_STEP_KM = 0.02


def _profile_stats(s_km: np.ndarray, profile: np.ndarray) -> dict:
    peak = int(np.argmax(profile))
    return {
        "max_average_m": float(profile.max()),
        "s_at_max_km": float(s_km[peak]),
    }


def _extend_to_edges(centres: np.ndarray, values: np.ndarray,
                     lo: float, hi: float) -> tuple[np.ndarray, np.ndarray]:
    """Carry the end columns flat to their true cell boundaries.

    An end column is a constant-value cell; stopping the curve at its
    centre would drop half the cell and draw a false boundary half a
    column inside the model's real extent.
    """
    return (np.concatenate(([lo], centres, [hi])),
            np.concatenate(([values[0]], values, [values[-1]])))


def _resample(x: np.ndarray, centres: np.ndarray, values: np.ndarray,
              lo: float, hi: float) -> np.ndarray:
    """Interpolate a column profile on the fine grid x, zero outside.

    Linear between column centres on the projected axis and flat from
    the end centres to the model edges lo/hi, so the stack and the
    per-event panels draw one shape. This display resampling is not
    volume-preserving: the payload's display block records the drawn
    areas next to the native-grid volumes, which are the quoted
    numbers.
    """
    c, v = _extend_to_edges(centres, values, lo, hi)
    return np.interp(x, c, v, left=0.0, right=0.0)


def compute_report(chain: str, fit_path: Path, data_dir: Path) -> tuple:
    """Rebuild the full report payload deterministically.

    mh26.validate calls this against the canonical fit and the tracked
    delivery and requires exact equality with the stored JSON.
    """
    three_epoch, context = _three_epoch_report(chain, fit_path, data_dir)
    delivery = context["delivery"]
    registration = context["registration"]
    epochs = context["epochs"]
    cum_tags = context["cum_tags"]
    geometry = context["geometry"]
    fit = context["fit"]

    s_km = registration["s_km"]
    widths = delivery["columns_width_km"]

    events = []
    fields = {"sept05": context["updated_field"]}
    for tag in cum_tags:
        fields[tag] = delivery["openings_m"][tag][:, :10]
    for tag in ("sept05", *cum_tags):
        profile = fields[tag].mean(axis=1)
        stats_key = "sept05_updated" if tag == "sept05" else tag
        volume = epochs[stats_key]["volume_above_10km_km3"]
        derived = float((profile * widths).sum()) * AVERAGE_DEPTH_KM * 1e-3
        if not np.isclose(derived, volume, rtol=1e-9, atol=0.0):
            raise ValueError(
                f"{tag}: profile integral {derived} km3 does not reproduce "
                f"the three-epoch above-10-km volume {volume} km3")
        events.append({
            "tag": tag,
            "label": EVENT_LABELS[tag],
            "file": (Path(UPDATED_SEPT05).name if tag == "sept05"
                     else f"dyke_only_{tag}_ex.txt"),
            "color": EVENT_COLORS[tag],
            "volume_above_10km_km3": volume,
            "average_opening_m": profile.tolist(),
            **_profile_stats(s_km, profile),
        })

    nz, ns = geometry["dike_shape"]
    opening_2026 = np.asarray(fit["dislocation_magnitudes_m"])[:nz * ns]
    opening_2026 = opening_2026.reshape(nz, ns)
    z_edges = np.asarray(geometry["dike_z_edges_m"], dtype=float) / 1e3
    if not (z_edges[0] == 0.0 and z_edges[-1] == AVERAGE_DEPTH_KM
            and np.allclose(np.diff(z_edges), 1.0)):
        raise ValueError("the 2026 plane no longer spans 0-10 km in one-km "
                         "rows; the depth average needs truncation handling")
    s_edges = np.asarray(geometry["dike_s_edges_m"], dtype=float) / 1e3
    s_2026 = 0.5 * (s_edges[:-1] + s_edges[1:])
    profile_2026 = opening_2026.mean(axis=0)
    volume_2026 = epochs["epoch_2026"]["volume_above_10km_km3"]
    derived = float((profile_2026 * np.diff(s_edges)).sum()) \
        * AVERAGE_DEPTH_KM * 1e-3
    if not np.isclose(derived, volume_2026, rtol=1e-9, atol=0.0):
        raise ValueError(
            f"2026 profile integral {derived} km3 does not reproduce the "
            f"three-epoch above-10-km volume {volume_2026} km3")
    epoch_2026 = {
        "label": LABEL_2026,
        "color": EVENT_COLORS["epoch_2026"],
        "reference_fit": three_epoch["epoch_2026"]["reference_fit"],
        "reference_fit_sha256":
            three_epoch["epoch_2026"]["reference_fit_sha256"],
        "geometry_sha256": three_epoch["epoch_2026"]["geometry_sha256"],
        "s_centres_km": s_2026.tolist(),
        "s_edges_km": s_edges.tolist(),
        "volume_above_10km_km3": volume_2026,
        "average_opening_m": profile_2026.tolist(),
        **_profile_stats(s_2026, profile_2026),
    }

    total = float(sum(e["volume_above_10km_km3"] for e in events)
                  + volume_2026)

    # The figures interpolate between column centres on the projected
    # axis, which redistributes area; record every drawn area beside
    # its native volume so the deviation is a stated, validated number
    # rather than a silent property of the rendering.
    x_max = max(registration["span_km"][1] + 0.6,
                registration["plane_span_km"][1])
    x = np.arange(0.0, x_max + 0.01, DISPLAY_STEP_KM)
    s_nwse = s_km[::-1]
    widths_nwse = widths[::-1]
    extent_ham = (float(s_nwse[0] - 0.5 * widths_nwse[0]),
                  float(s_nwse[-1] + 0.5 * widths_nwse[-1]))
    extent_2026 = (float(s_edges[0]), float(s_edges[-1]))
    drawn_area = {}
    stack = np.zeros_like(x)
    worst = ("", 0.0)
    for event in events:
        y = _resample(x, s_nwse,
                      np.asarray(event["average_opening_m"])[::-1],
                      *extent_ham)
        drawn_area[event["tag"]] = \
            float(np.trapezoid(y, x)) * AVERAGE_DEPTH_KM * 1e-3
        relative = abs(drawn_area[event["tag"]]
                       / event["volume_above_10km_km3"] - 1.0)
        if relative > worst[1]:
            worst = (event["tag"], relative)
        stack += y
    y_2026 = _resample(x, s_2026, profile_2026, *extent_2026)
    drawn_area["epoch_2026"] = \
        float(np.trapezoid(y_2026, x)) * AVERAGE_DEPTH_KM * 1e-3
    relative = abs(drawn_area["epoch_2026"] / volume_2026 - 1.0)
    if relative > worst[1]:
        worst = ("epoch_2026", relative)
    stack += y_2026
    drawn_area["cumulative"] = \
        float(np.trapezoid(stack, x)) * AVERAGE_DEPTH_KM * 1e-3
    se_2005 = registration["se_limit_2005_km"]
    y_sept05 = _resample(x, s_nwse,
                         np.asarray(events[0]["average_opening_m"])[::-1],
                         *extent_ham)
    past_edge = x[(x > se_2005) & (y_sept05 > 0.0)]
    display = {
        "note": "linear interpolation between column centres on the "
                "projected axis, carried flat from the end centres to "
                "the model edges (extent_km), for display only; it is "
                "not volume-preserving, so quoted volumes come from "
                "the native-grid tables",
        "fine_grid_step_km": DISPLAY_STEP_KM,
        "x_max_km": x_max,
        "extent_km": {"hamling_grid": list(extent_ham),
                      "epoch_2026": list(extent_2026)},
        "drawn_area_km3": drawn_area,
        "worst_drawn_vs_native": {"tag": worst[0],
                                  "relative_error": worst[1]},
        "sept05_drawn_at_plane_edge_m":
            float(_resample(np.array([se_2005]), s_nwse,
                            np.asarray(events[0]["average_opening_m"])
                            [::-1], *extent_ham)[0]),
        "sept05_drawn_positive_to_km":
            float(past_edge.max()) if past_edge.size else se_2005,
    }

    # The historical post-2005 files floor small opening while the
    # 2005 model and the 2026 fit keep all of it; state that asymmetry
    # with this product instead of only pointing at the three-epoch
    # report (the 0.5-m values are one uniform illustrative floor,
    # mirrored from that report, not a matched reconstruction of the
    # papers' local-uncertainty thresholds).
    floors = [float(delivery["openings_m"][tag]
                    [delivery["openings_m"][tag] > 0].min())
              for tag in cum_tags]
    threshold_sensitivity = {
        "post2005_positive_floors_m": [min(floors), max(floors)],
        "volume_2026_km3_at_0p5m_floor":
            epochs["epoch_2026"]["volume_at_0p5m_floor_km3"],
        "cumulative_2006_2010_km3_at_0p5m_floor":
            epochs["cumulative_2006_2010"]["volume_at_0p5m_floor_km3"],
        "note": "every canonical post-2005 file carries a positive "
                "floor, consistent with the papers' discard of opening "
                "below the local 2-sigma uncertainty, while the "
                "September 2005 model and the 2026 fit keep all fitted "
                "opening; magnitude and extent comparisons inherit "
                "that asymmetry",
    }

    suffix = "_gslc" if chain == "gslc" else ""
    payload = {
        "generated_by": "mh26.opening_profiles",
        "nisar_source": chain,
        "style": "after Fig. 13 of Hamling et al. (2009): top-10-km "
                 "depth-averaged opening per intrusion and the cumulative "
                 "stack, versus distance along the segment",
        "average_depth_km": AVERAGE_DEPTH_KM,
        "source": three_epoch["source"],
        "sept05_model_file": three_epoch["sept05_updated_model"]["file"],
        "sept05_model_sha256":
            three_epoch["sept05_updated_model"]["sha256"],
        "s_axis": {
            "frame": registration["frame"],
            "s_km": s_km.tolist(),
            "order": "as registered (the delivery's S->N column order, "
                     "so s decreases along the list)",
            "se_limit_2005_km": registration["se_limit_2005_km"],
            "plane_span_km": registration["plane_span_km"],
        },
        "events": events,
        "epoch_2026": epoch_2026,
        "display": display,
        "threshold_sensitivity": threshold_sensitivity,
        "cumulative": {
            "stack_order": [e["tag"] for e in events] + ["epoch_2026"],
            "volume_above_10km_km3_total": total,
        },
        "notes": [
            "the depth average is over 0-10 km everywhere, so a profile "
            "integrates exactly to the epoch's above-10-km volume "
            "(cross-checked against mh26.hamling on every run); the "
            "September 2005 profile averages Ian Hamling's 10-km model, "
            "the 2006-2010 profiles the upper 10 km of the delivered "
            "fields, and the 2026 profile this repository's fixed-plane "
            "fit, whose plane spans exactly 0-10 km",
            "the delivered inversions carry no uncertainties, so the "
            "panels have no error bars (the published figure's error "
            "bars came from the original inversions)",
            "panel labels are the intrusion months of the Hamling et "
            "al. (2009, 2010) sequence: the july09 file is the "
            "re-estimate of the June 2009 intrusion; exact onset dates "
            "for captions belong to the source catalogue",
            "input decoding, the canonical-13 selection, and the "
            "comparability caveats are those of "
            "results/uperp_three_epoch*.md, reused unchanged",
        ],
        "figures": [f"opening_profiles{suffix}.png",
                    f"opening_profiles_cumulative{suffix}.png"],
    }
    context = {**context, "three_epoch": three_epoch, "events": events,
               "fields": fields, "profile_2026": profile_2026,
               "s_2026": s_2026, "s_edges_2026": s_edges}
    return payload, context


def _threshold_caveat(payload: dict) -> str:
    floors = payload["threshold_sensitivity"]["post2005_positive_floors_m"]
    return (f"2006–2010 files floor small opening "
            f"({floors[0]:.1f}–{floors[1]:.1f} m);\n"
            "2005 and 2026 keep all fitted opening")


def _cumulative_panel(ax, payload, context, xlim, show_caveat=True):
    """Stack the per-event averages into the cumulative profile."""
    x = np.arange(xlim[0], xlim[1] + 0.01, DISPLAY_STEP_KM)
    centres_ham = np.asarray(payload["s_axis"]["s_km"])[::-1]  # NW->SE
    extent_ham = payload["display"]["extent_km"]["hamling_grid"]
    extent_2026 = payload["display"]["extent_km"]["epoch_2026"]
    bottom = np.zeros_like(x)
    for event in payload["events"]:
        values = np.asarray(event["average_opening_m"])[::-1]
        top = bottom + _resample(x, centres_ham, values, *extent_ham)
        ax.fill_between(x, bottom, top, color=event["color"], lw=0,
                        label=event["label"])
        # A hairline gap keeps neighbouring layers of similar hue apart.
        ax.plot(x, top, color="white", lw=0.35)
        bottom = top
    epoch_2026 = payload["epoch_2026"]
    chain_note = ("30-m GSLC" if payload["nisar_source"] == "gslc"
                  else "GUNW")
    top = bottom + _resample(
        x, np.asarray(epoch_2026["s_centres_km"]),
        np.asarray(epoch_2026["average_opening_m"]), *extent_2026)
    ax.fill_between(x, bottom, top, color=epoch_2026["color"], lw=0,
                    label=f"{epoch_2026['label']} ({chain_note} chain)")
    ax.plot(x, top, color="0.15", lw=0.7)
    # Beyond the plane's SE edge the stack holds historical models
    # only; grey the region as on the July 2026 panel so it does not
    # read as zero 2026 opening.
    plane_end = epoch_2026["s_edges_km"][-1]
    if xlim[1] > plane_end:
        ax.axvspan(plane_end, xlim[1], color="0.88", zorder=0)
        ax.text(plane_end + 0.4, 0.30,
                f"no 2026 model SE of {plane_end:.0f} km",
                transform=ax.get_xaxis_transform(), rotation=90,
                ha="left", va="center", fontsize=7.5, color="0.35",
                zorder=5)
    ax.set_xlim(*xlim)
    ax.set_ylim(0.0, 11.5)
    ax.set_title("Cumulative opening, September 2005 – July 2026",
                 loc="left", fontsize=10.5, fontweight="bold", pad=5)
    ax.set_title("after Fig. 13 of Hamling et al. (2009)", loc="right",
                 fontsize=8, color="0.30", pad=6)
    se_2005 = payload["s_axis"]["se_limit_2005_km"]
    ax.axvline(se_2005, color="0.25", lw=0.9, ls=(0, (4, 3)), zorder=4)
    ax.text(se_2005 + 0.4, 0.30, "Wright 2006 plane edge",
            transform=ax.get_xaxis_transform(), rotation=90, ha="left",
            va="center", fontsize=7.5, color="0.20",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 1.2}, zorder=5)
    total = payload["cumulative"]["volume_above_10km_km3_total"]
    n_events = len(payload["events"]) + 1
    label = f"{total:.2f} km$^3$ above 10 km · {n_events} intrusions"
    if show_caveat:
        label += "\n" + _threshold_caveat(payload).replace("\n", " ")
    ax.text(0.24, 0.93, label,
            transform=ax.transAxes, ha="left", va="top", fontsize=8.5,
            linespacing=1.5,
            bbox={"facecolor": "white", "edgecolor": "none",
                  "alpha": 0.8, "pad": 1.6})
    ax.set_xlabel("Distance along the 2005 dike plane, NW→SE (km)")


def profiles_figure(payload: dict, context: dict, chain: str) -> Path:
    xlim = (0.0, payload["display"]["x_max_km"])
    small_max = max(
        max(e["max_average_m"] for e in payload["events"][1:]),
        payload["epoch_2026"]["max_average_m"])
    small_top = 0.25 * np.ceil(small_max * 1.1 / 0.25)

    n_small = len(payload["events"])  # 13 small events + the 2026 panel
    heights = [1.7] + [1.0] * n_small + [3.4]
    fig, axes = plt.subplots(len(heights), 1,
                             figsize=(9.0, 0.62 * sum(heights) + 1.2),
                             sharex=True, constrained_layout=True,
                             gridspec_kw={"height_ratios": heights})

    se_2005 = payload["s_axis"]["se_limit_2005_km"]
    centres_ham = np.asarray(payload["s_axis"]["s_km"])[::-1]  # NW->SE
    extent_ham = payload["display"]["extent_km"]["hamling_grid"]
    extent_2026 = payload["display"]["extent_km"]["epoch_2026"]

    def draw(ax, s, profile, event, top):
        ax.fill_between(s, 0.0, profile, color=event["color"], alpha=0.30,
                        lw=0)
        ax.plot(s, profile, color=event["color"], lw=1.4)
        ax.set_ylim(0.0, top)
        ax.set_xlim(*xlim)
        ax.text(0.988, 0.86, event["label"], transform=ax.transAxes,
                ha="right", va="top", fontsize=8.6, color="0.10",
                fontweight="bold")
        ax.tick_params(labelsize=8)

    sept05 = payload["events"][0]
    draw(axes[0], *_extend_to_edges(
             centres_ham, np.asarray(sept05["average_opening_m"])[::-1],
             *extent_ham), sept05,
         1.12 * sept05["max_average_m"])
    axes[0].axvline(se_2005, color="0.25", lw=0.8, ls=(0, (4, 3)))
    axes[0].set_title("Depth-averaged opening, 0–10 km", loc="left",
                      fontsize=10.5, fontweight="bold", pad=5)
    axes[0].set_title("Hamling et al. (2009, 2010) models · "
                      "no uncertainties delivered", loc="right",
                      fontsize=8, color="0.30", pad=6)

    for ax, event in zip(axes[1:], payload["events"][1:]):
        draw(ax, *_extend_to_edges(
                 centres_ham,
                 np.asarray(event["average_opening_m"])[::-1],
                 *extent_ham), event, small_top)
        ax.set_yticks([0.0, small_top])

    ax_2026 = axes[1 + len(payload["events"][1:])]
    epoch_2026 = payload["epoch_2026"]
    draw(ax_2026, *_extend_to_edges(
             np.asarray(epoch_2026["s_centres_km"]),
             np.asarray(epoch_2026["average_opening_m"]),
             *extent_2026), epoch_2026, small_top)
    ax_2026.set_yticks([0.0, small_top])
    plane_end = epoch_2026["s_edges_km"][-1]
    if xlim[1] > plane_end:
        # Beyond the plane's SE edge there is no 2026 model; grey it so
        # it does not read as zero opening.
        ax_2026.axvspan(plane_end, xlim[1], color="0.88", zorder=0)
    ax_2026.axvline(se_2005, color="0.25", lw=0.8, ls=(0, (4, 3)))
    source_note = ("this study, 30-m GSLC chain" if chain == "gslc"
                   else "this study, GUNW chain")
    ax_2026.set_title(source_note, loc="right", fontsize=8, color="0.30",
                      pad=3)

    _cumulative_panel(axes[-1], payload, context, xlim)
    fig.supylabel("Average opening (m)", fontsize=10)

    suffix = "_gslc" if chain == "gslc" else ""
    out = RESULTS / f"opening_profiles{suffix}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=250)
    plt.close(fig)
    return Path(f"{out}.png")


def cumulative_figure(payload: dict, context: dict, chain: str) -> Path:
    """The cumulative stack alone, sized to sit under the three-epoch
    comparison figure (same width and x axis)."""
    xlim = (0.0, payload["display"]["x_max_km"])
    fig, ax = plt.subplots(figsize=(9.0, 3.4), constrained_layout=True)
    _cumulative_panel(ax, payload, context, xlim, show_caveat=False)
    ax.set_ylabel("Average opening (m)")
    # Identity needs a legend here: the labelled per-event panels of the
    # full profiles figure are absent. The threshold caveat rides as
    # its title so the panel carries the asymmetry statement on its own.
    legend = ax.legend(loc="upper right", fontsize=6.8, ncols=3,
                       frameon=True, framealpha=0.9, edgecolor="none",
                       columnspacing=1.0, handlelength=1.2,
                       handleheight=0.9, borderaxespad=0.2,
                       title=_threshold_caveat(payload),
                       title_fontsize=6.8)
    legend.get_title().set_color("0.30")
    suffix = "_gslc" if chain == "gslc" else ""
    out = RESULTS / f"opening_profiles_cumulative{suffix}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=250)
    plt.close(fig)
    return Path(f"{out}.png")


def run(chain: str, fit_path: Path | None, data_dir: Path) -> None:
    fit_path = Path(fit_path or _find_fit(chain)).resolve()
    payload, context = compute_report(chain, fit_path, data_dir)
    figure = profiles_figure(payload, context, chain)
    cumulative = cumulative_figure(payload, context, chain)

    suffix = "_gslc" if chain == "gslc" else ""
    json_path = RESULTS / f"opening_profiles{suffix}.json"
    json_path.write_text(json.dumps(payload, indent=1, sort_keys=True))

    lines = [
        "# Depth-averaged opening per intrusion "
        "(after Hamling et al. 2009, Fig. 13)",
        "",
        f"NISAR chain of the 2026 fit: {chain} "
        f"({payload['epoch_2026']['reference_fit']}).",
        "",
        "Top-10-km depth-averaged opening of every intrusion of the "
        "2005–2010 Dabbahu–Manda Hararo sequence and of the July 2026 "
        "dike, versus distance along the Wright et al. (2006) plane "
        "(NW→SE s axis of the 2026 contract), in the style of Fig. 13 "
        "of Hamling et al. (2009). The 2005 profile averages Ian "
        "Hamling's 10-km September 2005 model, the 2006–2010 profiles "
        "the upper 10 km of his canonical thirteen files, and the 2026 "
        "profile this repository's fixed-plane fit (its plane spans "
        "exactly 0–10 km). On the native grid a profile integrates "
        "exactly to the epoch's above-10-km volume: average × 10 km × "
        "column width = column volume (cross-checked against "
        "`uperp_three_epoch*.json` on every run). The drawn curves "
        "interpolate linearly between column centres on the projected "
        "axis and are not volume-preserving: the largest drawn-area "
        "deviation is "
        f"{EVENT_LABELS.get(payload['display']['worst_drawn_vs_native']['tag'], LABEL_2026)}"
        " at "
        f"{payload['display']['worst_drawn_vs_native']['relative_error']:.1%}"
        " and the cumulative stack draws as "
        f"{payload['display']['drawn_area_km3']['cumulative']:.2f} km3 "
        "against the native "
        f"{payload['cumulative']['volume_above_10km_km3_total']:.2f} km3 "
        "(every drawn area is in "
        "the JSON's display block), so quoted volumes always come from "
        "the native-grid tables. The delivered inversions carry no "
        "uncertainties, so unlike the published figure the panels have "
        "no error bars. Magnitude and extent comparisons also inherit "
        "a threshold asymmetry, stated on the cumulative panels: every "
        "canonical post-2005 file floors small opening ("
        f"{payload['threshold_sensitivity']['post2005_positive_floors_m'][0]:.2f}"
        "–"
        f"{payload['threshold_sensitivity']['post2005_positive_floors_m'][1]:.2f}"
        " m, consistent with the papers' discard of opening below the "
        "local 2σ uncertainty), while the September 2005 model and the "
        "2026 fit keep all fitted opening. As one uniform illustrative "
        "floor — a sensitivity, not a reconstruction of the historical "
        "thresholds — 0.5 m reduces the July 2026 volume from "
        f"{payload['epoch_2026']['volume_above_10km_km3']:.3f} to "
        f"{payload['threshold_sensitivity']['volume_2026_km3_at_0p5m_floor']:.3f}"
        " km3 (about "
        f"{1 - payload['threshold_sensitivity']['volume_2026_km3_at_0p5m_floor'] / payload['epoch_2026']['volume_above_10km_km3']:.0%}"
        "); the matching cumulative sensitivity is in "
        "`uperp_three_epoch*.md`, whose remaining comparability caveats "
        "(different regularization and greens functions) apply "
        "unchanged.",
        "",
        "| intrusion | file | volume above 10 km (km3) | "
        "max average opening (m) | at s (km) |",
        "|---|---|---|---|---|",
    ]
    for event in payload["events"]:
        lines.append(
            f"| {event['label']} | {event['file']} | "
            f"{event['volume_above_10km_km3']:.3f} | "
            f"{event['max_average_m']:.2f} | {event['s_at_max_km']:.1f} |")
    epoch_2026 = payload["epoch_2026"]
    lines.append(
        f"| {epoch_2026['label']} ({chain}) | "
        f"{Path(epoch_2026['reference_fit']).name} | "
        f"{epoch_2026['volume_above_10km_km3']:.3f} | "
        f"{epoch_2026['max_average_m']:.2f} | "
        f"{epoch_2026['s_at_max_km']:.1f} |")
    total = payload["cumulative"]["volume_above_10km_km3_total"]
    lines += [
        "",
        f"The cumulative stack (bottom panel of `{figure.name}`, also "
        f"standalone as `{cumulative.name}`) sums the fifteen profiles "
        f"in event order: {total:.2f} km3 above 10 km in total. Panel "
        "labels are the intrusion months of the Hamling et al. (2009, "
        "2010) sequence — the july09 file is the re-estimate of the "
        "June 2009 intrusion — and exact onset dates for captions "
        "belong to the source catalogue. Beyond the SE end of the 2026 "
        "model plane "
        f"({epoch_2026['s_edges_km'][-1]:.1f} km) there is no 2026 "
        "model, and that region is shaded grey in the July 2026 panel "
        "and in both cumulative panels (the stack there holds "
        "historical models only). The dashed line is the Wright et "
        "al. (2006) model-plane SE edge: no 2005 model column "
        "southeast of it holds opening, and the drawn 2005 taper "
        "crosses it only by the interpolation span between column "
        "centres ("
        f"{payload['display']['sept05_drawn_at_plane_edge_m']:.2f} m "
        "at the line, positive to "
        f"{payload['display']['sept05_drawn_positive_to_km']:.1f} km).",
        "",
    ]
    md_path = RESULTS / f"opening_profiles{suffix}.md"
    md_path.write_text("\n".join(lines))

    print(f"wrote {figure}, {cumulative.name}, {json_path.name}, "
          f"{md_path.name}")
    print(f"events: {len(payload['events'])} + 2026; cumulative "
          f"{total:.2f} km3 above 10 km; 2026 max average "
          f"{epoch_2026['max_average_m']:.2f} m at "
          f"s={epoch_2026['s_at_max_km']:.1f} km")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gslc", action="store_true",
                        help="use the 30-m GSLC-chain fixed-plane fit")
    parser.add_argument("--fit", type=Path,
                        help="explicit Wright-plane fit result JSON")
    parser.add_argument("--dir", type=Path, default=DATA_DIR,
                        help="directory with the Hamling opening files")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    run("gslc" if args.gslc else "gunw", args.fit, args.dir)
