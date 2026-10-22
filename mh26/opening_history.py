"""Single-figure opening history: the three-epoch U_perp sections
over the cumulative depth-averaged profile.

Combines the three cross-sections of mh26.hamling (September 2005,
cumulative 2006-2010, July 2026) and the cumulative stack of
mh26.opening_profiles into one manuscript-figure candidate on the
shared along-plane axis. Every drawn quantity, panel note and caveat
comes from those two reports (mh26.validate rebuilds both payloads);
this module only lays their panels out together and computes nothing
new.

The common-scale variant draws the three sections on one colour scale,
topped at the September 2005 panel's 16.5 m, so the amplitudes compare
directly; the July 2026 section then spans only the lowest colour
steps and reads mainly through its contour lines, which is why the
default variant keeps the per-section scales.

The manuscript variant (Figure 4 of the GRL letter) is the common-scale
figure laid out two columns wide (190 x 118 mm): the four panels at the
left, the colour bar and the event legend in a strip at the right, and
none of the in-panel notes, volume labels or source lines. Those go into
the caption, which `manuscript_caption` writes beside the figure from the
same report numbers (volumes above 10 km, maxima, contour levels, the
plane edge, the extent of the 2026 model, the opening floors).

Figures:  python -m mh26.opening_history          (GUNW chain)
          python -m mh26.opening_history --gslc   (30-m GSLC chain)
Outputs: results/opening_history[_gslc].{png,pdf},
results/opening_history_common_scale[_gslc].{png,pdf},
results/opening_history_manuscript[_gslc].{png,pdf} and
results/opening_history_manuscript[_gslc]_caption.txt.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .hamling import DATA_DIR, SECTION_CONTOURS, _find_fit, draw_epoch_sections
from .opening_profiles import (_cumulative_legend, _cumulative_panel,
                               compute_report)
from .params import RESULTS

# The September 2005 panel's scale in mh26.hamling; the shared scale
# must cover the largest drawn opening.
COMMON_VMAX = 16.5


def history_figure(payload: dict, context: dict, chain: str,
                   common_scale: bool) -> Path:
    fig, axes = plt.subplots(4, 1, figsize=(9.0, 11.3), sharex=True,
                             constrained_layout=True,
                             gridspec_kw={"height_ratios":
                                          [1.0, 1.0, 1.0, 1.25]})
    xlim = draw_epoch_sections(
        fig, axes[:3], context["registration"], context["fit"],
        context["geometry"], chain, context["epochs"], context["cum"],
        context["updated_field"],
        common_vmax=COMMON_VMAX if common_scale else None)
    if xlim[1] != payload["display"]["x_max_km"]:
        raise ValueError(
            f"section xlim {xlim[1]} != profile x_max_km "
            f"{payload['display']['x_max_km']}; the two modules no "
            "longer derive the same axis")

    ax_cum = axes[3]
    # The legend and its caveat title fill the upper right here (the
    # colorbar strip narrows the axis), so the total goes upper left,
    # above the September 2005 peak.
    _cumulative_panel(ax_cum, payload, context, xlim, show_caveat=False,
                      label_xy=(0.012, 0.975))
    ax_cum.set_ylabel("Average opening (m)")
    _cumulative_legend(ax_cum, payload)
    # Reserve the same right-hand strip the section rows give to their
    # colorbars, or the shared x axis drifts out of column alignment.
    spacer = fig.colorbar(
        matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(0.0, 1.0)),
        ax=ax_cum, shrink=0.94, pad=0.015)
    spacer.ax.set_visible(False)

    suffix = "_gslc" if chain == "gslc" else ""
    name = ("opening_history_common_scale" if common_scale
            else "opening_history")
    out = RESULTS / f"{name}{suffix}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=250)
    plt.close(fig)
    return Path(f"{out}.png")


MANUSCRIPT_MM = {"width": 190.0, "height": 118.0, "left": 12.0, "panel_w": 139.0,
                 "section_h": 19.5, "cumulative_h": 31.0, "gap": 5.0, "bottom": 9.0,
                 "strip_x": 160.0, "bar_w": 3.6}


def _levels(levels) -> str:
    text = [f"{v:g}" for v in levels]
    return ", ".join(text[:-1]) + " and " + text[-1]


def manuscript_caption(payload: dict, context: dict) -> str:
    """The Figure 4 caption, with every number taken from the report."""
    epochs = context["epochs"]
    sept05 = epochs["sept05_updated"]
    cum = epochs["cumulative_2006_2010"]
    e2026 = epochs["epoch_2026"]
    floors = payload["threshold_sensitivity"]["post2005_positive_floors_m"]
    floor_2005 = payload["threshold_sensitivity"]["sept05_model_min_positive_m"]
    n_post = len(payload["events"]) - 1
    names = {13: "thirteen"}
    se_edge = context["registration"]["se_limit_2005_km"]
    plane_end = payload["epoch_2026"]["s_edges_km"][-1]
    max_offplane_km = context["registration"]["max_offplane_m"] / 1000
    return (
        "Figure 4. Opening on the Dabbahu–Manda Hararo dike plane, 2005 to 2026. "
        "(a–c) Opening of the September 2005 dike, in the revised model supplied by "
        "I. Hamling in 2026 that updates Hamling et al. (2009); of the "
        f"{names.get(n_post, n_post)} dikes of 2006–2010 summed, in the models supplied "
        "by I. Hamling (Hamling et al., 2009, 2010, and the May 2010 dike); and of the "
        "July 2026 dike (this study); viewed from the south-west, "
        "upper 10 km, on one colour scale. The horizontal axis is distance along the "
        "2005 dike plane of Wright et al. (2006); the 2005–2010 models stay on the grid "
        f"of Hamling et al. (2009) and are projected onto that axis (within "
        f"{max_offplane_km:.1f} km of the plane). Contours are at "
        f"{_levels(SECTION_CONTOURS[0])} m in (a), {_levels(SECTION_CONTOURS[1])} m "
        f"in (b) and {_levels(SECTION_CONTOURS[2])} m in (c). The dashed line is the "
        f"south-eastern end of the 2005 model plane ({se_edge:.0f} km); the grey band "
        f"lies beyond the 2026 model plane ({plane_end:.0f} km). Volumes above 10 km "
        f"depth are {sept05['volume_above_10km_km3']:.1f} km³ (2005), "
        f"{cum['volume_above_10km_km3']:.1f} km³ (2006–2010) and "
        f"{e2026['volume_above_10km_km3']:.2f} km³ (2026), with maximum openings of "
        f"{sept05['max_opening_m']:.0f}, {cum['max_opening_m']:.0f} and "
        f"{e2026['max_opening_m']:.1f} m. (d) Opening averaged over the upper 10 km, "
        "stacked by intrusion in the order of the legend, in the layout of Fig. 13 of "
        "Hamling et al. (2009); the curves are interpolated for display, and the volumes "
        "quoted are from the models' own grids. The 2006–2010 models omit opening below "
        f"{floors[0]:.1f}–{floors[1]:.1f} m; the 2005 model, as supplied, has no positive "
        f"opening below {floor_2005:.1f} m; the 2026 model keeps all fitted opening.\n")


def manuscript_figure(payload: dict, context: dict, chain: str) -> Path:
    mm = MANUSCRIPT_MM
    w, h = mm["width"], mm["height"]

    def axes_mm(x, y, aw, ah):
        return fig.add_axes([x / w, y / h, aw / w, ah / h])

    with plt.rc_context({"font.size": 7, "axes.labelsize": 7.5, "xtick.labelsize": 7,
                         "ytick.labelsize": 7, "axes.linewidth": 0.6,
                         "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                         "xtick.major.size": 2.2, "ytick.major.size": 2.2,
                         "pdf.fonttype": 42}):
        fig = plt.figure(figsize=(w / 25.4, h / 25.4))
        y = mm["bottom"]
        ax_cum = axes_mm(mm["left"], y, mm["panel_w"], mm["cumulative_h"])
        y += mm["cumulative_h"] + mm["gap"]
        sections = []
        for _ in range(3):
            sections.append(axes_mm(mm["left"], y, mm["panel_w"], mm["section_h"]))
            y += mm["section_h"] + mm["gap"]
        sections = sections[::-1]            # top row first: 2005, 2006-2010, 2026
        y_bar0 = sections[1].get_position().y0 * h
        y_bar1 = sections[0].get_position().y1 * h
        cax = axes_mm(mm["strip_x"], y_bar0, mm["bar_w"], y_bar1 - y_bar0)
        xlim = draw_epoch_sections(
            fig, sections, context["registration"], context["fit"],
            context["geometry"], chain, context["epochs"], context["cum"],
            context["updated_field"], common_vmax=COMMON_VMAX, manuscript=True,
            cax=cax)
        for ax in sections:
            ax.tick_params(labelbottom=False)
        _cumulative_panel(ax_cum, payload, context, xlim, manuscript=True)
        ax_cum.set_ylabel("Average opening (m)")
        handles, labels = ax_cum.get_legend_handles_labels()
        labels = [label.replace(" (30-m GSLC chain)", "").replace(" (GUNW chain)", "")
                  for label in labels]
        y_leg_top = sections[2].get_position().y1
        fig.legend(handles, labels, loc="upper left",
                   bbox_to_anchor=(mm["strip_x"] / w - 0.012, y_leg_top + 0.004),
                   fontsize=6.6, frameon=False, handlelength=1.3, handleheight=0.9,
                   labelspacing=0.42, borderaxespad=0.0)
        suffix = "_gslc" if chain == "gslc" else ""
        out = RESULTS / f"opening_history_manuscript{suffix}"
        for ext in ("png", "pdf"):
            fig.savefig(f"{out}.{ext}", dpi=400)
        plt.close(fig)
    Path(f"{out}_caption.txt").write_text(manuscript_caption(payload, context))
    return Path(f"{out}.png")


def run(chain: str, fit_path: Path | None, data_dir: Path) -> None:
    fit_path = Path(fit_path or _find_fit(chain)).resolve()
    payload, context = compute_report(chain, fit_path, data_dir)
    for common_scale in (False, True):
        figure = history_figure(payload, context, chain, common_scale)
        print(f"wrote {figure}")
    print(f"wrote {manuscript_figure(payload, context, chain)}")


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
