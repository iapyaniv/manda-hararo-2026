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

Figures:  python -m mh26.opening_history          (GUNW chain)
          python -m mh26.opening_history --gslc   (30-m GSLC chain)
Outputs: results/opening_history[_gslc].{png,pdf} and
results/opening_history_common_scale[_gslc].{png,pdf}.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .hamling import DATA_DIR, _find_fit, draw_epoch_sections
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


def run(chain: str, fit_path: Path | None, data_dir: Path) -> None:
    fit_path = Path(fit_path or _find_fit(chain)).resolve()
    payload, context = compute_report(chain, fit_path, data_dir)
    for common_scale in (False, True):
        figure = history_figure(payload, context, chain, common_scale)
        print(f"wrote {figure}")


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
