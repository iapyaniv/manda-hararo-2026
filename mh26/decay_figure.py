"""Supplement figure for the post-event decay fit.

Pure layout of the pinned Sentinel-1 projection report
(results/s1_tracks.json) and the pinned NISAR continuation product
(results/continuation.json): each window's projection onto the
co-event dike pattern drawn over its date span with its 95%
cluster-bootstrap interval, the nine pre-event null windows, the
fitted window-integrated values from the exponential decay, and the
independent NISAR window that the fit does not use. Windows are
selected by their recorded role, never by value. No quantity is
recomputed here except the fitted model evaluated per window.

Run after mh26.s1_tracks and mh26.continuation:
    python -m mh26.decay_figure
Writes results/decay_fit_windows.{png,pdf}.
"""

from __future__ import annotations

import datetime as dt
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from .params import RESULTS

T0 = dt.date(2026, 7, 18)
NISAR_WINDOW = "2026-07-19→2026-08-12"
COLORS = {"ascending": "#2166ac", "descending": "#d95f02",
          "pre": "0.45", "nisar": "#1b7837", "fit": "k"}


def window_dates(row: dict) -> tuple[dt.date, dt.date]:
    ref, sec = row["pair"].split("_")
    to_date = lambda s: dt.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
    return to_date(ref), to_date(sec)


def draw_window(ax, start, end, value, interval, color, lw=2.6, label=None):
    ax.plot([start, end], [value, value], color=color, lw=lw,
            solid_capstyle="butt", zorder=3, label=label)
    if interval is not None:
        mid = start + (end - start) / 2
        ax.plot([mid, mid], list(interval), color=color, lw=0.9, zorder=2)


def main() -> None:
    report = json.loads((RESULTS / "s1_tracks.json").read_text())
    continuation = json.loads((RESULTS / "continuation.json").read_text())
    if continuation["window"] != NISAR_WINDOW:
        raise SystemExit(f"continuation window {continuation['window']} "
                         f"is not the expected {NISAR_WINDOW}")
    decay = report["post_event_decay"]
    null = report["null_summary"]
    amplitude = decay["amplitude_total_pattern_fraction"]
    tau = decay["tau_days"]

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    band = null["max_abs"]
    ax.axhspan(-band, band, color="0.88", zorder=0,
               label=f"pre-event null range (|max| {band:.3f}, "
                     f"σ {null['standard_deviation']:.3f})")
    ax.axhline(0, color="0.6", lw=0.6, zorder=1)
    ax.axvline(T0, color="0.4", lw=0.8, ls=":", zorder=1)
    ax.annotate("18 July", (T0, 0.185), fontsize=8, color="0.3",
                ha="right", va="top", rotation=90, xytext=(-3, 0),
                textcoords="offset points")

    seen = set()
    for row in report["projections"]:
        if row["role"] not in ("pre-event", "post-event"):
            continue  # the two co-event windows are off this scale
        start, end = window_dates(row)
        interval = row.get("pattern_fraction_cluster_wild_95")
        if row["role"] == "pre-event":
            key, color, lw = "pre", COLORS["pre"], 2.2
            label = "pre-event window"
        else:
            key = row["geometry"]
            color, lw = COLORS[key], 2.6
            label = f"Sentinel-1 {key} window"
        draw_window(ax, start, end, row["pattern_fraction"], interval,
                    color, lw=lw, label=None if key in seen else label)
        seen.add(key)
        if row["role"] == "post-event":
            t1, t2 = (start - T0).days, (end - T0).days
            fitted = amplitude * (np.exp(-t1 / tau) - np.exp(-t2 / tau))
            # The fitted value spans the same window as its observation,
            # so the solid/dashed gap is the fit residual.
            ax.plot([start, end], [fitted, fitted], color=COLORS["fit"],
                    lw=1.1, ls=(0, (4, 3)), solid_capstyle="butt", zorder=4,
                    label=None if "fit" in seen else
                    "decay fit, window-integrated")
            seen.add("fit")

    n_start = dt.date(2026, 7, 19)
    n_end = dt.date(2026, 8, 12)
    draw_window(ax, n_start, n_end, continuation["pattern_fraction"],
                continuation["pattern_fraction_cluster_wild_95"],
                COLORS["nisar"], label="NISAR window (not used in the fit)")
    ax.plot([n_start, n_end],
            [report["nisar_cross_sensor_check"]["predicted_from_s1_decay"]] * 2,
            color=COLORS["fit"], lw=1.1, ls=(0, (4, 3)),
            solid_capstyle="butt", zorder=4)

    ax.set_ylabel("projection on the co-event dike pattern\n"
                  "(fraction of co-event signal)")
    ax.set_title(
        "Post-event motion with the co-event dike pattern and the "
        f"exponential decay fit\n(total {amplitude:.3f} of the co-event "
        f"signal, decay time {tau:.1f} days; a descriptive fit)",
        fontsize=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.set_xlim(dt.date(2026, 6, 3), dt.date(2026, 8, 20))
    ax.set_ylim(-0.045, 0.19)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = RESULTS / f"decay_fit_windows.{ext}"
        fig.savefig(out, dpi=110)
        print(f"figure: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
