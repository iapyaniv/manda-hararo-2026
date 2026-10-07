"""Supplement figure for the post-event decay fit.

Pure layout of the pinned Sentinel-1 projection report
(results/s1_tracks.json): each window's projection onto the co-event
dike pattern drawn over its date span, the nine pre-event null windows,
the fitted window-integrated values from the exponential decay, and the
independent NISAR window that the fit does not use. No quantity is
recomputed here except the fitted model evaluated per window.

Run after mh26.s1_tracks:  python -m mh26.decay_figure
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
COEVENT_FRACTION_MIN = 0.5


def window_dates(row: dict) -> tuple[dt.date, dt.date]:
    ref, sec = row["pair"].split("_")
    to_date = lambda s: dt.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
    return to_date(ref), to_date(sec)


def main() -> None:
    report = json.loads((RESULTS / "s1_tracks.json").read_text())
    decay = report["post_event_decay"]
    null = report["null_summary"]
    nisar = report["nisar_cross_sensor_check"]
    amplitude = decay["amplitude_total_pattern_fraction"]
    tau = decay["tau_days"]

    pre, post = [], []
    for row in report["projections"]:
        start, end = window_dates(row)
        if row["pattern_fraction"] > COEVENT_FRACTION_MIN:
            continue  # the two co-event windows are off this scale
        (post if end > T0 else pre).append((start, end,
                                            row["pattern_fraction"]))

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    band = null["max_abs"]
    ax.axhspan(-band, band, color="0.88", zorder=0,
               label=f"pre-event null range (|max| {band:.3f}, "
                     f"σ {null['standard_deviation']:.3f})")
    ax.axhline(0, color="0.6", lw=0.6, zorder=1)
    ax.axvline(T0, color="0.4", lw=0.8, ls=":", zorder=1)
    ax.annotate("18 July", (T0, 0.158), fontsize=8, color="0.3",
                ha="right", rotation=90, xytext=(-3, 0),
                textcoords="offset points")

    for start, end, value in pre:
        ax.plot([start, end], [value, value], color="0.45", lw=2.2,
                solid_capstyle="butt", zorder=2)
    for k, (start, end, value) in enumerate(post):
        ax.plot([start, end], [value, value], color="#2166ac", lw=2.6,
                solid_capstyle="butt", zorder=3,
                label="Sentinel-1 post-event window" if k == 0 else None)
        t1, t2 = (start - T0).days, (end - T0).days
        fitted = amplitude * (np.exp(-t1 / tau) - np.exp(-t2 / tau))
        mid = start + (end - start) / 2
        ax.plot([mid], [fitted], marker="o", mfc="none", mec="#2166ac",
                ms=6, zorder=4,
                label="decay fit, window-integrated" if k == 0 else None)

    n_start = dt.date(2026, 7, 19)
    n_end = dt.date(2026, 8, 12)
    ax.plot([n_start, n_end], [nisar["observed_nisar_pattern_fraction"]] * 2,
            color="#1b7837", lw=2.6, solid_capstyle="butt", zorder=3,
            label="NISAR window (not used in the fit)")
    mid = n_start + (n_end - n_start) / 2
    ax.plot([mid], [nisar["predicted_from_s1_decay"]], marker="o",
            mfc="none", mec="#1b7837", ms=6, zorder=4)

    ax.set_ylabel("projection on the co-event dike pattern\n"
                  "(fraction of co-event signal)")
    ax.set_title(
        "Post-event motion with the co-event dike pattern and the "
        f"exponential decay fit\n(total {amplitude:.3f} of the co-event "
        f"signal, decay time {tau:.1f} days; a descriptive fit)",
        fontsize=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.set_xlim(dt.date(2026, 6, 3), dt.date(2026, 8, 20))
    ax.set_ylim(-0.045, 0.175)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = RESULTS / f"decay_fit_windows.{ext}"
        fig.savefig(out, dpi=110)
        print(f"figure: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
