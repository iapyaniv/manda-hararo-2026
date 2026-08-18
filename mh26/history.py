"""Step 4 — pre-event probe history across all extracted pairs.

Reports the 1-km-radius-disk mean of ionosphere-corrected, far-field-referenced
LOS at three probe points (the two co-event flank peaks and a southern
candidate-precursor patch) for every .npz in data/nisar, grouped by track.
The sixteen pre-event windows (November 2025 – mid-July 2026) are the
"no deformation above the noise" null result.

Run from the repo root:  python -m mh26.history

Writes results/probe_history.md.
"""

from __future__ import annotations

import json
import re

import numpy as np

from .gunw import corrected_los, load_pair, pair_window, reference_farfield
from .params import DATA_NISAR, PROBES, RESULTS

PROBE_RADIUS = 1000.0  # m


def probe_mean(los, x, y, px, py, r=PROBE_RADIUS):
    """Mean LOS (cm) in a disk of radius r around (px, py); NaN if under 6 pixels."""
    X, Y = np.meshgrid(x, y)
    m = (np.hypot(X - px, Y - py) < r) & np.isfinite(los)
    return float(np.nanmean(los[m])) * 100 if m.sum() > 5 else np.nan


def track_of(tag: str) -> str:
    """Group tag: strip a trailing _YYYYMMDD reference date if present."""
    return tag.rsplit("_", 1)[0] if re.search(r"_\d{8}$", tag) else tag


def main() -> None:
    rows = []
    for npz in sorted(DATA_NISAR.glob("*.npz")):
        tag = npz.stem
        d = load_pair(tag)
        los, _ = corrected_los(d)
        ref = reference_farfield(los, d["x"], d["y"])
        span = pair_window(d) or tag
        vals = {p: probe_mean(ref, d["x"], d["y"], *xy) for p, xy in PROBES.items()}
        start, end = span.split("→") if "→" in span else ("", "")
        if end and end <= "2026-07-15":
            phase = "pre-event"
        elif start and start >= "2026-07-19":
            phase = "post-event"
        else:
            phase = "co-event"
        rows.append((track_of(tag), tag, span, phase, vals))

    rows.sort(key=lambda r: (r[0], r[2]))
    pre_event = [abs(value) for _, _, _, phase, values in rows if phase == "pre-event"
                 for value in values.values() if np.isfinite(value)]
    pre_count = sum(phase == "pre-event" for _, _, _, phase, _ in rows)
    lines = [
        "# Probe history (cm LOS, 1-km-radius mean, +=toward satellite)",
        "",
        f"- pre-event windows: {pre_count}; maximum absolute probe mean: "
        f"{max(pre_event):.1f} cm",
        "- screening: ionosphere correction, far-field plane, coherence, connected",
        "  component, and GUNW valid-sample mask",
        "",
        "| track/pair | window | phase | EAST flank | WEST flank | SOUTH patch |",
        "|---|---|---|---:|---:|---:|",
    ]
    for track, _, span, phase, v in rows:
        lines.append(
            f"| {track} | {span} | {phase} | {v['EAST']:+.1f} | "
            f"{v['WEST']:+.1f} | {v['SOUTH']:+.1f} |"
        )
    out = "\n".join(lines)
    print(out)
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "probe_history.md").write_text(out + "\n")
    (RESULTS / "probe_history.json").write_text(json.dumps({
        "probe_radius_m": PROBE_RADIUS,
        "pre_event_windows": pre_count,
        "pre_event_max_abs_probe_cm": max(pre_event),
        "screening": (
            "ionosphere + far-field plane + coherence + connected component "
            "+ GUNW valid-sample mask"
        ),
        "rows": [
            {
                "track": track,
                "tag": tag,
                "window": span,
                "phase": phase,
                "probe_means_cm": {
                    key: float(value) for key, value in values.items()
                },
            }
            for track, tag, span, phase, values in rows
        ],
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
