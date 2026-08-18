"""Processing-choice sensitivity for the uniform-opening dike inversion.

This is not a probabilistic uncertainty analysis. It repeats the inversion
with different block sizes, far-field reference radii, LOS interpolation
heights, geometry weights, and fit-domain radii. The output reports the range
of fitted parameters across these processing choices.

Run from the repository root after ``python -m mh26.invert``::

    python -m mh26.sensitivity

Writes ``results/model_sensitivity.md`` and ``.json``.
"""

from __future__ import annotations

import json

import numpy as np

from .invert import fit, load_observations
from .params import RESULTS


VARIANTS = (
    ("default", {}, False),
    ("0.96 km blocks", {"block": 12}, False),
    ("1.92 km blocks", {"block": 24}, False),
    ("25 km far field", {"farfield_radius": 25_000.0}, False),
    ("32 km far field", {"farfield_radius": 32_000.0}, False),
    ("LOS height 0 m", {"los_height": 0.0}, False),
    ("LOS height 1000 m", {"los_height": 1_000.0}, False),
    ("equal geometry weight", {}, True),
    ("fit within 30 km", {"max_radius": 30_000.0}, False),
)


def summarize(label: str, p: np.ndarray, rms: float, n_obs: int) -> dict:
    volume = p[3] * (p[5] - p[4]) * p[6] / 1e9
    return {
        "variant": label,
        "n_obs": int(n_obs),
        "rms_cm": float(100 * rms),
        "strike_deg": float(p[2]),
        "length_km": float(p[3] / 1e3),
        "top_km": float(p[4] / 1e3),
        "bottom_km": float(p[5] / 1e3),
        "opening_m": float(p[6]),
        "volume_km3": float(volume),
    }


def main() -> None:
    fit_path = RESULTS / "okada_fit.json"
    if not fit_path.exists():
        raise SystemExit("results/okada_fit.json missing — run `python -m mh26.invert` first")
    default = json.loads(fit_path.read_text())
    base = np.asarray(default["values"], dtype=float)

    rows = [summarize("default", base, default["rms_m"], default["n_obs"])]
    for label, options, equal_weight in VARIANTS[1:]:
        print(f"sensitivity: {label}")
        data = load_observations(**options)
        p, _, rms, n_obs, _ = fit(
            data,
            starts=[base],
            equal_geometry_weight=equal_weight,
        )
        rows.append(summarize(label, p, rms, n_obs))

    keys = ("length_km", "top_km", "bottom_km", "opening_m", "volume_km3")
    envelope = {
        key: {
            "min": min(row[key] for row in rows),
            "max": max(row[key] for row in rows),
        }
        for key in keys
    }
    payload = {
        "description": "processing-choice sensitivity; not confidence intervals",
        "variants": rows,
        "envelope": envelope,
    }
    (RESULTS / "model_sensitivity.json").write_text(json.dumps(payload, indent=2) + "\n")

    lines = [
        "# Uniform-dike model sensitivity",
        "",
        "These runs change one processing choice at a time. They are **not**",
        "confidence intervals and do not account for spatially correlated atmosphere,",
        "faulting above the dike, nonuniform opening, or elastic-structure error.",
        "",
        "| variant | n | RMS (cm) | strike | length (km) | depth (km) | opening (m) | volume (km³) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['variant']} | {row['n_obs']} | {row['rms_cm']:.2f} | "
            f"{row['strike_deg']:.1f}° | {row['length_km']:.2f} | "
            f"{row['top_km']:.2f}–{row['bottom_km']:.2f} | "
            f"{row['opening_m']:.2f} | {row['volume_km3']:.3f} |"
        )
    lines.extend([
        "",
        "Envelope across these runs:",
        "",
        f"- length {envelope['length_km']['min']:.1f}–{envelope['length_km']['max']:.1f} km",
        f"- top depth {envelope['top_km']['min']:.1f}–{envelope['top_km']['max']:.1f} km; "
        f"bottom depth {envelope['bottom_km']['min']:.1f}–{envelope['bottom_km']['max']:.1f} km",
        f"- opening {envelope['opening_m']['min']:.2f}–{envelope['opening_m']['max']:.2f} m",
        f"- geometric volume {envelope['volume_km3']['min']:.3f}–"
        f"{envelope['volume_km3']['max']:.3f} km³",
        "",
        "The volume varies less than depth and opening separately. A statistical uncertainty",
        "analysis requires a spatial covariance model and comparison with nonuniform-opening",
        "and faulting models.",
    ])
    report = "\n".join(lines) + "\n"
    (RESULTS / "model_sensitivity.md").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
