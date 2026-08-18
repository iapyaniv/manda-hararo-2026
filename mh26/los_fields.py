"""Step 1 — LOS displacement panels and the asc/desc E-U decomposition.

Per pair: ionosphere-corrected, far-field-referenced LOS field with peak
location/extent stats. Then the ascending+descending pair is decomposed on
their common grid into quasi-east and quasi-up components using the
product LOS unit vectors.

Run from the repo root:  python -m mh26.los_fields

Writes results/los_panels.png and results/decomposition.png.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .gunw import corrected_los, load_pair, los_unit_neu, pair_window, reference_farfield
from .params import CENTER_E, CENTER_N, PAIRS, RESULTS


def pair_fields() -> dict:
    """Corrected+referenced LOS per pair with peak/extent stats near the centre."""
    results = {}
    for tag, label in PAIRS.items():
        d = load_pair(tag)
        if d is None:
            print(f"{tag}: npz missing, skip")
            continue
        los, _ = corrected_los(d)
        x, y = d["x"], d["y"]
        ref = reference_farfield(los, x, y)
        X, Y = np.meshgrid(x, y)
        R = np.hypot(X - CENTER_E, Y - CENTER_N)
        near = (R < 18_000) & np.isfinite(ref)
        if near.sum() < 50:
            print(f"{tag}: no valid data near target")
            continue
        idx = np.nanargmax(np.abs(np.where(near, ref, np.nan)))
        iy, ix = np.unravel_index(idx, ref.shape)
        cell_area_km2 = abs(float(np.median(np.diff(x))) *
                            float(np.median(np.diff(y)))) / 1e6
        area_1cm = ((np.abs(ref) > 0.01) & (R < 25_000)).sum() * cell_area_km2
        area_3cm = ((np.abs(ref) > 0.03) & (R < 25_000)).sum() * cell_area_km2
        robust_peak_cm = 100 * np.nanpercentile(np.abs(ref[near]), 99.9)
        results[tag] = dict(label=label, los=ref, x=x, y=y,
                            peak_cm=100 * ref[iy, ix], px=x[ix], py=y[iy],
                            robust_peak_cm=robust_peak_cm,
                            area_1cm_km2=area_1cm,
                            area_3cm_km2=area_3cm,
                            window=pair_window(d),
                            valid_frac=float(np.isfinite(ref[R < 18_000]).mean()))
        r = results[tag]
        print(f"{tag} [{label}]: peak {r['peak_cm']:+.1f} cm LOS at "
              f"E{r['px']:.0f} N{r['py']:.0f}; |signal|>1cm area {area_1cm:.0f} km²; "
              f"valid near-field {r['valid_frac']:.0%}")
    return results


def write_los_summary(results: dict) -> None:
    """Persist the numeric values that otherwise only appear on stdout."""
    serializable = {}
    lines = [
        "# NISAR LOS field summary",
        "",
        "All values use ionosphere correction, far-field plane referencing,",
        "coherence/connected-component screening, and the GUNW valid-sample mask.",
        "",
        "| tag | window | peak (cm) | p99.9 | E/N peak (m) | area >1 cm | area >3 cm | valid within 18 km |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for tag, result in results.items():
        values = {
            key: (float(value) if isinstance(value, (float, np.floating)) else value)
            for key, value in result.items()
            if key not in {"los", "x", "y"}
        }
        serializable[tag] = values
        lines.append(
            f"| {tag} | {result['window'] or '?'} | {result['peak_cm']:+.1f} | "
            f"{result['robust_peak_cm']:.1f} | {result['px']:.0f}/{result['py']:.0f} | "
            f"{result['area_1cm_km2']:.0f} km² | {result['area_3cm_km2']:.0f} km² | "
            f"{result['valid_frac']:.0%} |"
        )
    (RESULTS / "los_summary.md").write_text("\n".join(lines) + "\n")
    (RESULTS / "los_summary.json").write_text(json.dumps(serializable, indent=2) + "\n")


def plot_los_panels(results: dict) -> None:
    n = len(results)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5.2 * ncols, 5.5 * nrows), squeeze=False
    )
    flat_axes = axes.ravel()
    for ax, (tag, r) in zip(flat_axes, results.items()):
        v = np.nanpercentile(np.abs(r["los"]), 99.5)
        im = ax.imshow(r["los"] * 100, cmap="RdBu_r", vmin=-100 * v, vmax=100 * v,
                       extent=[r["x"][0] / 1e3, r["x"][-1] / 1e3,
                               r["y"][-1] / 1e3, r["y"][0] / 1e3])
        ax.plot(CENTER_E / 1e3, CENTER_N / 1e3, "k+", ms=12)
        ax.set_title(f"NISAR {r['label']}\npeak {r['peak_cm']:+.0f} cm", fontsize=10)
        ax.set_xlabel("UTM 37N easting (km)", fontsize=8)
        ax.set_ylabel("UTM 37N northing (km)", fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.7,
                     label="LOS displacement (cm, + = toward satellite)")
    for ax in flat_axes[n:]:
        ax.remove()
    fig.suptitle("Manda Hararo segment, Afar (≈12.2°N 40.7°E) — ionosphere-corrected, "
                 "far-field-referenced; + marks 12.18°N 40.66°E", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = RESULTS / "los_panels.png"
    fig.savefig(out, dpi=110)
    print(f"figure: {out}")


def decompose(results: dict) -> None:
    """Solve per-pixel [asc; desc] LOS = M [dE; dU] on the common grid."""
    if not {"172_008_A", "165_082_D"} <= results.keys():
        print("decomposition skipped: need both 172_008_A and 165_082_D")
        return
    a, dsc = results["172_008_A"], results["165_082_D"]
    cxs = a["x"][np.isin(a["x"], dsc["x"])]
    cys = a["y"][np.isin(a["y"], dsc["y"])]
    if len(cxs) <= 10 or len(cys) <= 10:
        print("decomposition skipped: no usable common grid")
        return

    def crop_to(r):
        xm = np.isin(r["x"], cxs)
        ym = np.isin(r["y"], cys)
        return r["los"][np.ix_(np.where(ym)[0], np.where(xm)[0])]

    la, ld = crop_to(a), crop_to(dsc)
    X, Y = np.meshgrid(cxs, cys)
    pts = np.column_stack([Y.ravel(), X.ravel()])
    ga = los_unit_neu(load_pair("172_008_A"), pts).reshape(*X.shape, 3)
    gd = los_unit_neu(load_pair("165_082_D"), pts).reshape(*X.shape, 3)
    gaE, gaU = ga[..., 1], ga[..., 2]
    gdE, gdU = gd[..., 1], gd[..., 2]
    det = gaE * gdU - gaU * gdE
    dE = (gdU * la - gaU * ld) / det
    dU = (-gdE * la + gaE * ld) / det
    print(f"LOS units mean (E,U): asc {np.mean(ga[..., 1:], axis=(0, 1))}, "
          f"desc {np.mean(gd[..., 1:], axis=(0, 1))}")
    print(f"decomposition (local geometry): dU range "
          f"[{100*np.nanmin(dU):+.0f}, {100*np.nanmax(dU):+.0f}] cm, "
          f"dE range [{100*np.nanmin(dE):+.0f}, {100*np.nanmax(dE):+.0f}] cm")

    near = (np.hypot(X - CENTER_E, Y - CENTER_N) < 18_000) & np.isfinite(dE) & np.isfinite(dU)
    if near.sum() < 100:
        raise ValueError("too few common near-field pixels for decomposition summary")
    summary = {
        "assumption": "north displacement is zero; solve east/up only",
        "n_pixels_within_18km": int(near.sum()),
        "east_cm_p01_p99": [float(v) for v in 100 * np.nanpercentile(dE[near], [1, 99])],
        "up_cm_p01_p99": [float(v) for v in 100 * np.nanpercentile(dU[near], [1, 99])],
        "east_cm_min_max": [float(100 * np.nanmin(dE[near])),
                            float(100 * np.nanmax(dE[near]))],
        "up_cm_min_max": [float(100 * np.nanmin(dU[near])),
                          float(100 * np.nanmax(dU[near]))],
    }
    (RESULTS / "decomposition_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (RESULTS / "decomposition_summary.md").write_text(
        "# Ascending/descending East–Up decomposition\n\n"
        "This two-component solution assumes zero north displacement and is reported as\n"
        "quasi-east/quasi-up rather than full 3-D displacement. Statistics\n"
        "below are within 18 km of the analysis centre.\n\n"
        f"- East p01–p99: {summary['east_cm_p01_p99'][0]:+.1f} to "
        f"{summary['east_cm_p01_p99'][1]:+.1f} cm\n"
        f"- Up p01–p99: {summary['up_cm_p01_p99'][0]:+.1f} to "
        f"{summary['up_cm_p01_p99'][1]:+.1f} cm\n"
    )

    # y descends north->south (row 0 = north); imshow extent = [l, r, bottom, top]
    ext = [cxs[0] / 1e3, cxs[-1] / 1e3, min(cys) / 1e3, max(cys) / 1e3]
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.6))
    for axx, F, name in ((ax[0], dU, "quasi-vertical (up +)"),
                         (ax[1], dE, "quasi-east (east +)")):
        v = np.nanpercentile(np.abs(F), 99.5)
        im = axx.imshow(F * 100, cmap="RdBu_r", vmin=-100 * v, vmax=100 * v,
                        extent=ext)
        axx.set_title(name, fontsize=11)
        axx.set_xlabel("UTM 37N easting (km)", fontsize=8)
        axx.set_ylabel("UTM 37N northing (km)", fontsize=8)
        plt.colorbar(im, ax=axx, shrink=0.7, label="displacement (cm)")
    fig.suptitle("Ascending+descending decomposition — Manda Hararo, Afar "
                 "(axial graben subsidence between opposing E–W motions)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = RESULTS / "decomposition.png"
    fig.savefig(out, dpi=110)
    print(f"figure: {out}")


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    results = pair_fields()
    if results:
        write_los_summary(results)
        plot_los_panels(results)
    decompose(results)


if __name__ == "__main__":
    main()
