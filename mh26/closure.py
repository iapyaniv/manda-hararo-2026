"""Loop-closure (phase-triplet) check of the reprocessed NISAR unwrapping.

For each NISAR track the three epochs A < B < C give three interferograms
(A→B, B→C, A→C) on the shared GSLC grid. The unwrapped closure

    C_unw = unw_AB + unw_BC - unw_AC

is compared with the wrapped closure of the same filtered phases

    S = wrap(phi_AB + phi_BC - phi_AC).

snaphu output is congruent to its input modulo 2 pi, so (C_unw - S) / 2 pi
is an exact integer field k. Multilooking and filtering make S itself small
but nonzero; atmosphere, ionosphere, and deformation are epoch differences
and cancel in the loop, so k isolates unwrap inconsistency alone. A global
integer offset of k is unconstrained (each unwrap has an arbitrary
reference); the far-field mode is subtracted before interpretation. This is
an internal-consistency test of the three unwraps against each other, not
an accuracy test of any single one.

Triangles (the B→C leg is the co-event or straddle pair the analysis uses):

    172_008_A          06-14 -> 07-08 -> 07-20   (co-event, inversion pair)
    165_082_D          06-13 -> 07-07 -> 07-19   (co-event, inversion pair)
    100_A_008_STRADDLE 07-03 -> 07-15 -> 07-27   (straddle, timing pair)

The cycle-023 GSLC epochs (06-25/06-26) are absent from the provisional
archive, so the ascending/descending pre legs span 24 days. Stages::

    python -m mh26.closure --crop      # fetch the one extra epoch per track
    python -m mh26.closure --ifg       # form the pre and skip interferograms
    python -m mh26.closure --check     # snaphu x3, closure, stats, figure

The pre and skip interferograms use the co-event chain's exact multilook,
filter, coherence, and snaphu configuration (mh26.gslc), and the co-event
leg is unwrapped from the saved interferogram file, so the checked unwraps
are the analysis unwraps.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from .gslc import (
    COHERENCE_WINDOW,
    DATA_GSLC,
    LOOKS,
    PAIRS,
    SNAPHU_PRIMARY,
    UNWRAP_COHERENCE_MIN,
    _axis_distance,
    _epoch_invalid,
    _extract_crop,
    _goldstein,
    _multilook,
    _remote_file,
    _run_snaphu,
    NEAR_AXIS_ALONG_M,
    NEAR_AXIS_PERP_M,
)
from .params import CENTER_E, CENTER_N, FARFIELD_RADIUS, RESULTS, WAVELENGTH_L

# Extra (earliest) GSLC epoch per track. The other two epochs are the
# co-event crops already on disk ({tag}_ref.npz, {tag}_sec.npz).
PRE_GRANULES = {
    "172_008_A": "NISAR_L2_PR_GSLC_022_172_A_008_2005_DHDH_A_"
                 "20260614T024620_20260614T024655_P05023_N_F_J_001",
    "165_082_D": "NISAR_L2_PR_GSLC_022_165_D_082_2005_SHSH_A_"
                 "20260613T154913_20260613T154947_P05023_N_F_J_001",
    "100_A_008_STRADDLE": "NISAR_L2_PR_GSLC_024_100_A_008_2005_DHDH_A_"
                          "20260703T025438_20260703T025513_P05023_N_F_J_001",
}
# Loop legs per track: (epoch pair, interferogram file suffix).
LEGS = ("pre", "co", "skip")
ONE_CYCLE_LOS_M = WAVELENGTH_L / 2


def _date(granule: str) -> str:
    for token in granule.split("_"):
        if len(token) > 8 and token[:8].isdigit() and "T" in token:
            return f"{token[0:4]}-{token[4:6]}-{token[6:8]}"
    raise ValueError(f"no acquisition date in granule name: {granule}")


def triangle_dates(tag: str) -> tuple[str, str, str]:
    ref, sec = PAIRS[tag]
    return _date(PRE_GRANULES[tag]), _date(ref), _date(sec)


def crop(tags: list[str] | None = None) -> None:
    import earthaccess

    earthaccess.login(strategy="netrc")
    for tag, granule in PRE_GRANULES.items():
        if tags and tag not in tags:
            continue
        out = DATA_GSLC / f"{tag}_pre.npz"
        if out.exists():
            print(f"{tag} pre: already extracted")
            continue
        print(f"{tag} pre: {granule}")
        with _remote_file(granule) as f:
            _extract_crop(f, granule, out)


def _form_pair(z1_npz, z2_npz, out, granules: tuple[str, str],
               looks: int = LOOKS) -> None:
    """One interferogram with the co-event chain's multilook/filter/coherence."""
    from scipy.ndimage import uniform_filter

    looks_y, looks_x = 2 * looks, looks
    if not (np.array_equal(z1_npz["x"], z2_npz["x"])
            and np.array_equal(z1_npz["y"], z2_npz["y"])):
        raise ValueError(f"{out.name}: epoch grids differ")
    z1, z2 = z1_npz["z"], z2_npz["z"]
    for z, npz in ((z1, z1_npz), (z2, z2_npz)):
        z[_epoch_invalid(npz, z.shape)] = 0
    ny, nx = z1.shape
    by, bx = ny // looks_y, nx // looks_x
    ifg = np.empty((by, bx), dtype=np.complex64)
    p1 = np.empty((by, bx), dtype=np.float32)
    p2 = np.empty((by, bx), dtype=np.float32)
    band_rows = 300
    for b0 in range(0, by, band_rows):
        b1 = min(b0 + band_rows, by)
        sl = np.s_[b0 * looks_y:b1 * looks_y, :bx * looks_x]
        a, b = z1[sl], z2[sl]
        ifg[b0:b1] = _multilook(a * np.conj(b), looks_y, looks_x)
        p1[b0:b1] = _multilook((a.real ** 2 + a.imag ** 2), looks_y, looks_x)
        p2[b0:b1] = _multilook((b.real ** 2 + b.imag ** 2), looks_y, looks_x)
    with np.errstate(invalid="ignore", divide="ignore"):
        denominator = np.sqrt(uniform_filter(p1, COHERENCE_WINDOW)
                              * uniform_filter(p2, COHERENCE_WINDOW))
        smoothed = (uniform_filter(ifg.real, COHERENCE_WINDOW)
                    + 1j * uniform_filter(ifg.imag, COHERENCE_WINDOW))
        coherence = np.abs(smoothed) / denominator
    coherence = np.clip(np.nan_to_num(coherence), 0.0, 1.0)
    filtered = _goldstein(ifg / np.maximum(np.abs(ifg), 1e-30))
    x = z1_npz["x"][:bx * looks_x].reshape(-1, looks_x).mean(axis=1)
    y = z1_npz["y"][:by * looks_y].reshape(-1, looks_y).mean(axis=1)
    np.savez_compressed(out, phase=np.angle(filtered).astype(np.float32),
                        coherence=coherence.astype(np.float32),
                        x=x, y=y, epsg=z1_npz["epsg"],
                        looks=np.array([looks_y, looks_x]),
                        granules=np.array(list(granules)))
    print(f"  saved {out.name} {filtered.shape}")


def make_ifgs(tags: list[str] | None = None) -> None:
    for tag in PRE_GRANULES:
        if tags and tag not in tags:
            continue
        pre = np.load(DATA_GSLC / f"{tag}_pre.npz")
        ref = np.load(DATA_GSLC / f"{tag}_ref.npz")
        sec = np.load(DATA_GSLC / f"{tag}_sec.npz")
        pre_granule = PRE_GRANULES[tag]
        ref_granule, sec_granule = PAIRS[tag]
        for suffix, first, second, granules in (
                ("pre", pre, ref, (pre_granule, ref_granule)),
                ("skip", pre, sec, (pre_granule, sec_granule))):
            out = DATA_GSLC / f"{tag}_loop{suffix}_ifg.npz"
            if out.exists():
                print(f"{tag} {suffix}: already formed")
                continue
            print(f"{tag} {suffix} leg:")
            _form_pair(first, second, out, granules)


def _leg_file(tag: str, leg: str):
    if leg == "co":
        return DATA_GSLC / f"{tag}_ifg.npz"
    return DATA_GSLC / f"{tag}_loop{leg}_ifg.npz"


def _mode(values: np.ndarray) -> int:
    ints, counts = np.unique(values, return_counts=True)
    return int(ints[np.argmax(counts)])


def check(tags: list[str] | None = None) -> None:
    model = np.asarray(
        json.loads((RESULTS / "okada_fit.json").read_text())["values"][:7])

    payload = {}
    panels = {}
    for tag in PRE_GRANULES:
        if tags and tag not in tags:
            continue
        legs = {}
        for leg in LEGS:
            f = np.load(_leg_file(tag, leg))
            legs[leg] = f
        x, y = legs["co"]["x"], legs["co"]["y"]
        for leg in ("pre", "skip"):
            if not (np.array_equal(legs[leg]["x"], x)
                    and np.array_equal(legs[leg]["y"], y)):
                raise ValueError(f"{tag}: {leg} grid differs from the co grid")

        unw, phases, coherences, unwrap_configs = {}, {}, {}, {}
        for leg in LEGS:
            phases[leg] = np.asarray(legs[leg]["phase"], dtype=np.float64)
            coherences[leg] = np.asarray(legs[leg]["coherence"])
            # Raw (unreferenced, no iono) unwraps are cached so the combined
            # report does not repeat the snaphu runs.
            cache = DATA_GSLC / f"{tag}_loop{leg}_unwraw.npz"
            if cache.exists():
                stored = np.load(cache)
                cached = np.asarray(stored["unw"], dtype=np.float64)
                if cached.shape != phases[leg].shape:
                    raise ValueError(f"{tag} {leg}: cached unwrap shape "
                                     "differs from the interferogram")
                unwrap_configs[leg] = (json.loads(str(stored["config"]))
                                       if "config" in stored.files
                                       else dict(SNAPHU_PRIMARY))
                print(f"{tag} {leg}: cached raw unwrap", flush=True)
                unw[leg] = cached
                continue
            print(f"{tag} {leg}: snaphu on {phases[leg].shape}", flush=True)
            try:
                config = dict(SNAPHU_PRIMARY)
                raw = _run_snaphu(phases[leg].astype(np.float32),
                                  coherences[leg], **config)
            except (RuntimeError, MemoryError) as error:
                # The 36-day skip legs can carry enough residues that the
                # MCF initialization exhausts memory; the tree
                # initialization is far lighter and the full DEFO
                # optimization still follows, so the output stays congruent
                # to the input modulo 2 pi (all the closure needs).
                print(f"{tag} {leg}: MCF initialization failed "
                      f"({str(error).splitlines()[0][:120]}); retrying with "
                      "MST", flush=True)
                config = {"stat_cost_mode": "DEFO", "init_method": "MST"}
                raw = _run_snaphu(phases[leg].astype(np.float32),
                                  coherences[leg], **config)
            unwrap_configs[leg] = config
            np.savez_compressed(
                cache, unw=raw, config=np.array(json.dumps(config)),
                note=np.array("raw snaphu radians; no ionosphere or "
                              "far-field referencing (loop-closure use)"))
            unw[leg] = raw.astype(np.float64)

        wrapped_sum = np.angle(
            np.exp(1j * (phases["pre"] + phases["co"] - phases["skip"])))
        closure = unw["pre"] + unw["co"] - unw["skip"]
        k = np.rint((closure - wrapped_sum) / (2 * np.pi))
        congruence = np.max(np.abs(closure - wrapped_sum - 2 * np.pi * k))

        mask = np.ones(closure.shape, dtype=bool)
        for leg in LEGS:
            mask &= coherences[leg] > UNWRAP_COHERENCE_MIN
        X, Y = np.meshgrid(x, y)
        radius = np.hypot(X - CENTER_E, Y - CENTER_N)
        far = mask & (radius > FARFIELD_RADIUS)
        k0 = _mode(k[far])
        k_rel = np.where(mask, k - k0, np.nan)
        cycles = np.where(mask, closure / (2 * np.pi) - k0, np.nan)

        perp, along = _axis_distance(X, Y, model)
        near_axis = (perp < NEAR_AXIS_PERP_M) & (along < NEAR_AXIS_ALONG_M)
        pixel_km2 = abs(float(np.median(np.diff(x)))
                        * float(np.median(np.diff(y)))) / 1e6

        def region_stats(region: np.ndarray) -> dict:
            sel = mask & region
            if not sel.any():
                return {"n_pixels": 0}
            kk = k_rel[sel]
            return {
                "n_pixels": int(sel.sum()),
                "fraction_error_ge_1_cycle": float(np.mean(np.abs(kk) >= 1)),
                "area_error_ge_1_cycle_km2": float(
                    np.sum(np.abs(kk) >= 1) * pixel_km2),
                "max_abs_error_cycles": int(np.max(np.abs(kk))),
                "rms_wrapped_nonclosure_cycles": float(np.sqrt(np.mean(
                    (wrapped_sum[sel] / (2 * np.pi)) ** 2))),
            }

        error_region = mask & (np.abs(np.where(mask, k - k0, 0)) >= 1)
        attribution = None
        if error_region.any():
            attribution = {
                "median_coherence_in_error_region": {
                    leg: float(np.median(coherences[leg][error_region]))
                    for leg in LEGS},
                "median_coherence_elsewhere": {
                    leg: float(np.median(coherences[leg][mask & ~error_region]))
                    for leg in LEGS},
                "median_abs_perp_km": float(
                    np.median(perp[error_region]) / 1e3),
            }

        dates = triangle_dates(tag)
        stats = {
            "triangle": {"epochs": dates,
                         "legs": {"pre": f"{dates[0]} -> {dates[1]}",
                                  "co": f"{dates[1]} -> {dates[2]}",
                                  "skip": f"{dates[0]} -> {dates[2]}"}},
            "unwrap_configuration": unwrap_configs,
            "error_region": attribution,
            "farfield_mode_cycles": k0,
            "congruence_max_abs_rad": float(congruence),
            "whole_scene": region_stats(np.ones_like(mask)),
            "near_axis": region_stats(near_axis),
            "off_axis_within_25km": region_stats(~near_axis
                                                 & (radius < 25_000)),
            "one_cycle_los_m": ONE_CYCLE_LOS_M,
        }
        payload[tag] = stats
        w = stats["whole_scene"]
        n = stats["near_axis"]
        print(f"{tag}: k0 {k0}; >=1-cycle fraction whole "
              f"{w['fraction_error_ge_1_cycle']:.2%} "
              f"({w['area_error_ge_1_cycle_km2']:.2f} km2), near-axis "
              f"{n['fraction_error_ge_1_cycle']:.2%} "
              f"({n['area_error_ge_1_cycle_km2']:.2f} km2), max "
              f"{w['max_abs_error_cycles']} cycles", flush=True)

        display = {}
        for leg in LEGS:
            u = np.where(coherences[leg] > UNWRAP_COHERENCE_MIN,
                         unw[leg], np.nan)
            display[leg] = u - np.nanmedian(u[far])
        panels[tag] = (display, np.abs(cycles), x, y)

    for tag, stats in payload.items():
        (RESULTS / f"loop_closure_{tag}.json").write_text(
            json.dumps(stats, indent=2) + "\n")

    if len(payload) == len(PRE_GRANULES):
        write_md(payload)
    if panels:
        figure(panels, payload)


def write_md(payload: dict) -> None:
    lines = [
        "# Loop-closure check of the reprocessed NISAR unwrapping", "",
        "Per track, the three epochs give three interferograms on the "
        "shared GSLC grid",
        "(pre, co-event, and the skip pair that spans both). The unwrapped "
        "closure",
        "`unw_pre + unw_co − unw_skip` minus the wrapped closure of the "
        "same filtered",
        "phases is an exact integer number of cycles per pixel (snaphu is "
        "congruent",
        "modulo 2π). Deformation, atmosphere, and ionosphere are epoch "
        "differences and",
        "cancel in the loop, so a nonzero integer isolates unwrap "
        "inconsistency among",
        "the three legs. The arbitrary per-unwrap reference makes a "
        "constant integer",
        "offset unobservable; the far-field mode is removed. One cycle is "
        f"λ/2 = {100 * ONE_CYCLE_LOS_M:.1f} cm of LOS. This is an "
        "internal-consistency test",
        "of the three unwraps against each other, not an accuracy test of "
        "any one of",
        "them, and the co-event legs are unwrapped from the exact saved "
        "interferograms",
        "the analysis uses (same snaphu configuration).", "",
    ]
    for tag, stats in payload.items():
        legs = stats["triangle"]["legs"]
        lines += [f"## {tag}", "",
                  f"- triangle: pre {legs['pre']}, co {legs['co']}, "
                  f"skip {legs['skip']}"]
        for region, label in (("whole_scene", "whole coherent scene"),
                              ("near_axis", "near-axis band (|⊥| < 3 km, "
                                            "|∥| < 20 km)"),
                              ("off_axis_within_25km", "off-axis within "
                                                       "25 km")):
            r = stats[region]
            if not r.get("n_pixels"):
                lines.append(f"- {label}: no pixels pass the mask")
                continue
            lines.append(
                f"- {label}: {r['fraction_error_ge_1_cycle']:.2%} of "
                f"{r['n_pixels']:,} pixels off by ≥ 1 cycle "
                f"({r['area_error_ge_1_cycle_km2']:.2f} km²; max "
                f"{r['max_abs_error_cycles']} cycles); wrapped nonclosure "
                f"RMS {r['rms_wrapped_nonclosure_cycles']:.3f} cycles")
        attribution = stats.get("error_region")
        if attribution:
            in_error = attribution["median_coherence_in_error_region"]
            elsewhere = attribution["median_coherence_elsewhere"]
            lines.append(
                "- error-region attribution: median |⊥ distance| "
                f"{attribution['median_abs_perp_km']:.1f} km; median "
                "coherence there "
                + ", ".join(f"{leg} {in_error[leg]:.2f}" for leg in LEGS)
                + " (elsewhere "
                + ", ".join(f"{elsewhere[leg]:.2f}" for leg in LEGS)
                + ") — the loop flags the triple, and cannot assign the "
                "error to one leg by itself")
        lines.append("")
    lines += [
        "## Scope", "",
        "- The mask requires coherence above the analysis threshold "
        f"({UNWRAP_COHERENCE_MIN}) in all three legs, so the checked area "
        "is the intersection.",
        "- The wrapped-closure RMS is the multilook/filter nonclosure "
        "floor; it is not an unwrap error.",
        "- A consistent triple can still share one common unwrapping "
        "error; the GUNW comparison (results/gslc_gunw_compare_*.json) "
        "and the Sentinel-1 cross-checks bound that case independently.",
        "- The pre legs span 24 days on tracks 172/165 (the 12-day "
        "cycle-023 GSLC epochs are absent from the provisional archive) "
        "and 12 days on track 100.",
    ]
    out = RESULTS / "loop_closure.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


def figure(panels: dict, payload: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    n_rows = len(panels)
    fig, axes = plt.subplots(n_rows, 4, figsize=(16, 4.6 * n_rows),
                             squeeze=False)
    phase_limit = 100.0
    cycles_limit = 5.0
    phase_cmap = matplotlib.colormaps["RdBu_r"].copy()
    cycles_cmap = matplotlib.colormaps["Reds"].copy()
    for cmap in (phase_cmap, cycles_cmap):
        cmap.set_bad("#b3b3b3")
    label = {"172_008_A": "NISAR asc 172", "165_082_D": "NISAR desc 165",
             "100_A_008_STRADDLE": "NISAR asc 100"}
    for row, (tag, (display, cycles, x, y)) in enumerate(panels.items()):
        extent = [x[0] / 1e3, x[-1] / 1e3, y[-1] / 1e3, y[0] / 1e3]
        legs = payload[tag]["triangle"]["legs"]
        for column, leg in enumerate(LEGS):
            ax = axes[row, column]
            im = ax.imshow(display[leg], cmap=phase_cmap, extent=extent,
                           vmin=-phase_limit, vmax=phase_limit)
            ax.set_title(f"{label[tag]}: {legs[leg]}", fontsize=10)
            if column == 0:
                ax.set_ylabel("UTM 37N northing (km)", fontsize=9)
        ax = axes[row, 3]
        imc = ax.imshow(cycles, cmap=cycles_cmap, extent=extent,
                        vmin=0, vmax=cycles_limit)
        stats = payload[tag]["whole_scene"]
        ax.set_title(
            f"loop closure: {stats['fraction_error_ge_1_cycle']:.2%} of "
            f"pixels ≥ 1 cycle", fontsize=10)
        if row == n_rows - 1:
            for column in range(4):
                axes[row, column].set_xlabel("UTM 37N easting (km)",
                                             fontsize=9)
    fig.subplots_adjust(left=0.05, right=0.90, top=0.94, bottom=0.07,
                        hspace=0.28, wspace=0.16)
    cax1 = fig.add_axes([0.915, 0.55, 0.013, 0.35])
    fig.colorbar(matplotlib.cm.ScalarMappable(
        norm=Normalize(-phase_limit, phase_limit), cmap="RdBu_r"),
        cax=cax1, label="unwrapped phase (rad; far-field median removed)")
    cax2 = fig.add_axes([0.915, 0.10, 0.013, 0.35])
    fig.colorbar(matplotlib.cm.ScalarMappable(
        norm=Normalize(0, cycles_limit), cmap="Reds"),
        cax=cax2, label="absolute loop-closure error (cycles)")
    fig.suptitle("NISAR 30-m reprocessed chain: unwrapped-phase loop closure "
                 "(gray = below the coherence threshold in any leg)",
                 fontsize=12)
    out = RESULTS / "loop_closure_nisar.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"figure: {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crop", action="store_true")
    ap.add_argument("--ifg", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--tags", nargs="*")
    args = ap.parse_args()
    if args.crop:
        crop(tags=args.tags)
    if args.ifg:
        make_ifgs(tags=args.tags)
    if args.check:
        check(tags=args.tags)
    if not any((args.crop, args.ifg, args.check)):
        ap.print_help()


if __name__ == "__main__":
    main()
