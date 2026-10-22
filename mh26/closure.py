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
    python -m mh26.closure --check --cached-only  # read-only saved replay
    python -m mh26.closure --check --cached-only --regenerate  # replace reports
    python -m mh26.closure --check --regenerate  # can run missing snaphu unwraps

The pre and skip interferograms use the co-event chain's multilook, filter,
and coherence settings (mh26.gslc). Legacy raw caches have no snaphu
configuration metadata; their settings are unknown. Report generation records
exact file hashes. Replay checks each raw unwrap against its phase input. Once an input
inventory exists, changed cached inputs are rejected. Default replay also
requires the saved source, settings, and reports to match before calculation.
Only --regenerate can replace those reports and their inventory. This
reproduces the closure calculation from saved arrays, not the original
snaphu execution.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .gslc import (
    earthdata_login,
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
from .params import (CENTER_E, CENTER_N, EPSG_UTM37N, FARFIELD_RADIUS,
                     REPO, RESULTS, WAVELENGTH_L)

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
OFF_AXIS_RADIUS_M = 25_000.0
CONGRUENCE_TOLERANCE_RAD = 1e-3


def file_record(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path.relative_to(REPO)),
            "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def analysis_parameters() -> dict:
    return {"coherence_min": UNWRAP_COHERENCE_MIN,
            "centre_e_m": CENTER_E, "centre_n_m": CENTER_N,
            "farfield_radius_m": FARFIELD_RADIUS,
            "near_axis_half_width_m": NEAR_AXIS_PERP_M,
            "near_axis_half_length_m": NEAR_AXIS_ALONG_M,
            "off_axis_radius_m": OFF_AXIS_RADIUS_M,
            "one_cycle_los_m": ONE_CYCLE_LOS_M,
            "congruence_tolerance_rad": CONGRUENCE_TOLERANCE_RAD,
            "ifg_epsg": EPSG_UTM37N, "ifg_looks": [2 * LOOKS, LOOKS],
            "ifg_posting_m": [10 * LOOKS, -10 * LOOKS],
            "triangle_granules": {tag: [pre, *PAIRS[tag]]
                                  for tag, pre in PRE_GRANULES.items()}}


def write_inventory(payload: dict) -> None:
    inputs = {record["path"]: record
              for stats in payload.values()
              for leg in stats["input_files"].values()
              for record in leg.values()}
    reports = [RESULTS / f"loop_closure_{tag}.json" for tag in payload]
    reports += [RESULTS / "loop_closure.md", RESULTS / "loop_closure_nisar.png"]
    inventory = {
        "schema_version": 1,
        "generated_by": "mh26.closure",
        "replay_command": "python -m mh26.closure --check --cached-only",
        "regenerate_command": "python -m mh26.closure --check --cached-only --regenerate",
        "scope": "Exact file identities for closure replay from saved arrays. "
                 "The closure data add-on contains the nine raw caches and "
                 "nine interferograms; other archives are not sufficient. Unknown "
                 "legacy snaphu settings prevent exact reconstruction of "
                 "the original unwrapping run from its configuration. Replay "
                 "is read-only and checks the saved baseline; regeneration "
                 "requires an explicit option. Legacy executable identities "
                 "and versions are also unknown.",
        "analysis_parameters": analysis_parameters(),
        "numpy_version": np.__version__,
        "source_files": [file_record(REPO / path) for path in (
            "mh26/closure.py", "mh26/gslc.py", "mh26/params.py",
            "results/okada_fit.json")],
        "input_files": [inputs[path] for path in sorted(inputs)],
        "reports": [file_record(path) for path in reports],
    }
    out = RESULTS / "loop_closure_inventory.json"
    out.write_text(json.dumps(inventory, indent=2, allow_nan=False) + "\n")
    checksum_path = REPO / "data" / "SHA256SUMS-closure.txt"
    checksum_path.write_text("\n".join(sorted(
        f"{record['sha256']}  {Path(record['path']).relative_to('data')}"
        for record in inputs.values())) + "\n")
    print(f"wrote {out}")


def _date(granule: str) -> str:
    for token in granule.split("_"):
        if len(token) > 8 and token[:8].isdigit() and "T" in token:
            return f"{token[0:4]}-{token[4:6]}-{token[6:8]}"
    raise ValueError(f"no acquisition date in granule name: {granule}")


def triangle_dates(tag: str) -> tuple[str, str, str]:
    ref, sec = PAIRS[tag]
    return _date(PRE_GRANULES[tag]), _date(ref), _date(sec)


def crop(tags: list[str] | None = None) -> None:
    earthdata_login()
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


def check(tags: list[str] | None = None, *, cached_only: bool = True,
          regenerate: bool = False) -> dict:
    """Verify the saved baseline and replay cached inputs without writing.

    ``regenerate=True`` replaces all reports and their inventory after an
    authorized analysis change. Existing input hashes remain pinned. Only
    regeneration with ``cached_only=False`` can run missing snaphu unwraps.
    """
    if tags and not set(tags).issubset(PRE_GRANULES):
        raise ValueError("unknown closure track")
    if regenerate:
        if tags and set(tags) != set(PRE_GRANULES):
            raise ValueError("regeneration requires all closure tracks")
        return _compute(tags, cached_only=cached_only, write_outputs=True)
    from .validate import validate_closure

    return validate_closure(require_inputs=True, tags=tags)["recomputed"]


def _compute(tags: list[str] | None = None, *, cached_only: bool = True,
             write_outputs: bool = False) -> dict:
    model = np.asarray(
        json.loads((RESULTS / "okada_fit.json").read_text())["values"][:7])

    payload = {}
    panels = {}
    complete_run = not tags or set(tags) == set(PRE_GRANULES)
    inventory_path = RESULTS / "loop_closure_inventory.json"
    recorded_inputs = {}
    if inventory_path.exists():
        inventory = json.loads(inventory_path.read_text())
        recorded_inputs = {record["path"]: record
                           for record in inventory["input_files"]}
    for tag in PRE_GRANULES:
        if tags and tag not in tags:
            continue
        legs = {}
        pre, (ref, sec) = PRE_GRANULES[tag], PAIRS[tag]
        expected_granules = {"pre": (pre, ref), "co": (ref, sec), "skip": (pre, sec)}
        for leg in LEGS:
            f = np.load(_leg_file(tag, leg))
            if not {"phase", "coherence", "x", "y"}.issubset(f.files):
                raise ValueError(f"{tag} {leg}: incomplete interferogram arrays")
            for axis, step in (("x", 10 * LOOKS), ("y", -10 * LOOKS)):
                coordinate = f[axis]
                if coordinate.ndim != 1 or coordinate.size < 2 \
                        or not np.all(np.isfinite(coordinate)) \
                        or not np.allclose(np.diff(coordinate), step, rtol=0, atol=1e-6):
                    raise ValueError(f"{tag} {leg}: invalid interferogram {axis} grid")
            shape = (len(f["y"]), len(f["x"]))
            for field in ("phase", "coherence"):
                if f[field].ndim != 2 or f[field].shape != shape:
                    raise ValueError(f"{tag} {leg}: interferogram {field} shape "
                                     "does not match its grid")
            for field, expected in (("granules", expected_granules[leg]),
                                    ("epsg", EPSG_UTM37N),
                                    ("looks", [2 * LOOKS, LOOKS])):
                if field not in f.files or not np.array_equal(f[field], expected):
                    raise ValueError(f"{tag} {leg}: interferogram {field} "
                                     "does not match the closure inputs")
            legs[leg] = f
        x, y = legs["co"]["x"], legs["co"]["y"]
        for leg in ("pre", "skip"):
            if not (np.array_equal(legs[leg]["x"], x)
                    and np.array_equal(legs[leg]["y"], y)):
                raise ValueError(f"{tag}: {leg} grid differs from the co grid")

        unw, phases, coherences, unwrap_configs = {}, {}, {}, {}
        config_sources, input_files, leg_congruence, snaphu_executions = {}, {}, {}, {}
        for leg in LEGS:
            phases[leg] = np.asarray(legs[leg]["phase"], dtype=np.float64)
            coherences[leg] = np.asarray(legs[leg]["coherence"])
            ifg_record = file_record(_leg_file(tag, leg))
            # Raw (unreferenced, no iono) unwraps are cached so the combined
            # report does not repeat the snaphu runs.
            cache = DATA_GSLC / f"{tag}_loop{leg}_unwraw.npz"
            if cache.exists():
                cache_record = file_record(cache)
                for record in (ifg_record, cache_record):
                    previous = recorded_inputs.get(record["path"])
                    if previous is not None and record != previous:
                        raise ValueError(f"{record['path']}: differs from "
                                         "the recorded closure input")
                    if inventory_path.exists() and previous is None:
                        raise ValueError(f"{record['path']}: missing from "
                                         "the closure input inventory")
                with np.load(cache, allow_pickle=False) as stored:
                    cached = np.asarray(stored["unw"], dtype=np.float64)
                    if cached.shape != phases[leg].shape:
                        raise ValueError(f"{tag} {leg}: cached unwrap shape "
                                         "differs from the interferogram")
                    if "ifg_sha256" in stored.files and str(
                            stored["ifg_sha256"]) != ifg_record["sha256"]:
                        raise ValueError(f"{tag} {leg}: cached unwrap "
                                         "records a different interferogram")
                    unwrap_configs[leg] = (json.loads(str(stored["config"]))
                                           if "config" in stored.files
                                           else None)
                    snaphu_executions[leg] = (
                        json.loads(str(stored["snaphu_execution"]))
                        if "snaphu_execution" in stored.files else None)
                config_sources[leg] = ("cache config metadata"
                                       if unwrap_configs[leg] is not None
                                       else "unknown: legacy cache has no "
                                            "config metadata")
                print(f"{tag} {leg}: cached raw unwrap", flush=True)
                unw[leg] = cached
                input_files[leg] = {"interferogram": ifg_record,
                                    "raw_unwrap_cache": cache_record}
                continue
            if cached_only or not write_outputs:
                raise FileNotFoundError(f"{cache}: saved raw unwrap required; "
                                        "snaphu is disabled for this replay")
            print(f"{tag} {leg}: snaphu on {phases[leg].shape}", flush=True)
            execution = {}
            try:
                config = dict(SNAPHU_PRIMARY)
                raw = _run_snaphu(phases[leg].astype(np.float32),
                                  coherences[leg], execution_identity=execution, **config)
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
                                  coherences[leg], execution_identity=execution, **config)
            unwrap_configs[leg] = config
            config_sources[leg] = "cache config metadata"
            snaphu_executions[leg] = execution
            np.savez_compressed(
                cache, unw=raw, config=np.array(json.dumps(config)),
                ifg_sha256=np.array(ifg_record["sha256"]),
                snaphu_execution=np.array(json.dumps(execution)),
                note=np.array("raw snaphu radians; no ionosphere or "
                              "far-field referencing (loop-closure use)"))
            unw[leg] = raw.astype(np.float64)
            input_files[leg] = {"interferogram": ifg_record,
                                "raw_unwrap_cache": file_record(cache)}

        for leg in LEGS:
            difference = unw[leg] - phases[leg]
            wrapped = difference - 2 * np.pi * np.rint(difference / (2 * np.pi))
            error = float(np.max(np.abs(wrapped)))
            if not np.isfinite(error) or error > CONGRUENCE_TOLERANCE_RAD:
                raise ValueError(f"{tag} {leg}: raw unwrap is not congruent "
                                 f"with its phase input ({error:.6g} rad)")
            leg_congruence[leg] = error

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
                "n_error_ge_1_cycle": int(np.count_nonzero(np.abs(kk) >= 1)),
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
            "triangle": {"epochs": list(dates),
                         "legs": {"pre": f"{dates[0]} -> {dates[1]}",
                                  "co": f"{dates[1]} -> {dates[2]}",
                                  "skip": f"{dates[0]} -> {dates[2]}"}},
            "unwrap_configuration": unwrap_configs,
            "unwrap_configuration_source": config_sources,
            "snaphu_execution": snaphu_executions,
            "input_files": input_files,
            "leg_congruence_max_abs_rad": leg_congruence,
            "error_region": attribution,
            "farfield_mode_cycles": k0,
            "congruence_max_abs_rad": float(congruence),
            "whole_scene": region_stats(np.ones_like(mask)),
            "near_axis": region_stats(near_axis),
            "off_axis_within_25km": region_stats(~near_axis
                                                 & (radius < OFF_AXIS_RADIUS_M)),
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

        if write_outputs and complete_run:
            display = {}
            for leg in LEGS:
                u = np.where(coherences[leg] > UNWRAP_COHERENCE_MIN,
                             unw[leg], np.nan)
                display[leg] = u - np.nanmedian(u[far])
            panels[tag] = (display, np.abs(cycles), x, y)

    if not write_outputs:
        return payload
    if not complete_run:
        print("Partial check: combined reports and input inventory were "
              "not replaced; check all tracks to save a complete result.")
        return payload
    for tag, stats in payload.items():
        (RESULTS / f"loop_closure_{tag}.json").write_text(
            json.dumps(stats, indent=2, allow_nan=False) + "\n")

    if len(payload) == len(PRE_GRANULES):
        write_md(payload)
    if panels:
        figure(panels, payload)
    if len(payload) == len(PRE_GRANULES):
        write_inventory(payload)
    return payload


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
        "them. The co-event phase inputs are the saved analysis "
        "interferograms.",
        "The raw caches are checked for phase congruence; their snaphu "
        "settings are reported only when cache metadata records them.",
        "The near-axis band "
        "has a 3-km half-width (6 km wide) and a 20-km along-axis "
        "half-length.", "",
    ]
    for tag, stats in payload.items():
        legs = stats["triangle"]["legs"]
        lines += [f"## {tag}", "",
                  f"- triangle: pre {legs['pre']}, co {legs['co']}, "
                  f"skip {legs['skip']}"]
        configurations = stats["unwrap_configuration"]
        lines.append("- cached snaphu settings: " + "; ".join(
            f"{leg} " + (json.dumps(configurations[leg], sort_keys=True)
                         if configurations[leg] is not None
                         else "unknown (no config metadata)")
            for leg in LEGS))
        lines.append("- cached snaphu executable identity and version: " + "; ".join(
            f"{leg} " + (json.dumps(stats['snaphu_execution'][leg], sort_keys=True)
                         if stats['snaphu_execution'][leg] is not None else "unknown")
            for leg in LEGS))
        lines.append("- maximum raw/phase congruence error: " + ", ".join(
            f"{leg} {stats['leg_congruence_max_abs_rad'][leg]:.3g} rad"
            for leg in LEGS))
        for region, label in (("whole_scene", "whole coherent scene"),
                              ("near_axis", "near-axis band (|⊥| < 3 km, "
                                            "|∥| < 20 km)"),
                              ("off_axis_within_25km", "off-axis within "
                                                       "25 km")):
            r = stats[region]
            if not r.get("n_pixels"):
                lines.append(f"- {label}: no pixels pass the mask")
                continue
            error_count = r["n_error_ge_1_cycle"]
            precision = 4 if region == "off_axis_within_25km" else 2
            lines.append(
                f"- {label}: {r['fraction_error_ge_1_cycle']:.{precision}%} "
                f"({error_count:,} of {r['n_pixels']:,} pixels) off by ≥ 1 "
                f"cycle ({r['area_error_ge_1_cycle_km2']:.{precision}f} km²; max "
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
        "- `loop_closure_inventory.json` records SHA-256 hashes, file "
        "sizes, analysis settings, source files, and report identities. "
        "Its 18 interferogram/cache inputs reproduce closure with "
        "`python -m mh26.closure --check --cached-only`, without snaphu or "
        "file writes. Replay first verifies the saved source, settings, "
        "reports, and input hashes, then checks each interferogram's "
        "ordered granules, EPSG, looks, array shape, 30-m grid, and raw-array "
        "phase congruence. Explicit "
        "`--regenerate` replaces reports and the inventory after an "
        "authorized analysis change; input hashes remain pinned. "
        "The original settings of legacy caches are unknown; the current "
        "default is not evidence of those settings. Missing optional raw "
        "inputs prevent this replay even when the archive checks pass.",
        "- The mask requires coherence above the analysis threshold "
        f"({UNWRAP_COHERENCE_MIN}) in all three legs, so the checked area "
        "is the intersection. Flagged pixels therefore pass that threshold; "
        "coherence screening does not remove every unwrap inconsistency. "
        "This check does not show that all flagged axial support was "
        "excluded from the model fits.",
        "- The wrapped-closure RMS is the multilook/filter nonclosure "
        "floor; it is not an unwrap error.",
        "- Closure inconsistency is nearly zero outside the near-axis "
        "band: less than 0.005% on each track within 25 km. Inside the "
        "band, inconsistencies reach 1–4 cycles (about 12–48 cm LOS). "
        "The loop cannot identify the faulty leg or certify the flank "
        "or axial displacement accuracy.",
        "- A triple can contain unwrapping errors that cancel in the loop. "
        "The GUNW comparison (results/gslc_gunw_compare_*.json) and the "
        "Sentinel-1 cross-checks add evidence, but do not exclude every "
        "consistent error.",
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
    ap.add_argument("--cached-only", action="store_true",
                    help="require saved raw unwraps; do not run snaphu")
    ap.add_argument("--regenerate", action="store_true",
                    help="explicitly replace all closure reports and their inventory")
    ap.add_argument("--tags", nargs="*", choices=tuple(PRE_GRANULES),
                    help="partial checks do not replace the combined reports")
    args = ap.parse_args()
    if args.cached_only and (not args.check or args.crop or args.ifg):
        ap.error("--cached-only requires --check without --crop or --ifg")
    if args.regenerate and not args.check:
        ap.error("--regenerate requires --check")
    if args.regenerate and args.tags and set(args.tags) != set(PRE_GRANULES):
        ap.error("--regenerate requires all closure tracks")
    if args.crop:
        crop(tags=args.tags)
    if args.ifg:
        make_ifgs(tags=args.tags)
    if args.check:
        check(tags=args.tags, cached_only=args.cached_only or not args.regenerate,
              regenerate=args.regenerate)
    if not any((args.crop, args.ifg, args.check)):
        ap.print_help()


if __name__ == "__main__":
    main()
