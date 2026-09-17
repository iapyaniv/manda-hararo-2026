"""Validate data files, source manifests, numerical results, and checksums.

Examples::

    python -m mh26.validate --profile core
    python -m mh26.validate --profile full

The core profile covers event, control, and continuation NISAR crops, the
two-track Sentinel-1 inputs, and reported results. The full profile also
requires all 15 pre-event crops and the GSLC unwrap source files (core
checks those only when present). SHA-256 verification is enabled by default.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import numpy as np

from .fetch import read_pairs
from .gunw import corrected_los, load_pair, los_unit_neu, pair_window
from .params import (
    CORE_NISAR_TAGS,
    DATA_NISAR,
    DATA_S1,
    EPSG_UTM37N,
    GRID_POSTING,
    PRE_EVENT_NISAR_TAGS,
    REPO,
    RESULTS,
    S1_FRAME,
    S1_FRAMES,
    S1_PAIR,
    S1_TRACK_PAIRS,
)
from .s1_check import assert_same_grid, coherence_fraction, read_crop


REQUIRED_NPZ_KEYS = {
    "unwrappedPhase",
    "coherence",
    "ionospherePhaseScreen",
    "connectedComponents",
    "mask",
    "x",
    "y",
    "epsg",
    "meta",
    "cube_heightAboveEllipsoid",
    "cube_xCoordinates",
    "cube_yCoordinates",
    "cube_losUnitVectorX",
    "cube_losUnitVectorY",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_checksum_file(path: Path) -> int:
    if not path.exists():
        raise FileNotFoundError(f"checksum manifest missing: {path}")
    count = 0
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match:
            raise ValueError(f"{path}:{line_number}: malformed checksum line")
        expected, relative = match.groups()
        target = path.parent / relative
        if not target.is_file():
            raise FileNotFoundError(f"manifest target missing: {target}")
        actual = sha256(target)
        if actual != expected:
            raise ValueError(f"checksum mismatch: {target} ({actual} != {expected})")
        count += 1
    if count == 0:
        raise ValueError(f"checksum manifest is empty: {path}")
    return count


def validate_npz(tag: str) -> dict:
    data = load_pair(tag)
    if data is None:
        raise FileNotFoundError(DATA_NISAR / f"{tag}.npz")
    missing = REQUIRED_NPZ_KEYS - data.keys()
    if missing:
        raise KeyError(f"{tag}: NPZ keys missing: {sorted(missing)}")
    shape = data["unwrappedPhase"].shape
    for key in ("coherence", "ionospherePhaseScreen", "connectedComponents", "mask"):
        if data[key].shape != shape:
            raise ValueError(f"{tag}: {key} shape {data[key].shape} != {shape}")
    if shape != (len(data["y"]), len(data["x"])):
        raise ValueError(f"{tag}: raster shape does not match coordinates")
    if int(data["epsg"]) != EPSG_UTM37N:
        raise ValueError(f"{tag}: EPSG {int(data['epsg'])} != {EPSG_UTM37N}")
    if not np.allclose(np.diff(data["x"]), GRID_POSTING):
        raise ValueError(f"{tag}: x coordinates are not +{GRID_POSTING:g} m")
    if not np.allclose(np.diff(data["y"]), -GRID_POSTING):
        raise ValueError(f"{tag}: y coordinates are not -{GRID_POSTING:g} m")
    granule_dates = re.findall(r"_(\d{8})T\d{6}", read_pairs()[tag])
    if len(granule_dates) != 4:
        raise ValueError(f"{tag}: could not parse four acquisition times from granule")
    expected_window = (
        f"{granule_dates[0][:4]}-{granule_dates[0][4:6]}-{granule_dates[0][6:]}"
        f"→{granule_dates[2][:4]}-{granule_dates[2][4:6]}-{granule_dates[2][6:]}"
    )
    if pair_window(data) != expected_window:
        raise ValueError(
            f"{tag}: embedded acquisition window {pair_window(data)!r} "
            f"!= manifest {expected_window!r}"
        )
    los, _ = corrected_los(data)
    valid = int(np.isfinite(los).sum())
    if valid < 10_000:
        raise ValueError(f"{tag}: only {valid} valid LOS pixels")
    center = np.array([[float(np.median(data["y"])), float(np.median(data["x"]))]])
    unit = los_unit_neu(data, center)[0]
    if not np.isclose(np.linalg.norm(unit), 1.0, atol=1e-6):
        raise ValueError(f"{tag}: LOS unit vector is not unit length: {unit}")
    return {"shape": shape, "valid": valid, "unit": unit}


def validate_granule_manifests() -> None:
    pairs = read_pairs()
    expected_tags = set(CORE_NISAR_TAGS) | set(PRE_EVENT_NISAR_TAGS)
    if set(pairs) != expected_tags:
        raise ValueError(
            "pairs.tsv tag set differs from configured core+pre-event tags: "
            f"missing={sorted(expected_tags - set(pairs))}, "
            f"extra={sorted(set(pairs) - expected_tags)}"
        )
    granule_lines = [line.strip() for line in
                     (REPO / "granules" / "nisar_gunw.txt").read_text().splitlines()
                     if line.strip()]
    if len(granule_lines) != len(set(granule_lines)):
        raise ValueError("granules/nisar_gunw.txt contains duplicates")
    if set(granule_lines) != set(pairs.values()):
        raise ValueError("nisar_gunw.txt and pairs.tsv list different granules")

    licsar_lines = [
        line.split("\t", 1)
        for line in (REPO / "granules" / "licsar.txt").read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    frame_values = [value for key, value in licsar_lines if key == "frame"]
    pair_values = [tuple(value.split())
                   for key, value in licsar_lines if key == "pair"]
    if sorted(frame_values) != sorted(S1_FRAMES):
        raise ValueError("LiCSAR manifest frames differ from configuration")
    if len(pair_values) != len(set(pair_values)):
        raise ValueError("LiCSAR manifest contains duplicate pair lines")
    if set(pair_values) != set(S1_TRACK_PAIRS):
        raise ValueError("LiCSAR manifest pairs differ from configuration")
    if (S1_FRAME, S1_PAIR, "co-event") not in S1_TRACK_PAIRS:
        raise ValueError("primary cross-check pair missing from S1_TRACK_PAIRS")
    source_names = {
        Path(value).name for key, value in licsar_lines if key == "file"
    }
    expected_names = (
        {f"{frame}.geo.{component}.tif"
         for frame in S1_FRAMES for component in "ENU"}
        | {f"{pair}.geo.{suffix}.tif"
           for _, pair, _ in S1_TRACK_PAIRS for suffix in ("unw", "cc")}
    )
    if source_names != expected_names:
        raise ValueError(
            "LiCSAR manifest files differ from configured frames and pairs: "
            f"missing={sorted(expected_names - source_names)}, "
            f"extra={sorted(source_names - expected_names)}"
        )
    local_names = {path.name for path in DATA_S1.glob("*.tif")}
    if source_names != local_names:
        raise ValueError(
            "LiCSAR manifest differs from local TIFF set: "
            f"missing={sorted(local_names - source_names)}, "
            f"extra={sorted(source_names - local_names)}"
        )

    from .gslc import PAIRS as gslc_pairs

    gslc_rows = [
        line.split("\t")
        for line in (REPO / "granules" / "gslc.tsv").read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if any(len(row) != 3 for row in gslc_rows):
        raise ValueError("granules/gslc.tsv: expected tag<TAB>role<TAB>granule")
    listed = {(tag, role): granule for tag, role, granule in gslc_rows}
    if len(listed) != len(gslc_rows):
        raise ValueError("granules/gslc.tsv contains duplicate tag/role rows")
    from .closure import PRE_GRANULES as closure_pre_granules

    expected_gslc = {
        (tag, role): granule
        for tag, granules in gslc_pairs.items()
        for role, granule in zip(("ref", "sec"), granules)
    }
    expected_gslc.update({
        (tag, "pre"): granule
        for tag, granule in closure_pre_granules.items()
    })
    if listed != expected_gslc:
        raise ValueError("granules/gslc.tsv differs from mh26.gslc PAIRS "
                         "plus mh26.closure PRE_GRANULES")

    expected_geotiffs = {
        f"{tag}.{suffix}.tif"
        for tag in expected_tags
        for suffix in ("coherence", "los_disp_m")
    }
    actual_geotiffs = {
        path.name for path in (REPO / "data" / "geotiff").glob("*.tif")
    }
    if actual_geotiffs != expected_geotiffs:
        raise ValueError(
            "GeoTIFF export set differs from configured NISAR tags: "
            f"missing={sorted(expected_geotiffs - actual_geotiffs)}, "
            f"extra={sorted(actual_geotiffs - expected_geotiffs)}"
        )


def validate_s1() -> None:
    for frame in S1_FRAMES:
        base = []
        for component in "ENU":
            path = DATA_S1 / f"{frame}.geo.{component}.tif"
            if not path.is_file():
                raise FileNotFoundError(path)
            array, x, y = read_crop(path)
            base.append((f"{frame} {component}", array, x, y))
        assert_same_grid(base)
        for pair_frame, pair, _ in S1_TRACK_PAIRS:
            if pair_frame != frame:
                continue
            entries = [base[0]]
            for suffix in ("unw", "cc"):
                path = DATA_S1 / f"{pair}.geo.{suffix}.tif"
                if not path.is_file():
                    raise FileNotFoundError(path)
                array, x, y = read_crop(path)
                entries.append((f"{pair} {suffix}", array, x, y))
            assert_same_grid(entries)
            coherence = coherence_fraction(entries[2][1])
            finite = coherence[np.isfinite(coherence)]
            if finite.min() < 0 or finite.max() > 1:
                raise ValueError(
                    f"{pair}: normalized Sentinel-1 coherence lies outside 0–1"
                )


def validate_results(require_continuation: bool = True) -> None:
    import json

    fit = json.loads((RESULTS / "okada_fit.json").read_text())
    p = np.asarray(fit["values"], dtype=float)
    if len(p) < 7 or not np.isfinite(p).all():
        raise ValueError("okada fit parameters are incomplete or non-finite")
    if not (p[3] > 0 and 0 <= p[4] < p[5] and p[6] >= 0):
        raise ValueError("okada fit geometry is internally invalid")
    expected_volume = p[3] * (p[5] - p[4]) * p[6] / 1e9
    if not np.isclose(fit["volume_km3"], expected_volume, rtol=1e-12):
        raise ValueError("okada fit volume does not match its parameters")
    optimizer = fit["optimizer"]
    if optimizer["status"] <= 0 \
            or optimizer["starts_completed"] != optimizer["starts_attempted"]:
        raise ValueError("okada optimizer did not complete all recorded starts")
    if fit["n_obs"] <= len(p) or fit["rms_m"] < 0:
        raise ValueError("okada fit observation count or RMS is invalid")

    s1 = json.loads((RESULTS / "s1_crosscheck.json").read_text())
    if s1["n_pixels"] <= 0 or s1["rms_m"] < 0 \
            or not np.isfinite([s1["correlation"], s1["amplitude_ratio"],
                                s1["rms_m"]]).all():
        raise ValueError("Sentinel-1 comparison contract is invalid")
    sensitivity = json.loads((RESULTS / "model_sensitivity.json").read_text())
    if not sensitivity["variants"]:
        raise ValueError("model sensitivity table is empty")
    history = json.loads((RESULTS / "probe_history.json").read_text())
    pre_rows = [row for row in history["rows"] if row["phase"] == "pre-event"]
    pre_values = [abs(value) for row in pre_rows
                  for value in row["probe_means_cm"].values()]
    if history["pre_event_windows"] != len(pre_rows) \
            or not np.isclose(history["pre_event_max_abs_probe_cm"],
                              max(pre_values), rtol=1e-12):
        raise ValueError("probe-history summary does not match its rows")
    if require_continuation:
        continuation = json.loads((RESULTS / "continuation.json").read_text())
        interval = continuation["pattern_fraction_cluster_wild_95"]
        envelope = continuation["pattern_fraction_processing_envelope"]
        variants = [row["pattern_fraction"]
                    for row in continuation["processing_sensitivity"]]
        if continuation["n_blocks"] <= 0 or interval[0] > interval[1] \
                or not np.isclose(envelope[0], min(variants)) \
                or not np.isclose(envelope[1], max(variants)):
            raise ValueError("continuation summary does not match its records")
        if continuation["model_fit_sha256"] != sha256(RESULTS / "okada_fit.json"):
            raise ValueError("continuation product predates the current okada fit")

    tracks = json.loads((RESULTS / "s1_tracks.json").read_text())
    expected_pairs = {(frame, pair, role) for frame, pair, role in S1_TRACK_PAIRS}
    recorded_pairs = {(row["frame"], row["pair"], row["role"])
                      for row in tracks["projections"]}
    if recorded_pairs != expected_pairs:
        raise ValueError("Sentinel-1 projections do not match the source manifest")
    if set(tracks["frame_signs"]) != set(S1_FRAMES) \
            or any(sign not in (-1.0, 1.0)
                   for sign in tracks["frame_signs"].values()):
        raise ValueError("Sentinel-1 frame signs are incomplete or invalid")
    for frame, result in tracks["coevent_crosschecks"].items():
        if result["applied_observation_sign"] != tracks["frame_signs"][frame]:
            raise ValueError(f"{frame}: co-event sign does not match frame sign")
    asc = tracks["coevent_crosschecks"][s1["frame"]]
    for key in ("pair", "n_pixels", "applied_observation_sign", "correlation",
                "amplitude_ratio", "rms_m"):
        if asc[key] != s1[key]:
            raise ValueError(f"s1_crosscheck.json disagrees on {key}")
    for row in tracks["projections"]:
        interval = row["pattern_fraction_cluster_wild_95"]
        values = [row["pattern_fraction"], row["zero_intercept_pattern_fraction"],
                  row["correlation"], row["rms_m"], *interval]
        if row["n_blocks"] <= 0 or interval[0] > interval[1] \
                or not np.isfinite(values).all():
            raise ValueError(f"{row['pair']}: invalid projection record")
    nulls = tracks["null_summary"]
    null_values = np.asarray([row["pattern_fraction"]
                              for row in tracks["projections"]
                              if row["role"] == "pre-event"])
    expected_nulls = (len(null_values), float(null_values.mean()),
                      float(null_values.std(ddof=1)),
                      float(np.abs(null_values).max()))
    recorded_nulls = (nulls["n_pairs"], nulls["mean"],
                      nulls["standard_deviation"], nulls["max_abs"])
    if not np.allclose(recorded_nulls, expected_nulls, rtol=1e-12, atol=1e-15):
        raise ValueError("pre-event null summary does not match its projections")
    decay = tracks["post_event_decay"]
    from .s1_tracks import decay_fit
    refit = decay_fit([row for row in tracks["projections"]
                       if row["role"] == "post-event"])
    for key in ("amplitude_total_pattern_fraction", "tau_days", "rms_of_fit"):
        if not np.isclose(decay[key], refit[key], rtol=1e-10, atol=1e-12):
            raise ValueError(f"post-event decay summary disagrees on {key}")


def validate_distributed_result(name: str, expected_source: str) -> None:
    """Check a distributed result against its complete saved solution."""
    import json

    from .distributed import geometry_from_json, geometry_sha256

    payload = json.loads((RESULTS / name).read_text())
    if payload.get("nisar_source") != expected_source:
        raise ValueError(f"{name}: nisar_source is not {expected_source!r}")
    geometry = geometry_from_json(payload["geometry"])
    if payload.get("geometry_sha256") != geometry_sha256(geometry):
        raise ValueError(f"{name}: geometry hash mismatch")
    magnitudes = np.asarray(payload["dislocation_magnitudes_m"], dtype=float)
    if magnitudes.shape != geometry["kinds"].shape \
            or not np.isfinite(magnitudes).all() or (magnitudes < 0).any():
        raise ValueError(f"{name}: invalid dislocation vector")
    dataset_names = [row["name"] for row in payload["datasets"]]
    offsets = payload["dataset_offsets_m"]
    if set(offsets) != set(dataset_names) \
            or not np.isfinite(list(offsets.values())).all() \
            or payload["n_parameters"] != len(magnitudes) + len(dataset_names):
        raise ValueError(f"{name}: parameter records are inconsistent")

    dike = geometry["kinds"] == "dike"
    volume = float(np.sum(geometry["areas"][dike] * magnitudes[dike]))
    if not np.isclose(payload["dike"]["dike_volume_km3"], volume / 1e9,
                      rtol=1e-12):
        raise ValueError(f"{name}: dike volume does not match patches")
    if not np.isclose(payload["dike"]["max_opening_m"],
                      float(magnitudes[dike].max()), rtol=1e-12):
        raise ValueError(f"{name}: max opening does not match patches")
    opening = magnitudes[dike].reshape(geometry["dike_shape"])
    if not np.allclose(payload["opening_distribution"], opening,
                       rtol=1e-12, atol=1e-14):
        raise ValueError(f"{name}: opening grid does not match solution")

    for side in ("east", "west"):
        selected = geometry["kinds"] == f"fault_{side}"
        if not selected.any():
            if side in payload["fault_slip"]:
                raise ValueError(f"{name}: disabled {side} fault has a slip grid")
            continue
        expected = np.full(geometry["fault_shape"], np.nan)
        expected[geometry["fault_masks"][side]] = magnitudes[selected]
        stored = np.asarray([
            [np.nan if value is None else value for value in row]
            for row in payload["fault_slip"][side]
        ])
        if not np.allclose(stored, expected, rtol=1e-12, atol=1e-14,
                           equal_nan=True):
            raise ValueError(f"{name}: {side} slip grid does not match solution")

    fault = ~dike
    fault_potency = float(np.sum(geometry["areas"][fault] * magnitudes[fault]))
    partition = payload["potency_partition"]
    if not np.allclose([partition["dike_m3"], partition["faults_m3"]],
                       [volume, fault_potency], rtol=1e-12):
        raise ValueError(f"{name}: potency partition does not match patches")
    counts = {row["name"]: row["n_blocks"] for row in payload["datasets"]}
    for model_name, values in payload["rms_m"].items():
        expected_all = np.sqrt(sum(counts[key] * values[key] ** 2
                                   for key in counts) / sum(counts.values()))
        if not np.isclose(values["all"], expected_all, rtol=1e-12):
            raise ValueError(f"{name}: {model_name} pooled RMS is inconsistent")
    if name == "distributed_fit.json":
        selection = payload["lambda_selection"]
        sweep = selection["sweep"]
        chosen = selection["chosen_index"]
        scores = [row["cv10_test_rms_mean_m"] for row in sweep]
        if chosen != int(np.argmin(scores)) or chosen in (0, len(sweep) - 1) \
                or not np.isclose(payload["lambda_chosen"],
                                  sweep[chosen]["lambda"]):
            raise ValueError(f"{name}: lambda selection is inconsistent")


def check_three_epoch_report(name: str, expected_source: str,
                             report: dict | None = None) -> None:
    """Check one three-epoch report against the tracked inputs.

    The report never chooses the paths: the delivery directory is the
    tracked one, its files are hashed unconditionally, and the whole
    payload must equal a recomputation from the canonical fit — so a
    stale, edited, or re-pointed report fails.
    """
    import json as _json

    from .distributed import _display_source_path
    from .hamling import DATA_DIR, _find_fit, compute_report

    if report is None:
        report = _json.loads((REPO / "results" / name).read_text())
    if report.get("nisar_source") != expected_source:
        raise ValueError(f"{name}: nisar_source is not {expected_source!r}")
    ref = report["epoch_2026"]
    ref_path = Path(ref["reference_fit"])
    if not ref_path.is_absolute():
        ref_path = REPO / ref_path
    if ref["reference_fit_sha256"] != sha256(ref_path):
        raise ValueError(f"{name}: comparison predates the current "
                         f"{ref['reference_fit']}")
    s_km = np.asarray(report["registration"]["s_km"], dtype=float)
    if s_km.shape != (75,) or not np.all(np.diff(s_km) < 0):
        raise ValueError(f"{name}: registration s axis is not the "
                         "monotone 75-column projection")
    if report["source"]["directory"] != _display_source_path(DATA_DIR):
        raise ValueError(f"{name}: source.directory is not the tracked "
                         f"delivery {_display_source_path(DATA_DIR)}")
    for fname, digest in report["source"]["files"].items():
        path = DATA_DIR / fname
        if not path.is_file():
            raise ValueError(f"{name}: {fname} is not in the tracked "
                             "delivery")
        if sha256(path) != digest:
            raise ValueError(f"{name}: {fname} does not match the "
                             "recorded delivery hash")
    expected, _ = compute_report(expected_source, _find_fit(expected_source),
                                 DATA_DIR)
    expected = _json.loads(_json.dumps(expected, sort_keys=True))
    if expected != report:
        raise ValueError(
            f"{name}: report does not match a full recomputation "
            f"from the canonical fit and the tracked delivery "
            f"(first difference at {_first_json_diff(expected, report)})")


def _first_json_diff(a, b, at=""):
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                return f"{at}.{key} (missing on one side)"
            found = _first_json_diff(a[key], b[key], f"{at}.{key}")
            if found:
                return found
        return None
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return f"{at} (length {len(a)} vs {len(b)})"
        for i, (x, y) in enumerate(zip(a, b)):
            found = _first_json_diff(x, y, f"{at}[{i}]")
            if found:
                return found
        return None
    return at if a != b else None


def check_opening_profiles_report(name: str, expected_source: str) -> None:
    """Check one depth-averaged-profiles report against the tracked inputs.

    mh26.opening_profiles rebuilds its payload through
    mh26.hamling.compute_report, so this check pins the same tracked
    delivery directory and canonical fit and requires exact equality
    with a full recomputation — every guard of the three-epoch check
    applies transitively.
    """
    import json as _json

    from .hamling import DATA_DIR, _find_fit
    from .opening_profiles import compute_report

    report = _json.loads((REPO / "results" / name).read_text())
    if report.get("nisar_source") != expected_source:
        raise ValueError(f"{name}: nisar_source is not {expected_source!r}")
    expected, _ = compute_report(expected_source, _find_fit(expected_source),
                                 DATA_DIR)
    expected = _json.loads(_json.dumps(expected, sort_keys=True))
    if expected != report:
        raise ValueError(
            f"{name}: report does not match a full recomputation "
            f"from the canonical fit and the tracked delivery "
            f"(first difference at {_first_json_diff(expected, report)})")


def validate_new_products(require_gslc_sources: bool = True) -> None:
    """Contract checks on the distributed, decomposition, and GSLC results."""
    import json as _json

    from .distributed import geometry_from_json, geometry_sha256

    for name, expected_source in (
        ("distributed_fit.json", "gunw"),
        ("distributed_fit_gslc.json", "gslc"),
        ("distributed_fit_mapped.json", "gunw"),
        ("distributed_fit_mapped_gslc.json", "gslc"),
        ("distributed_fit_geometry_4a7b2386c8b5.json", "gunw"),
        ("distributed_fit_geometry_4a7b2386c8b5_gslc.json", "gslc"),
    ):
        validate_distributed_result(name, expected_source)

    for dec_name, source in (("decomposition_joint.json", "gunw"),
                             ("decomposition_joint_gslc.json", "gslc")):
        decomposition = _json.loads(
            (REPO / "results" / dec_name).read_text())
        if decomposition.get("nisar_source") != source:
            raise ValueError(f"{dec_name}: nisar_source is not {source!r}")
        for key in ("rift_perpendicular", "zero_north_variant", "free_3d"):
            if key not in decomposition:
                raise ValueError(f"{dec_name}: {key} missing")
            if "closure_probes" not in decomposition[key] \
                    and key != "zero_north_variant":
                raise ValueError(f"{dec_name}: {key} closure_probes missing")
        counts = decomposition["rift_perpendicular"]
        if not (counts["n_cells"] >= counts["n_cells_ge3"] > 0
                and counts["n_cells_ge3"]
                == counts["n_cells_3geo"] + counts["n_cells_4geo"]):
            raise ValueError(f"{dec_name}: cell counts inconsistent")
        free = decomposition["free_3d"]
        if free["n_cells"] != free["n_cells_3geo"] + free["n_cells_4geo"] \
                or free["n_cells"] != counts["n_cells_ge3"]:
            raise ValueError(f"{dec_name}: free-3D cell counts inconsistent")

    from .params import INVERSION_PAIRS

    for tag in INVERSION_PAIRS:
        path = REPO / "results" / f"gslc_gunw_compare_{tag}.json"
        stats = _json.loads(path.read_text())
        if stats["tag"] != tag or not -1 <= stats["correlation"] <= 1 \
                or not np.isfinite(stats["slope_gslc_on_gunw"]):
            raise ValueError(f"{path.name}: comparison identity or fit is invalid")
        for region in ("whole_scene", "near_axis", "off_axis_within_25km"):
            values = stats[region]
            if values["n_cells"] < 0 or values.get("rms_m", 0) < 0 \
                    or values.get("p95_abs_m", 0) < 0:
                raise ValueError(f"{path.name}: {region} summary is invalid")

    unwrap_sensitivity = _json.loads(
        (REPO / "results" / "gslc_unwrap_sensitivity_172_008_A.json").read_text()
    )
    if unwrap_sensitivity["tag"] != "172_008_A":
        raise ValueError("unwrap sensitivity has the wrong pair identity")
    for source_key in ("ifg", "primary_unwrap"):
        source = unwrap_sensitivity[source_key]
        source_path = REPO / source["path"]
        if not source_path.exists() and not require_gslc_sources:
            # The unwrap sources live in the GSLC archive, a full-profile
            # input; a core-profile install still validates the stored
            # contract values without those files.
            continue
        if sha256(source_path) != source["sha256"]:
            raise ValueError(f"unwrap sensitivity {source_key} hash mismatch")
    for variant in unwrap_sensitivity["variants"].values():
        for region in variant["regions"].values():
            fraction = region["fraction_over_half_fringe"]
            if region["n_pixels"] <= 0 or region["rms_difference_m"] < 0 \
                    or region["p95_abs_difference_m"] < 0 \
                    or not 0 <= fraction <= 1:
                raise ValueError("unwrap sensitivity region summary is invalid")

    # Every published result must parse as strict RFC 8259 JSON (a NaN
    # literal in fault_traces.json once broke strict readers).
    def _no_constants(name: str):
        raise ValueError(f"non-finite JSON literal {name!r}")

    for path in sorted((REPO / "results").glob("*.json")):
        _json.loads(path.read_text(), parse_constant=_no_constants)

    # Ridge-anchored fits must carry their complete geometry and trace identity.
    from .distributed import _strike_vectors

    for name in ("distributed_fit_mapped.json",
                 "distributed_fit_mapped_gslc.json"):
        payload = _json.loads((REPO / "results" / name).read_text())
        geom = geometry_from_json(payload["geometry"])
        if payload.get("geometry_sha256") != geometry_sha256(geom):
            raise ValueError(f"{name}: geometry hash mismatch")
        if geom.get("fault_trace_source") != "mapped":
            raise ValueError(f"{name}: fault_trace_source is not 'mapped'")
        for key in ("fault_trace_provenance", "trace_dip_sweep"):
            if key not in payload:
                raise ValueError(f"{name}: {key} missing")
        provenance = payload["fault_trace_provenance"]
        source_path = Path(provenance["path"])
        if not source_path.is_absolute():
            source_path = REPO / source_path
        if source_path.exists() and sha256(source_path) != provenance["sha256"]:
            raise ValueError(f"{name}: trace source hash mismatch")
        if provenance["modeled_fault_sides"] != geom["enabled_fault_sides"]:
            raise ValueError(f"{name}: modeled trace sides disagree")
        supported = [
            side for side in ("east", "west")
            if provenance["trace_support"][side][
                "consistent_axis_side_down_step"
            ]
        ]
        if geom["enabled_fault_sides"] != supported:
            raise ValueError(f"{name}: unsupported trace side is modeled")
        _, perp = _strike_vectors(geom["strike"])
        e_c, n_c = geom["trace_centre"]
        dip_r = np.radians(geom["fault_dip"])
        for side in geom["enabled_fault_sides"]:
            sign = 1.0 if side == "east" else -1.0
            selected = geom["kinds"] == f"fault_{side}"
            patches = geom["patches"][selected]
            offsets = ((patches[:, 0] - n_c) * perp[0]
                       + (patches[:, 1] - e_c) * perp[1])
            deepest = sign * offsets - patches[:, 8] * np.cos(dip_r)
            if not (sign * offsets > 0).all() or deepest.min() < -1.0:
                raise ValueError(f"{name}: {side} fault reaches past the dike plane")
        for side in ("east", "west"):
            stored = int(geom["fault_masks"][side].sum())
            expected = geom["fault_retained_patches"][side]
            if stored != expected:
                raise ValueError(f"{name}: {side} fault mask does not match "
                                 "the retained patch count")

    for name, expected_source in (("fault_depth_mapped.json", "gunw"),
                                  ("fault_depth_mapped_gslc.json", "gslc")):
        sweep = _json.loads((REPO / "results" / name).read_text())
        if sweep.get("nisar_source") != expected_source:
            raise ValueError(f"{name}: nisar_source is not {expected_source!r}")
        ref_name = sweep["reference_fit"]
        if sweep["reference_fit_sha256"] != sha256(REPO / "results" / ref_name):
            raise ValueError(f"{name}: sweep predates the current {ref_name}")
        caps = [row["depth_cap_km"] for row in sweep["rows"]]
        if not caps or caps != sorted(caps):
            raise ValueError(f"{name}: depth caps are not increasing")
        base = sweep["baseline"]
        if not np.isclose(base["fault_moment_Nm"], sweep["shear_modulus_Pa"]
                          * base["fault_potency_m3"]):
            raise ValueError(f"{name}: baseline moment is not mu x potency")

    for name, expected_source in (("uperp_three_epoch.json", "gunw"),
                                  ("uperp_three_epoch_gslc.json", "gslc")):
        check_three_epoch_report(name, expected_source)

    for name, expected_source in (("opening_profiles.json", "gunw"),
                                  ("opening_profiles_gslc.json", "gslc")):
        check_opening_profiles_report(name, expected_source)

    from .export_2026_opening import (VALUE_FORMAT, _find_fit as _export_fit,
                                      compute_export)
    from .hamling import DATA_DIR as _hamling_dir

    for name, expected_source in (("opening_2026_hamling_grid", "gunw"),
                                  ("opening_2026_hamling_grid_gslc", "gslc")):
        report = _json.loads((REPO / "results" / f"{name}.json").read_text())
        if report.get("nisar_source") != expected_source:
            raise ValueError(f"{name}.json: nisar_source is not "
                             f"{expected_source!r}")
        values, expected = compute_export(
            expected_source, _export_fit(expected_source), _hamling_dir)
        expected = _json.loads(_json.dumps(expected, sort_keys=True))
        if expected != report:
            raise ValueError(
                f"{name}.json: report does not match a full recomputation "
                f"(first difference at {_first_json_diff(expected, report)})")
        rendered = "\n".join(VALUE_FORMAT % v
                             for v in values.reshape(-1)) + "\n"
        if (REPO / "results" / f"{name}.txt").read_text() != rendered:
            raise ValueError(f"{name}.txt does not match a full "
                             "recomputation from the canonical fit")

    traces = _json.loads((REPO / "results" / "fault_traces.json").read_text())
    for key in ("los_step_across_ridges", "trace_support",
                "threshold_sensitivity", "cross_geometry_alignment_null"):
        if key not in traces:
            raise ValueError(f"fault_traces.json: {key} missing")

    stress = _json.loads((REPO / "results" / "stress_change.json").read_text())
    source_path = Path(stress["source_fit"]["path"])
    if not source_path.is_absolute():
        source_path = REPO / source_path
    if sha256(source_path) != stress["source_fit"]["sha256"]:
        raise ValueError("stress_change.json: source fit hash mismatch")
    source_fit = _json.loads(source_path.read_text())
    stress_geometry = geometry_from_json(source_fit["geometry"])
    if stress["geometry_sha256"] != geometry_sha256(stress_geometry):
        raise ValueError("stress_change.json: source geometry hash mismatch")
    from .stress import (
        DIKE_STANDOFF_M,
        MAIN_DIP,
        PROFILE_DZ_M,
        fault_patch_receivers,
        fault_profile_receivers,
    )

    base_receivers = fault_patch_receivers(
        stress_geometry, stress_geometry["fault_dip"]
    )
    for side, receivers in base_receivers.items():
        saved = stress_geometry["patches"][
            stress_geometry["kinds"] == f"fault_{side}"
        ]
        if not np.allclose(receivers["points"], saved[:, :3],
                           rtol=0, atol=1e-9):
            raise ValueError(
                f"stress_change.json: {side} base receivers leave saved plane"
            )
    for dip_key, sides in stress["per_dip"].items():
        receivers = fault_patch_receivers(stress_geometry, float(dip_key))
        for side, records in sides.items():
            keep = np.abs(receivers[side]["offsets_m"]) > DIKE_STANDOFF_M
            if records["overall"]["n_receivers"] != int(keep.sum()) \
                    or records["n_excluded_near_dike"] != int((~keep).sum()):
                raise ValueError(
                    f"stress_change.json: dip {dip_key} {side} receiver count mismatch"
                )
            kept_rows = receivers[side]["rows"][keep]
            kept_points = receivers[side]["points"][keep]
            for row in records["rows"]:
                selected = kept_rows == row["row"]
                expected_depth = float(np.mean(kept_points[selected, 2])) / 1e3
                if row["n_receivers"] != int(selected.sum()) \
                        or not np.isclose(row["mean_receiver_depth_km"],
                                          expected_depth, rtol=1e-12):
                    raise ValueError(
                        f"stress_change.json: dip {dip_key} {side} row mismatch"
                    )

    map_receivers = fault_profile_receivers(stress_geometry, MAIN_DIP)
    stress_maps = stress["stress_map_dip65"]
    if set(stress_maps) != set(map_receivers):
        raise ValueError("stress_change.json: spatial-map sides are inconsistent")
    max_depth = float(stress_geometry["fault_z_edges_m"][-1])
    expected_depths = np.arange(PROFILE_DZ_M / 2, max_depth, PROFILE_DZ_M)
    expected_along = np.asarray(stress_geometry["fault_s_centres_m"])
    for side, receivers in map_receivers.items():
        record = stress_maps[side]
        if not np.allclose(record["depth_km"], expected_depths / 1e3) \
                or not np.allclose(record["along_strike_km"],
                                   expected_along / 1e3):
            raise ValueError(f"stress_change.json: {side} spatial axes mismatch")
        keep = np.abs(receivers["offsets_m"]) > DIKE_STANDOFF_M
        shape = (len(expected_depths), len(expected_along))
        expected_offsets = np.full(shape, np.nan)
        rows = np.rint(
            (receivers["depths_m"][keep] - PROFILE_DZ_M / 2) / PROFILE_DZ_M
        ).astype(int)
        columns = receivers["columns"][keep]
        expected_offsets[rows, columns] = receivers["offsets_m"][keep] / 1e3

        def stored_grid(key: str) -> np.ndarray:
            return np.asarray([
                [np.nan if value is None else value for value in row]
                for row in record[key]
            ], dtype=float)

        offsets = stored_grid("horizontal_offset_from_dike_km")
        sigma_n = stored_grid("sigma_n_mpa")
        cfs = stored_grid("cfs_mpa")
        expected_mask = np.isfinite(expected_offsets)
        if offsets.shape != shape or not np.allclose(
                offsets, expected_offsets, rtol=0, atol=1e-12, equal_nan=True
        ):
            raise ValueError(f"stress_change.json: {side} receiver grid mismatch")
        if not np.array_equal(np.isfinite(sigma_n), expected_mask) \
                or not np.array_equal(np.isfinite(cfs), expected_mask) \
                or record["n_receivers"] != int(keep.sum()) \
                or record["n_excluded_near_dike"] != int((~keep).sum()):
            raise ValueError(f"stress_change.json: {side} stress-map mask mismatch")

    unc = _json.loads(
        (REPO / "results" / "distributed_uncertainty.json").read_text())
    for key in ("target_fit", "target_fit_sha256", "geometry_sha256",
                "covariances", "best_fit", "noiseless_refit",
                "interval_method", "intervals_95", "reference_projection"):
        if key not in unc:
            raise ValueError(f"distributed_uncertainty.json: {key} missing")
    uncertainty_target = RESULTS / "distributed_fit.json"
    target_payload = _json.loads(uncertainty_target.read_text())
    if unc["target_fit_sha256"] != sha256(uncertainty_target) \
            or unc["geometry_sha256"] != target_payload["geometry_sha256"]:
        raise ValueError("distributed_uncertainty.json: target fit identity mismatch")
    if unc["n_realizations"] != len(next(iter(unc["draws"].values()))) \
            or any(len(values) != unc["n_realizations"]
                   for values in unc["draws"].values()):
        raise ValueError("distributed_uncertainty.json: draw counts disagree")
    applied = set(unc["reference_projection"]["applied_to"])
    if applied != {"best-model prediction", "every noise draw"}:
        raise ValueError("distributed_uncertainty.json: reference projection incomplete")


def validate_closure(*, require_inputs: bool = False,
                     tags: list[str] | None = None) -> dict:
    """Verify the saved baseline and replay available tracks without writing.

    The published data archives omit most closure inputs. Missing optional
    files are reported, while any present file must match its saved hash.
    ``require_inputs`` requires all selected tracks, or all tracks when
    ``tags`` is omitted. Source, settings, and report checks always apply
    to the whole saved baseline before any raw-array calculation.
    """
    import json
    from . import closure

    inventory = json.loads((RESULTS / "loop_closure_inventory.json").read_text())
    if not isinstance(inventory, dict) or type(inventory.get("schema_version")) is not int \
            or inventory["schema_version"] != 1:
        raise ValueError("loop_closure_inventory.json: unknown schema version")
    text_fields = {"generated_by", "replay_command", "regenerate_command",
                   "scope", "numpy_version"}
    if set(inventory) != text_fields | {"schema_version", "analysis_parameters",
                                       "source_files", "input_files", "reports"} \
            or any(not isinstance(inventory[field], str) or not inventory[field]
                   for field in text_fields):
        raise ValueError("loop_closure_inventory.json: malformed schema")
    if inventory.get("analysis_parameters") != closure.analysis_parameters():
        raise ValueError("loop_closure_inventory.json: analysis settings changed")
    selected = set(tags or closure.PRE_GRANULES)
    if not selected.issubset(closure.PRE_GRANULES):
        raise ValueError("unknown closure track")

    def index_records(records, expected_paths: set[str], label: str) -> dict:
        if not isinstance(records, list):
            raise ValueError(f"closure inventory: malformed {label} records")
        indexed = {}
        for record in records:
            if not isinstance(record, dict) or set(record) != {"path", "sha256", "size_bytes"}:
                raise ValueError(f"closure inventory: malformed {label} record")
            name = record["path"]
            if not isinstance(name, str) or not name:
                raise ValueError("closure inventory: invalid path")
            relative = Path(name)
            if relative.is_absolute() or ".." in relative.parts \
                    or "\\" in name or relative.as_posix() != name:
                raise ValueError("closure inventory paths must be safe repository-relative paths")
            digest = record["sha256"]
            if not isinstance(digest, str) or len(digest) != 64 \
                    or any(char not in "0123456789abcdef" for char in digest) \
                    or type(record["size_bytes"]) is not int or record["size_bytes"] < 0:
                raise ValueError(f"closure inventory: malformed identity for {name}")
            if name in indexed:
                raise ValueError(f"closure inventory: duplicate {label} path {name}")
            indexed[name] = record
        if set(indexed) != expected_paths:
            raise ValueError(f"closure inventory: missing or unexpected {label} paths")
        return indexed

    def verify(record: dict, *, optional: bool = False) -> bool:
        relative = Path(record["path"])
        path = REPO / relative
        if optional and not path.exists():
            return False
        if path.stat().st_size != record["size_bytes"] \
                or sha256(path) != record["sha256"]:
            raise ValueError(f"{relative}: closure content hash mismatch")
        return True

    expected_sources = {"mh26/closure.py", "mh26/gslc.py", "mh26/params.py",
                        "results/okada_fit.json"}
    expected_reports = {f"results/loop_closure_{tag}.json"
                        for tag in closure.PRE_GRANULES}
    expected_reports |= {"results/loop_closure.md", "results/loop_closure_nisar.png"}
    expected_inputs = {
        str(path.relative_to(REPO)) for tag in closure.PRE_GRANULES for leg in closure.LEGS
        for path in (closure._leg_file(tag, leg),
                     closure.DATA_GSLC / f"{tag}_loop{leg}_unwraw.npz")}
    sources = index_records(inventory.get("source_files"), expected_sources, "source")
    reports = index_records(inventory.get("reports"), expected_reports, "report")
    recorded = index_records(inventory.get("input_files"), expected_inputs, "input")
    for record in [*sources.values(), *reports.values()]:
        verify(record)

    saved = {}
    input_records = {}
    track_inputs = {}
    for tag in closure.PRE_GRANULES:
        stats = json.loads((RESULTS / f"loop_closure_{tag}.json").read_text())
        saved[tag] = stats
        if set(stats["input_files"]) != set(closure.LEGS):
            raise ValueError(f"{tag}: incomplete closure leg inventory")
        track_inputs[tag] = []
        for leg in closure.LEGS:
            records = stats["input_files"][leg]
            expected_paths = {
                "interferogram": str(closure._leg_file(tag, leg).relative_to(REPO)),
                "raw_unwrap_cache": str((closure.DATA_GSLC / f"{tag}_loop{leg}_unwraw.npz")
                                        .relative_to(REPO)),
            }
            if not isinstance(records, dict) or set(records) != set(expected_paths):
                raise ValueError(f"{tag} {leg}: incorrect closure input inventory")
            index_records(list(records.values()), set(expected_paths.values()), f"{tag} {leg}")
            if any(records[role]["path"] != path for role, path in expected_paths.items()):
                raise ValueError(f"{tag} {leg}: incorrect closure input roles")
            configuration = stats["unwrap_configuration"][leg]
            source = stats["unwrap_configuration_source"][leg]
            if configuration is None:
                if source != "unknown: legacy cache has no config metadata":
                    raise ValueError(f"{tag} {leg}: unknown settings lack an "
                                     "explicit provenance statement")
            elif not isinstance(configuration, dict) or source != "cache config metadata":
                raise ValueError(f"{tag} {leg}: unsupported settings provenance")
            execution = stats["snaphu_execution"][leg]
            if execution is not None:
                if not isinstance(execution, dict) or set(execution) != {
                        "resolved_path", "sha256", "reported_version", "version_exit_code"} \
                        or not isinstance(execution["resolved_path"], str) \
                        or not Path(execution["resolved_path"]).is_absolute() \
                        or not isinstance(execution["sha256"], str) \
                        or len(execution["sha256"]) != 64 \
                        or any(char not in "0123456789abcdef" for char in execution["sha256"]) \
                        or type(execution["version_exit_code"]) is not int \
                        or (execution["reported_version"] is not None
                            and not isinstance(execution["reported_version"], str)):
                    raise ValueError(f"{tag} {leg}: malformed snaphu executable identity")
            error = stats["leg_congruence_max_abs_rad"][leg]
            if not np.isfinite(error) or error > closure.CONGRUENCE_TOLERANCE_RAD:
                raise ValueError(f"{tag} {leg}: phase congruence check failed")
            for record in records.values():
                input_records[record["path"]] = record
                track_inputs[tag].append(record["path"])

    if recorded != input_records:
        raise ValueError("closure inventory and per-track input identities disagree")
    checksums = "\n".join(sorted(
        f"{record['sha256']}  {Path(record['path']).relative_to('data')}"
        for record in input_records.values())) + "\n"
    if (REPO / "data" / "SHA256SUMS-closure.txt").read_text() != checksums:
        raise ValueError("SHA256SUMS-closure.txt differs from the closure inventory")
    missing = [path for path, record in input_records.items()
               if not verify(record, optional=True)]
    required_missing = set(missing).intersection(
        path for tag in selected for path in track_inputs[tag])
    if required_missing and require_inputs:
        raise FileNotFoundError(f"{len(required_missing)} optional closure inputs missing; "
                                "see results/loop_closure_inventory.json")
    available = [tag for tag, paths in track_inputs.items()
                 if tag in selected and not set(paths).intersection(missing)]
    recomputed = {}
    if available:
        recomputed = closure._compute(tags=available, cached_only=True,
                                      write_outputs=False)

        def same_statistics(actual, expected) -> bool:
            if isinstance(expected, dict):
                return (isinstance(actual, dict) and actual.keys() == expected.keys()
                        and all(same_statistics(actual[key], value)
                                for key, value in expected.items()))
            if isinstance(expected, float):
                return (isinstance(actual, float)
                        and bool(np.isclose(actual, expected, rtol=1e-10, atol=1e-12)))
            if isinstance(expected, int):
                return isinstance(actual, int) and actual == expected
            return actual == expected

        derived = {"error_region", "congruence_max_abs_rad",
                   "leg_congruence_max_abs_rad", "whole_scene", "near_axis",
                   "off_axis_within_25km"}
        for tag, stats in recomputed.items():
            exact = {key: value for key, value in stats.items() if key not in derived}
            saved_exact = {key: value for key, value in saved[tag].items()
                           if key not in derived}
            if exact != saved_exact or any(
                    not same_statistics(stats[key], saved[tag][key])
                    for key in derived):
                raise ValueError(f"loop_closure_{tag}.json: saved values differ "
                                 "from raw-array reproduction")
    return {"n_inputs": len(input_records),
            "n_inputs_checked": len(input_records) - len(missing),
            "reproduced_tracks": available, "missing_inputs": missing,
            "recomputed": recomputed}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("core", "full"), default="core")
    parser.add_argument("--skip-checksums", action="store_true")
    parser.add_argument("--require-closure-inputs", action="store_true",
                        help="fail unless optional raw closure inputs are "
                             "present and all three tracks can be replayed")
    args = parser.parse_args()

    validate_granule_manifests()
    print("PASS NISAR and LiCSAR source manifests agree")
    tags = list(CORE_NISAR_TAGS)
    if args.profile == "full":
        tags += list(PRE_EVENT_NISAR_TAGS)
    for tag in tags:
        info = validate_npz(tag)
        print(f"PASS {tag}: {info['shape']}, {info['valid']} valid pixels")
    validate_s1()
    print("PASS Sentinel-1 grids and coherence encoding")
    validate_results(require_continuation=True)
    print("PASS result contracts and derived-summary consistency")
    validate_new_products(require_gslc_sources=args.profile == "full")
    print("PASS distributed, decomposition, GSLC-comparison, ridge/stress/"
          "uncertainty contracts and strict-JSON results")
    closure = validate_closure(require_inputs=args.require_closure_inputs)
    print("PASS closure report/source hashes and "
          f"{closure['n_inputs_checked']}/{closure['n_inputs']} input hashes")
    if closure["reproduced_tracks"]:
        print(f"PASS closure raw-array reproduction: "
              f"{len(closure['reproduced_tracks'])}/3 tracks")
    if closure["missing_inputs"]:
        print("SKIP remaining closure reproduction: optional raw inputs "
              "are missing. The full archive profile alone does not "
              "reproduce closure; use --require-closure-inputs to require it.")
    if not args.skip_checksums:
        # The s1tracks and gslc manifests stay separate so the 2026-08-18
        # core archive verifies unchanged next to the later add-ons.
        if args.profile == "full":
            manifest_counts = (
                ("SHA256SUMS.txt", 104),
                ("SHA256SUMS-gslc.txt", 4),
            )
        else:
            manifest_counts = (
                ("SHA256SUMS-core.txt", 50),
                ("SHA256SUMS-s1tracks.txt", 35),
            )
        for name, expected_count in manifest_counts:
            manifest = REPO / "data" / name
            count = verify_checksum_file(manifest)
            if count != expected_count:
                raise ValueError(
                    f"{name} has {count} entries; expected {expected_count}"
                )
            print(f"PASS {count} SHA-256 checksums from {name}")
    print(f"VALIDATION PASSED ({args.profile} profile)")


if __name__ == "__main__":
    main()
