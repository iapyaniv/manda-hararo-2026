"""Validate data files, source manifests, numerical results, and checksums.

Examples::

    python -m mh26.validate --profile core
    python -m mh26.validate --profile full

The core profile covers event, control, and continuation NISAR crops,
Sentinel-1 inputs, and reported results. The full profile also requires all 15
pre-event crops. SHA-256 verification is enabled by default.
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
    S1_PAIR,
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
    pair_values = [value for key, value in licsar_lines if key == "pair"]
    if frame_values != [S1_FRAME] or pair_values != [S1_PAIR]:
        raise ValueError("LiCSAR manifest frame/pair differs from configuration")
    source_names = {
        Path(value).name for key, value in licsar_lines if key == "file"
    }
    local_names = {path.name for path in DATA_S1.glob("*.tif")}
    if source_names != local_names:
        raise ValueError(
            "LiCSAR manifest differs from local TIFF set: "
            f"missing={sorted(local_names - source_names)}, "
            f"extra={sorted(source_names - local_names)}"
        )

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
    inputs = []
    for name, path in (
        ("unwrapped", DATA_S1 / f"{S1_PAIR}.geo.unw.tif"),
        ("coherence", DATA_S1 / f"{S1_PAIR}.geo.cc.tif"),
        ("E", DATA_S1 / f"{S1_FRAME}.geo.E.tif"),
        ("N", DATA_S1 / f"{S1_FRAME}.geo.N.tif"),
        ("U", DATA_S1 / f"{S1_FRAME}.geo.U.tif"),
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
        array, x, y = read_crop(path)
        inputs.append((name, array, x, y))
    assert_same_grid(inputs)
    coherence = coherence_fraction(inputs[1][1])
    finite = coherence[np.isfinite(coherence)]
    if finite.min() < 0 or finite.max() > 1:
        raise ValueError("normalized Sentinel-1 coherence lies outside 0–1")


def validate_results(require_continuation: bool = True) -> None:
    import json

    fit = json.loads((RESULTS / "okada_fit.json").read_text())
    p = fit["values"]
    checks = (
        (10_000 < p[3] < 30_000, "dike length outside sanity range"),
        (0 <= p[4] < p[5] < 12_000, "dike depth bounds invalid"),
        (0.5 < p[6] < 5.0, "dike opening outside sanity range"),
        (0.1 < fit["volume_km3"] < 0.3, "dike volume outside sanity range"),
        (fit["rms_m"] < 0.10, "NISAR model RMS exceeds 10 cm"),
        (fit["optimizer"]["optimality"] < 0.1, "optimizer optimality is too large"),
    )
    for passed, message in checks:
        if not passed:
            raise ValueError(message)

    s1 = json.loads((RESULTS / "s1_crosscheck.json").read_text())
    if s1["correlation"] < 0.90 or s1["rms_m"] > 0.05:
        raise ValueError("Sentinel-1 comparison is outside the defined result bounds")
    sensitivity = json.loads((RESULTS / "model_sensitivity.json").read_text())
    if len(sensitivity["variants"]) < 5:
        raise ValueError("model sensitivity table is unexpectedly incomplete")
    history = json.loads((RESULTS / "probe_history.json").read_text())
    if history["pre_event_windows"] != 16:
        raise ValueError("probe history does not contain 16 pre-event windows")
    if history["pre_event_max_abs_probe_cm"] > 5.0:
        raise ValueError("pre-event probe history exceeds the 5 cm sanity bound")
    if require_continuation:
        continuation = json.loads((RESULTS / "continuation.json").read_text())
        if continuation["n_blocks"] < 100:
            raise ValueError("continuation test has too few blocks")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("core", "full"), default="core")
    parser.add_argument("--skip-checksums", action="store_true")
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
    print("PASS reported result bounds")
    if not args.skip_checksums:
        manifest = REPO / "data" / (
            "SHA256SUMS.txt" if args.profile == "full" else "SHA256SUMS-core.txt"
        )
        count = verify_checksum_file(manifest)
        expected_count = 65 if args.profile == "full" else 50
        if count != expected_count:
            raise ValueError(
                f"{manifest.name} has {count} entries; expected {expected_count}"
            )
        print(f"PASS {count} SHA-256 checksums from {manifest.name}")
    print(f"VALIDATION PASSED ({args.profile} profile)")


if __name__ == "__main__":
    main()
