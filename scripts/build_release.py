#!/usr/bin/env python3
"""Build deterministic, checksum-verified code and data archives.

Release workflow:

1. Run ``python scripts/build_release.py --manifests-only`` and commit the
   resulting data manifests with all code/results changes.
2. From that clean commit, run ``python scripts/build_release.py --archives-only``.

The core and pre-event archives have separate checksum manifests, so each
archive can be verified without the other.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mh26.params import (  # noqa: E402
    CORE_NISAR_TAGS,
    DATA_NISAR,
    DATA_S1,
    PRE_EVENT_NISAR_TAGS,
    S1_FRAME,
    S1_FRAMES,
    S1_PAIR,
    S1_TRACK_PAIRS,
)
from mh26.gslc import DATA_GSLC  # noqa: E402


DATA = ROOT / "data"
GEOTIFF = DATA / "geotiff"
DIST = ROOT / "dist"
CODE_ARCHIVE = DIST / "manda-hararo-2026_code.zip"
CORE_ARCHIVE = DIST / "manda-hararo-2026_data-core.zip"
PREEVENT_ARCHIVE = DIST / "manda-hararo-2026_data-preevent.zip"
S1TRACKS_ARCHIVE = DIST / "manda-hararo-2026_data-s1tracks.zip"
GSLC_ARCHIVE = DIST / "manda-hararo-2026_data-gslc.zip"


def run_git(*args: str) -> str:
    result = subprocess.run(
        ("git", *args), cwd=ROOT, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_files(paths: list[Path]) -> None:
    missing = [path for path in paths if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  {path}" for path in missing)
        raise SystemExit(f"release inputs missing:\n{formatted}")


def core_files() -> list[Path]:
    files = [DATA_NISAR / f"{tag}.npz" for tag in CORE_NISAR_TAGS]
    files += [
        DATA_S1 / f"{S1_PAIR}.geo.unw.tif",
        DATA_S1 / f"{S1_PAIR}.geo.cc.tif",
        DATA_S1 / f"{S1_FRAME}.geo.E.tif",
        DATA_S1 / f"{S1_FRAME}.geo.N.tif",
        DATA_S1 / f"{S1_FRAME}.geo.U.tif",
    ]
    # The core archive includes all GeoTIFFs. The pre-event NPZ files remain in
    # the separate add-on because they account for most of its size.
    tags = CORE_NISAR_TAGS + PRE_EVENT_NISAR_TAGS
    files += [
        GEOTIFF / f"{tag}.{suffix}.tif"
        for tag in tags
        for suffix in ("coherence", "los_disp_m")
    ]
    return files


def preevent_files() -> list[Path]:
    return [DATA_NISAR / f"{tag}.npz" for tag in PRE_EVENT_NISAR_TAGS]


def s1tracks_files() -> list[Path]:
    """Two-track LiCSAR add-on: everything not already in the core archive."""
    in_core = {
        f"{S1_PAIR}.geo.unw.tif",
        f"{S1_PAIR}.geo.cc.tif",
    } | {f"{S1_FRAME}.geo.{component}.tif" for component in "ENU"}
    names = {
        f"{frame}.geo.{component}.tif"
        for frame in S1_FRAMES for component in "ENU"
    } | {
        f"{pair}.geo.{suffix}.tif"
        for _, pair, _ in S1_TRACK_PAIRS for suffix in ("unw", "cc")
    }
    return [DATA_S1 / name for name in sorted(names - in_core)]


def gslc_files() -> list[Path]:
    """Derived 30-m interferograms and unwrapped fields for the two inversion
    pairs; other mh26.gslc pairs and the 1.6-GB GSLC crops are regenerable
    from ASF and are not shipped."""
    from mh26.params import INVERSION_PAIRS

    return [DATA_GSLC / f"{tag}_{kind}.npz"
            for tag in INVERSION_PAIRS for kind in ("ifg", "unw")]


def manifest_text(paths: list[Path]) -> str:
    lines = [f"{sha256(path)}  {path.relative_to(DATA).as_posix()}" for path in paths]
    return "\n".join(sorted(lines)) + "\n"


def write_if_changed(path: Path, content: str) -> None:
    if path.exists() and path.read_text() == content:
        return
    path.write_text(content)
    print(f"wrote {path.relative_to(ROOT)}")


def write_manifests() -> None:
    core = core_files()
    preevent = preevent_files()
    s1tracks = s1tracks_files()
    gslc = gslc_files()
    require_files(core + preevent + s1tracks + gslc)
    write_if_changed(DATA / "SHA256SUMS-core.txt", manifest_text(core))
    write_if_changed(DATA / "SHA256SUMS-preevent.txt", manifest_text(preevent))
    write_if_changed(DATA / "SHA256SUMS-s1tracks.txt", manifest_text(s1tracks))
    write_if_changed(DATA / "SHA256SUMS-gslc.txt", manifest_text(gslc))
    write_if_changed(DATA / "SHA256SUMS.txt",
                     manifest_text(core + preevent + s1tracks + gslc))


def commit_datetime() -> tuple[int, int, int, int, int, int]:
    stamp = int(run_git("show", "-s", "--format=%ct", "HEAD"))
    value = datetime.fromtimestamp(stamp, tz=timezone.utc)
    # ZIP timestamps have two-second precision.
    return value.year, value.month, value.day, value.hour, value.minute, value.second // 2 * 2


def tracked_files() -> list[Path]:
    output = subprocess.run(
        ("git", "ls-files", "-z"), cwd=ROOT, check=True,
        stdout=subprocess.PIPE,
    ).stdout
    paths = [ROOT / item.decode() for item in output.split(b"\0") if item]
    require_files(paths)
    return paths


def add_file(archive: zipfile.ZipFile, source: Path, name: str,
             timestamp: tuple[int, int, int, int, int, int]) -> None:
    info = zipfile.ZipInfo(name, date_time=timestamp)
    info.compress_type = zipfile.ZIP_STORED
    info.create_system = 3
    mode = 0o100755 if os.access(source, os.X_OK) else 0o100644
    info.external_attr = mode << 16
    with source.open("rb") as reader, archive.open(info, "w") as writer:
        shutil.copyfileobj(reader, writer, length=8 * 1024 * 1024)


def write_archive(output: Path, entries: list[tuple[str, Path]], comment: str,
                  timestamp: tuple[int, int, int, int, int, int]) -> None:
    require_files([source for _, source in entries])
    output.parent.mkdir(exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    with zipfile.ZipFile(temporary, "w", allowZip64=True) as archive:
        archive.comment = comment.encode()
        for name, source in sorted(entries):
            add_file(archive, source, name, timestamp)
    with zipfile.ZipFile(temporary) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise RuntimeError(f"CRC verification failed for {bad} in {temporary}")
    temporary.replace(output)
    print(f"wrote {output.relative_to(ROOT)} ({output.stat().st_size / 1e6:.1f} MB)")


def build_archives() -> None:
    dirty = run_git("status", "--porcelain", "--untracked-files=no")
    if dirty:
        raise SystemExit(
            "tracked worktree is dirty; commit the manifest changes before "
            "building archives"
        )
    commit = run_git("rev-parse", "HEAD")
    timestamp = commit_datetime()

    code_entries = [
        (f"manda-hararo-2026/{path.relative_to(ROOT).as_posix()}", path)
        for path in tracked_files()
    ]
    write_archive(CODE_ARCHIVE, code_entries, commit, timestamp)

    core = core_files()
    core_entries = [(path.relative_to(DATA).as_posix(), path) for path in core]
    core_entries += [
        ("README.md", DATA / "README.md"),
        ("SHA256SUMS-core.txt", DATA / "SHA256SUMS-core.txt"),
    ]
    write_archive(CORE_ARCHIVE, core_entries, commit, timestamp)

    preevent = preevent_files()
    preevent_entries = [(path.relative_to(DATA).as_posix(), path) for path in preevent]
    preevent_entries += [
        ("README.md", DATA / "README.md"),
        ("SHA256SUMS-preevent.txt", DATA / "SHA256SUMS-preevent.txt"),
    ]
    write_archive(PREEVENT_ARCHIVE, preevent_entries, commit, timestamp)

    s1tracks = s1tracks_files()
    s1tracks_entries = [(path.relative_to(DATA).as_posix(), path)
                        for path in s1tracks]
    s1tracks_entries += [
        ("README.md", DATA / "README.md"),
        ("SHA256SUMS-s1tracks.txt", DATA / "SHA256SUMS-s1tracks.txt"),
    ]
    write_archive(S1TRACKS_ARCHIVE, s1tracks_entries, commit, timestamp)

    gslc = gslc_files()
    gslc_entries = [(path.relative_to(DATA).as_posix(), path) for path in gslc]
    gslc_entries += [
        ("README.md", DATA / "README.md"),
        ("SHA256SUMS-gslc.txt", DATA / "SHA256SUMS-gslc.txt"),
    ]
    write_archive(GSLC_ARCHIVE, gslc_entries, commit, timestamp)

    write_if_changed(
        DIST / "SHA256SUMS.txt",
        "".join(
            f"{sha256(path)}  {path.name}\n"
            for path in (CODE_ARCHIVE, CORE_ARCHIVE, PREEVENT_ARCHIVE,
                         S1TRACKS_ARCHIVE, GSLC_ARCHIVE)
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--manifests-only", action="store_true")
    group.add_argument("--archives-only", action="store_true")
    args = parser.parse_args()
    if not args.archives_only:
        write_manifests()
    if not args.manifests_only:
        build_archives()


if __name__ == "__main__":
    main()
