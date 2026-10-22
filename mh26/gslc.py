"""full-resolution interferograms from the 10 m NISAR GSLC products.

The GUNW analysis chain works at the product's 80 m unwrapped posting and
inherits the product's unwrapping (subswath-edge artifacts required the
valid-sample mask). This module rebuilds the two co-event interferograms from
the L2 GSLC products instead: GSLC epochs of one track share the same UTM
grid and are DEM-flattened, so the interferogram is a per-pixel product with
no co-registration or topographic step. We control multilooking, filtering,
and unwrapping (snaphu). The GSLC result keeps near-axis coherent cells that
the GUNW result omits, but its fringe count remains processing-sensitive.

Stages (each writes files the next stage reads):

    python -m mh26.gslc --probe            # print the remote HDF5 layout, no data
    python -m mh26.gslc --crop [--download]  # extract 90x90 km HH crops (data/gslc/)
    python -m mh26.gslc --ifg              # form, multilook, filter, coherence
    python -m mh26.gslc --unwrap           # snaphu, iono reuse, referencing

`--crop` streams ranged reads through earthaccess/fsspec by default;
`--download` falls back to whole-granule download-extract-delete (the GSLC
granules are 6-12 GB each). Earthdata Login: a user token in
EARTHDATA_TOKEN or ~/.config/earthdata/token (one line, mode 600) is used
when present, otherwise ~/.netrc.

The ionosphere screen cannot be re-derived from single-band GSLC data. The
GUNW `ionospherePhaseScreen` for the same acquisition pair (data/nisar/) is
therefore interpolated onto the interferogram grid and subtracted. Both
products use the same L-band carrier, so no wavelength conversion is needed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .gunw import load_pair
from .params import CENTER_E, CENTER_N, FARFIELD_RADIUS, REPO, RESULTS, WAVELENGTH_L

DATA_GSLC = REPO / "data" / "gslc"
TMP = REPO / "data" / "h5tmp"
CROP_HALF = 45_000.0            # m, same window as the GUNW crops
# The provisional GSLC grid is 10 m in x and 5 m in y, so square multilooked
# pixels need twice the looks in y: (6, 3) looks -> 30 x 30 m posting.
LOOKS = 3                       # x looks; y looks = 2x this
GOLDSTEIN_ALPHA = 0.6
GOLDSTEIN_TILE = 64
COHERENCE_WINDOW = 3            # multilooked pixels (90 m estimation window)
UNWRAP_COHERENCE_MIN = 0.25     # snaphu sees everything; analysis mask only
SNAPHU = "snaphu"
SNAPHU_PRIMARY = {"stat_cost_mode": "DEFO", "init_method": "MCF"}
SNAPHU_SENSITIVITY_VARIANTS = {
    "smooth_mcf": {"stat_cost_mode": "SMOOTH", "init_method": "MCF"},
}

GRIDS = "/science/LSAR/GSLC/grids/frequencyA"
IDENT = "/science/LSAR/identification"
EARTHDATA_TOKEN_FILE = Path.home() / ".config" / "earthdata" / "token"

# tag -> (reference GSLC granule, secondary GSLC granule). The tag names the
# GUNW crop of the same pair, which supplies the ionosphere screen. The two
# inversion pairs ship in the data-gslc release archive; other pairs rebuild
# on demand.
PAIRS = {
    "172_008_A": (
        "NISAR_L2_PR_GSLC_024_172_A_008_2005_DHDH_A_20260708T024619_20260708T024654_P05023_N_F_J_001",
        "NISAR_L2_PR_GSLC_025_172_A_008_2005_DHDH_A_20260720T024618_20260720T024653_P05023_N_F_J_001",
    ),
    "165_082_D": (
        "NISAR_L2_PR_GSLC_024_165_D_082_2005_SHSH_A_20260707T154911_20260707T154946_P05023_N_F_J_001",
        "NISAR_L2_PR_GSLC_025_165_D_082_2005_DHDH_A_20260719T154911_20260719T154945_P05023_N_F_J_001",
    ),
    "100_A_008_STRADDLE": (
        "NISAR_L2_PR_GSLC_025_100_A_008_2005_DHDH_A_20260715T025438_20260715T025512_P05023_N_F_J_001",
        "NISAR_L2_PR_GSLC_026_100_A_008_2005_DHDH_A_20260727T025437_20260727T025512_P05023_N_F_J_001",
    ),
}
COLLECTION = "NISAR_L2_GSLC_PROVISIONAL_V1"


# --- Earthdata Login ----------------------------------------------------------

def earthdata_login():
    """earthaccess login, with an Earthdata user token (EARTHDATA_TOKEN or
    ~/.config/earthdata/token) preferred over the ~/.netrc password."""
    import earthaccess

    token = os.environ.get("EARTHDATA_TOKEN", "").strip()
    if not token and EARTHDATA_TOKEN_FILE.is_file():
        if EARTHDATA_TOKEN_FILE.stat().st_mode & 0o077:
            raise PermissionError(f"{EARTHDATA_TOKEN_FILE} must have mode 600")
        token = EARTHDATA_TOKEN_FILE.read_text().strip()
    if token:
        os.environ["EARTHDATA_TOKEN"] = token
        return earthaccess.login(strategy="environment")
    return earthaccess.login(strategy="netrc")


# --- stage 0: layout probe ----------------------------------------------------

def _remote_file(granule: str):
    import earthaccess
    import h5py

    earthdata_login()
    hits = earthaccess.search_data(short_name=COLLECTION, granule_name=granule)
    if not hits:
        raise SystemExit(f"granule not in archive: {granule}")
    url = [u for u in hits[0].data_links() if u.endswith(".h5")][0]
    fs = earthaccess.get_fsspec_https_session()
    fobj = fs.open(url, mode="rb", cache_type="blockcache", block_size=32 * 1024 * 1024)
    return h5py.File(fobj, "r", rdcc_nbytes=64 * 1024 * 1024)


def probe(granule: str) -> None:
    with _remote_file(granule) as f:
        def walk(name, obj):
            import h5py

            if isinstance(obj, h5py.Dataset):
                print(f"{name}  {obj.shape} {obj.dtype} chunks={obj.chunks}")
        f[GRIDS].visititems(lambda n, o: walk(f"{GRIDS}/{n}", o))
        print("identification:")
        for k in sorted(f[IDENT]):
            v = f[IDENT][k][()]
            print(f"  {k} = {v.decode() if isinstance(v, bytes) else v}")


# --- stage 1: crops -----------------------------------------------------------

def _detect_pol(f) -> str:
    for pol in ("HH", "VV"):
        if pol in f[GRIDS]:
            return pol
    raise KeyError(f"no co-pol dataset under {GRIDS}: {list(f[GRIDS].keys())}")


ANALYSIS_WINDOW = (CENTER_E - CROP_HALF, CENTER_E + CROP_HALF,
                   CENTER_N - CROP_HALF, CENTER_N + CROP_HALF)


def _extract_crop(f, granule: str, out: Path,
                  window: tuple[float, float, float, float] = ANALYSIS_WINDOW) -> None:
    """Crop the co-pol layer over `window` = (x0, x1, y0, y1) in product UTM metres."""
    pol = _detect_pol(f)
    g = f[GRIDS]
    x = g["xCoordinates"][:]
    y = g["yCoordinates"][:]
    x0, x1, y0, y1 = window
    ix = np.where((x >= x0) & (x <= x1))[0]
    iy = np.where((y >= y0) & (y <= y1))[0]
    if len(ix) == 0 or len(iy) == 0:
        raise SystemExit(f"{granule}: crop window outside the product grid")
    ds = g[pol]
    # Read chunk-aligned row bands to keep ranged reads sequential.
    rows = []
    step = max(1, (ds.chunks[0] if ds.chunks else 512) * 8)
    for r0 in range(iy[0], iy[-1] + 1, step):
        r1 = min(r0 + step, iy[-1] + 1)
        rows.append(ds[r0:r1, ix[0]:ix[-1] + 1])
        print(f"  rows {r0}-{r1} of {iy[0]}-{iy[-1] + 1}", flush=True)
    z = np.concatenate(rows, axis=0).astype(np.complex64)

    # Subswath/validity rasters: the provisional GUNW subswath-edge artifacts
    # motivated the valid-sample mask there; keep the GSLC equivalents so the
    # interferogram stage can screen the same defect at full resolution.
    aux = {}
    for name in ("mask", "inputDataExceptionMask"):
        if name in g:
            aux[name] = g[name][iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]

    meta = {}
    for k in ("orbitPassDirection", "zeroDopplerStartTime", "trackNumber",
              "frameNumber"):
        if k in f[IDENT]:
            v = f[IDENT][k][()]
            meta[k] = v.decode() if isinstance(v, bytes) else (
                v.item() if isinstance(v, np.generic) else v)
    epsg = None
    for key in ("projection", "epsg"):
        if key in g:
            epsg = int(np.array(g[key]).item())
            break
    if epsg is None:
        raise KeyError(f"{granule}: projection dataset missing")

    DATA_GSLC.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    with tmp.open("wb") as stream:
        np.savez(stream, z=z, x=x[ix[0]:ix[-1] + 1], y=y[iy[0]:iy[-1] + 1],
                 epsg=np.array(epsg), pol=np.array(pol),
                 granule=np.array(granule),
                 meta=np.array(json.dumps(meta, sort_keys=True)), **aux)
    tmp.replace(out)
    print(f"  saved {out.name} ({out.stat().st_size / 1e6:.0f} MB, "
          f"{z.shape[0]} x {z.shape[1]})")


def crop(download: bool = False, tags: list[str] | None = None) -> None:
    import earthaccess

    earthdata_login()
    for tag, granules in PAIRS.items():
        if tags and tag not in tags:
            continue
        for role, granule in zip(("ref", "sec"), granules):
            out = DATA_GSLC / f"{tag}_{role}.npz"
            if out.exists():
                print(f"{tag} {role}: already extracted")
                continue
            print(f"{tag} {role}: {granule}")
            if download:
                import h5py

                hits = earthaccess.search_data(short_name=COLLECTION,
                                               granule_name=granule)
                if not hits:
                    raise SystemExit(f"granule not in archive: {granule}")
                TMP.mkdir(parents=True, exist_ok=True)
                files = earthaccess.download(hits[:1], str(TMP))
                h5path = Path(files[0])
                with h5py.File(h5path, "r") as f:
                    _extract_crop(f, granule, out)
                h5path.unlink()
            else:
                with _remote_file(granule) as f:
                    _extract_crop(f, granule, out)


# --- stage 2: interferogram ---------------------------------------------------

def _multilook(a: np.ndarray, looks_y: int, looks_x: int) -> np.ndarray:
    ny, nx = a.shape
    by, bx = ny // looks_y, nx // looks_x
    return (a[:by * looks_y, :bx * looks_x]
            .reshape(by, looks_y, bx, looks_x).mean(axis=(1, 3)))


def _goldstein(ifg: np.ndarray, alpha: float = GOLDSTEIN_ALPHA,
               tile: int = GOLDSTEIN_TILE) -> np.ndarray:
    """Goldstein-Werner spectral filter on the complex interferogram."""
    from numpy.fft import fft2, ifft2
    from scipy.ndimage import uniform_filter

    out = np.zeros_like(ifg)
    weight = np.zeros(ifg.shape, dtype=np.float32)
    step = tile // 2
    window = np.hanning(tile)[:, None] * np.hanning(tile)[None, :]
    ny, nx = ifg.shape
    for r0 in range(0, max(1, ny - step), step):
        for c0 in range(0, max(1, nx - step), step):
            r1, c1 = min(r0 + tile, ny), min(c0 + tile, nx)
            block = np.zeros((tile, tile), dtype=np.complex64)
            block[:r1 - r0, :c1 - c0] = np.nan_to_num(ifg[r0:r1, c0:c1])
            spectrum = fft2(block)
            H = uniform_filter(np.abs(spectrum), size=3)
            H = (H / max(H.max(), 1e-30)) ** alpha
            filtered = ifft2(spectrum * H)
            out[r0:r1, c0:c1] += (filtered[:r1 - r0, :c1 - c0]
                                  * window[:r1 - r0, :c1 - c0])
            weight[r0:r1, c0:c1] += window[:r1 - r0, :c1 - c0]
    good = weight > 1e-6
    out[good] /= weight[good]
    return out


def _epoch_invalid(npz, shape: tuple[int, int]) -> np.ndarray:
    """Invalid-sample mask for one epoch; prints the encoding it sees.

    The provisional GSLC mask raster uses 255 as fill; 0 marks samples in no
    subswath. `inputDataExceptionMask` flags degraded input samples. Any
    invalid epoch sample zeroes the interferogram pixel, which drives block
    coherence to zero rather than silently keeping subswath-edge phase — the
    failure mode the valid-sample mask guards against in the GUNW products.
    """
    invalid = np.zeros(shape, dtype=bool)
    if "mask" in npz.files:
        m = npz["mask"]
        values, counts = np.unique(m, return_counts=True)
        print(f"    mask values: {dict(zip(values.tolist(), counts.tolist()))}")
        invalid |= (m == 255) | (m == 0)
    if "inputDataExceptionMask" in npz.files:
        invalid |= npz["inputDataExceptionMask"] != 0
    return invalid


def form_interferogram(ref, sec, tag: str, looks: int = LOOKS) -> dict:
    """Multilooked, Goldstein-filtered interferogram of two crops on one grid:
    filtered phase, boxcar coherence and the multilooked x, y."""
    from scipy.ndimage import uniform_filter

    looks_y, looks_x = 2 * looks, looks
    if not (np.array_equal(ref["x"], sec["x"]) and np.array_equal(ref["y"], sec["y"])):
        raise ValueError(f"{tag}: reference and secondary grids differ")
    z1, z2 = ref["z"], sec["z"]
    print(f"{tag}: forming {z1.shape} interferogram, "
          f"{looks_y} x {looks_x} looks")
    for z, npz in ((z1, ref), (z2, sec)):
        z[_epoch_invalid(npz, z.shape)] = 0
    ny, nx = z1.shape
    by, bx = ny // looks_y, nx // looks_x
    ifg = np.empty((by, bx), dtype=np.complex64)
    p1 = np.empty((by, bx), dtype=np.float32)
    p2 = np.empty((by, bx), dtype=np.float32)
    # Band-wise accumulation keeps peak memory at the two input rasters
    # plus one band of products (the machine has 15 GB total).
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
    x = ref["x"][:bx * looks_x].reshape(-1, looks_x).mean(axis=1)
    y = ref["y"][:by * looks_y].reshape(-1, looks_y).mean(axis=1)
    return {"phase": np.angle(filtered).astype(np.float32),
            "coherence": coherence.astype(np.float32), "x": x, "y": y,
            "epsg": ref["epsg"], "looks": np.array([looks_y, looks_x])}


def make_ifg(tags: list[str] | None = None, looks: int = LOOKS) -> None:
    for tag in PAIRS:
        if tags and tag not in tags:
            continue
        ref = np.load(DATA_GSLC / f"{tag}_ref.npz")
        sec = np.load(DATA_GSLC / f"{tag}_sec.npz")
        product = form_interferogram(ref, sec, tag, looks)
        out = DATA_GSLC / f"{tag}_ifg.npz"
        np.savez_compressed(out, granules=np.array(list(PAIRS[tag])), **product)
        print(f"  saved {out.name} {product['phase'].shape}")


# --- stage 3: unwrap ----------------------------------------------------------

def _run_snaphu(phase: np.ndarray, coherence: np.ndarray,
                *, stat_cost_mode: str = "DEFO",
                init_method: str = "MCF",
                execution_identity: dict | None = None) -> np.ndarray:
    if stat_cost_mode not in {"DEFO", "SMOOTH"}:
        raise ValueError(f"unsupported snaphu cost mode: {stat_cost_mode}")
    if init_method not in {"MCF", "MST"}:
        raise ValueError(f"unsupported snaphu initialization: {init_method}")
    located = shutil.which(SNAPHU)
    if located is None:
        raise FileNotFoundError(f"snaphu executable not found: {SNAPHU}")
    executable = Path(located).resolve(strict=True)
    if execution_identity is not None:
        # snaphu 2.0.6 prints its version with a nonzero exit status. Python
        # package metadata does not identify this command-line executable.
        version = subprocess.run([str(executable), "--version"],
                                 capture_output=True, text=True, check=False,
                                 timeout=10)
        match = re.search(r"\bsnaphu\s+v[^\s]+", version.stdout + "\n" + version.stderr)
        execution_identity.clear()
        execution_identity.update({"resolved_path": str(executable),
                                   "sha256": _sha256(executable),
                                   "reported_version": match.group(0) if match else None,
                                   "version_exit_code": version.returncode})
    ny, nx = phase.shape
    with tempfile.TemporaryDirectory(prefix="snaphu_", dir=str(DATA_GSLC)) as td:
        td = Path(td)
        (td / "phase.in").write_bytes(np.nan_to_num(phase).astype(np.float32).tobytes())
        (td / "corr.in").write_bytes(np.clip(np.nan_to_num(coherence), 0, 1)
                                     .astype(np.float32).tobytes())
        config = td / "snaphu.conf"
        # DEFO cost + MCF init handle the near-axis deformation discontinuity;
        # a single tile avoids the quadrant ambiguity offsets that 2x2 tiling
        # produced (visible against the GUNW field as tile-seam steps).
        config.write_text("\n".join([
            f"STATCOSTMODE {stat_cost_mode}",
            f"INITMETHOD {init_method}",
            "INFILEFORMAT FLOAT_DATA",
            "CORRFILEFORMAT FLOAT_DATA",
            "OUTFILEFORMAT FLOAT_DATA",
            f"CORRFILE {td / 'corr.in'}",
            f"OUTFILE {td / 'unw.out'}",
            f"LINELENGTH {nx}",
        ]) + "\n")
        # LINELENGTH comes from the config file; a positional width argument
        # would be parsed as a second input file.
        result = subprocess.run(
            [str(executable), "-f", str(config), str(td / "phase.in")],
            capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"snaphu failed:\n{result.stdout}\n{result.stderr}")
        unw = np.frombuffer((td / "unw.out").read_bytes(), dtype=np.float32)
    return unw.reshape(ny, nx).copy()


def _iono_on_grid(tag: str, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """GUNW ionosphere screen of the same pair, interpolated to (y, x)."""
    from scipy.interpolate import RegularGridInterpolator

    d = load_pair(tag)
    if d is None:
        raise SystemExit(f"{tag}.npz missing from data/nisar (needed for the "
                         "ionosphere screen)")
    iono = np.asarray(d["ionospherePhaseScreen"], dtype=np.float64)
    gx, gy = d["x"], d["y"]
    iono = np.nan_to_num(iono, nan=0.0)
    interpolator = RegularGridInterpolator(
        (gy[np.argsort(gy)], gx), iono[np.argsort(gy), :],
        bounds_error=False, fill_value=0.0)
    X, Y = np.meshgrid(x, y)
    return interpolator(np.column_stack([Y.ravel(), X.ravel()])).reshape(X.shape)


def _referenced_los(tag: str, unwrapped_phase: np.ndarray,
                    coherence: np.ndarray, x: np.ndarray,
                    y: np.ndarray) -> np.ndarray:
    phase = unwrapped_phase - _iono_on_grid(tag, x, y)
    los = -phase * WAVELENGTH_L / (4 * np.pi)
    X, Y = np.meshgrid(x, y)
    far = (np.hypot(X - CENTER_E, Y - CENTER_N) > FARFIELD_RADIUS) \
        & (coherence > UNWRAP_COHERENCE_MIN)
    if far.sum() < 3:
        raise ValueError(f"{tag}: too few far-field pixels for plane referencing")
    A = np.c_[np.ones(far.sum()), (X[far] - CENTER_E) / 1e4,
              (Y[far] - CENTER_N) / 1e4]
    coef, *_ = np.linalg.lstsq(A, los[far], rcond=None)
    return los - (coef[0] + coef[1] * (X - CENTER_E) / 1e4
                  + coef[2] * (Y - CENTER_N) / 1e4)


def unwrap(tags: list[str] | None = None) -> None:
    for tag in PAIRS:
        if tags and tag not in tags:
            continue
        f = np.load(DATA_GSLC / f"{tag}_ifg.npz")
        phase, coherence, x, y = f["phase"], f["coherence"], f["x"], f["y"]
        print(f"{tag}: snaphu on {phase.shape}")
        unw = _run_snaphu(phase, coherence, **SNAPHU_PRIMARY)
        los = _referenced_los(tag, unw, coherence, x, y)
        X, Y = np.meshgrid(x, y)
        out = DATA_GSLC / f"{tag}_unw.npz"
        np.savez_compressed(out, los=los.astype(np.float32),
                            coherence=coherence, x=x, y=y, epsg=f["epsg"],
                            note=np.array("LOS m, +toward satellite; GUNW iono "
                                          "subtracted; far-field plane removed"))
        print(f"  saved {out.name}; LOS p1/p99 within 25 km: ", end="")
        near = np.hypot(X - CENTER_E, Y - CENTER_N) < 25_000
        ok = near & (coherence > UNWRAP_COHERENCE_MIN)
        print(f"{100 * np.percentile(los[ok], 1):+.1f} / "
              f"{100 * np.percentile(los[ok], 99):+.1f} cm")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unwrap_sensitivity(tags: list[str] | None = None) -> None:
    """Compare SMOOTH costs with the saved DEFO-cost primary unwrap.

    This is an internal processing sensitivity. It does not provide an
    independent unwrap because all cases use the same interferogram,
    coherence, ionosphere screen, and far-field reference.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .params import INVERSION_PAIRS

    selected = list(INVERSION_PAIRS) if not tags else tags
    unknown = sorted(set(selected) - PAIRS.keys())
    if unknown:
        raise ValueError("unknown GSLC tag(s): " + ", ".join(unknown))
    model = np.asarray(json.loads((RESULTS / "okada_fit.json").read_text())["values"][:7])
    cycle_m = WAVELENGTH_L / 2

    for tag in selected:
        ifg_path = DATA_GSLC / f"{tag}_ifg.npz"
        primary_path = DATA_GSLC / f"{tag}_unw.npz"
        f = np.load(ifg_path)
        primary = np.load(primary_path)
        phase = f["phase"]
        coherence = f["coherence"]
        x, y = f["x"], f["y"]
        if primary["los"].shape != phase.shape:
            raise ValueError(f"{tag}: saved unwrap and interferogram shapes differ")
        baseline = np.asarray(primary["los"], dtype=float)
        X, Y = np.meshgrid(x, y)
        perp, along = _axis_distance(X, Y, model)
        radius = np.hypot(X - CENTER_E, Y - CENTER_N)
        common = (coherence > UNWRAP_COHERENCE_MIN) & np.isfinite(baseline)
        regions = {
            "whole_scene": np.ones_like(common),
            "near_axis": ((perp < NEAR_AXIS_PERP_M)
                          & (along < NEAR_AXIS_ALONG_M)),
            "off_axis_within_25km": ((perp >= NEAR_AXIS_PERP_M)
                                     | (along >= NEAR_AXIS_ALONG_M))
                                    & (radius < 25_000),
        }

        fig, axes = plt.subplots(1, len(SNAPHU_SENSITIVITY_VARIANTS),
                                 figsize=(7 * len(SNAPHU_SENSITIVITY_VARIANTS), 6),
                                 squeeze=False)
        results = {}
        for column, (name, config) in enumerate(
                SNAPHU_SENSITIVITY_VARIANTS.items()):
            print(f"{tag}: snaphu sensitivity {name} on {phase.shape}")
            alternative_phase = _run_snaphu(phase, coherence, **config)
            alternative = _referenced_los(
                tag, alternative_phase, coherence, x, y
            )
            valid = common & np.isfinite(alternative)
            difference = alternative - baseline

            def region_stats(region: np.ndarray) -> dict:
                chosen = valid & region
                values = difference[chosen]
                if not len(values):
                    return {"n_pixels": 0}
                return {
                    "n_pixels": int(len(values)),
                    "mean_difference_m": float(np.mean(values)),
                    "rms_difference_m": float(np.sqrt(np.mean(values ** 2))),
                    "p95_abs_difference_m": float(
                        np.percentile(np.abs(values), 95)
                    ),
                    "fraction_over_half_fringe": float(
                        np.mean(np.abs(values) > cycle_m / 2)
                    ),
                }

            results[name] = {
                "configuration": config,
                "regions": {region_name: region_stats(region)
                            for region_name, region in regions.items()},
            }
            ax = axes[0, column]
            image = ax.imshow(
                np.where(valid, difference, np.nan) * 100,
                cmap="RdBu_r", vmin=-100 * cycle_m, vmax=100 * cycle_m,
                extent=[x[0] / 1e3, x[-1] / 1e3,
                        y[-1] / 1e3, y[0] / 1e3],
            )
            plt.colorbar(image, ax=ax, shrink=0.78, label="alternative − primary (cm)")
            ax.set_title(f"{name}: {config['stat_cost_mode']} + "
                         f"{config['init_method']}", fontsize=10)
            ax.set_xlabel("UTM 37N easting (km)")
            if column == 0:
                ax.set_ylabel("UTM 37N northing (km)")

        payload = {
            "tag": tag,
            "ifg": {"path": str(ifg_path.relative_to(REPO)),
                    "sha256": _sha256(ifg_path)},
            "primary_unwrap": {
                "path": str(primary_path.relative_to(REPO)),
                "sha256": _sha256(primary_path),
                "configuration": SNAPHU_PRIMARY,
            },
            "one_fringe_los_m": cycle_m,
            "variants": results,
            "scope": "Same interferogram and corrections in every case. This "
                     "tests processing sensitivity and is not an independent "
                     "unwrap validation.",
        }
        json_path = RESULTS / f"gslc_unwrap_sensitivity_{tag}.json"
        json_path.write_text(json.dumps(payload, indent=2) + "\n")
        fig.suptitle(f"{tag}: snaphu sensitivity relative to the saved "
                     "DEFO + MCF single-tile unwrap", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        figure_path = RESULTS / f"gslc_unwrap_sensitivity_{tag}.png"
        fig.savefig(figure_path, dpi=110)
        plt.close(fig)
        print(f"sensitivity: {json_path}")
        print(f"figure: {figure_path}")


# Acceptance limits for the GUNW comparison. The near-axis band genuinely
# differs between the processings (that difference is the reprocessing's
# purpose), so its limit is loose; everywhere else the two fields must agree.
COMPARE_LIMITS = {
    "slope": (0.90, 1.10),
    "r_min": 0.97,
    "rms_all_max_m": 0.03,
    "rms_near_axis_max_m": 0.15,
}
NEAR_AXIS_PERP_M = 3_000.0
NEAR_AXIS_ALONG_M = 20_000.0


def _axis_distance(X: np.ndarray, Y: np.ndarray, model: np.ndarray):
    """(|perpendicular|, |along-strike|) distance to the fitted dike trace."""
    e_c, n_c, strike = model[0], model[1], model[2]
    s = np.radians(strike)
    d_n, d_e = Y - n_c, X - e_c
    along = d_n * np.cos(s) + d_e * np.sin(s)
    perp = d_n * np.cos(s - np.pi / 2) + d_e * np.sin(s - np.pi / 2)
    return np.abs(perp), np.abs(along)


def compare(tags: list[str] | None = None) -> None:
    """Validate each unwrapped GSLC field against the GUNW product field.

    Both fields are binned to the shared 1.28-km lattice. The regression
    slope catches a phase-sign error immediately (expected slope ≈ +1). The
    statistics are reported per region — whole scene, the near-axis band
    (where the processings are expected to differ), and off-axis within
    25 km — written to JSON, and checked against acceptance limits. GSLC-only
    cells (kept by this processing, absent from GUNW) are counted; they have
    no independent check here.
    """
    import json as _json

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .decompose3d import bin_cells, lattice_axes
    from .gunw import corrected_los, reference_farfield
    from .params import INVERSION_PAIRS, RESULTS

    fit_path = RESULTS / "okada_fit.json"
    model = np.asarray(_json.loads(fit_path.read_text())["values"][:7])

    for tag in PAIRS:
        if tags and tag not in tags:
            continue
        f = np.load(DATA_GSLC / f"{tag}_unw.npz")
        los_g = np.where(f["coherence"] > UNWRAP_COHERENCE_MIN, f["los"], np.nan)
        Xg, Yg = np.meshgrid(f["x"], f["y"])
        cells_g = bin_cells(Xg, Yg, {"los": np.nan_to_num(los_g)},
                            np.isfinite(los_g), 200)["los"]

        d = load_pair(tag)
        los_u, _ = corrected_los(d)
        los_u = reference_farfield(los_u, d["x"], d["y"])
        Xu, Yu = np.meshgrid(d["x"], d["y"])
        cells_u = bin_cells(Xu, Yu, {"los": np.nan_to_num(los_u)},
                            np.isfinite(los_u), 60)["los"]

        centres_e, centres_n = lattice_axes()
        CEg, CNg = np.meshgrid(centres_e, centres_n)
        perp, along = _axis_distance(CEg, CNg, model)
        radius = np.hypot(CEg - CENTER_E, CNg - CENTER_N)
        near_axis = (perp < NEAR_AXIS_PERP_M) & (along < NEAR_AXIS_ALONG_M)

        both = np.isfinite(cells_g) & np.isfinite(cells_u)
        gslc_only = np.isfinite(cells_g) & ~np.isfinite(cells_u)
        slope = float(np.sum(cells_g[both] * cells_u[both])
                      / np.sum(cells_u[both] ** 2))
        r = float(np.corrcoef(cells_g[both], cells_u[both])[0, 1])
        diff = cells_g - cells_u

        def region_stats(mask: np.ndarray) -> dict:
            sel = both & mask
            if sel.sum() < 5:
                return {"n_cells": int(sel.sum())}
            vals = diff[sel]
            return {"n_cells": int(sel.sum()),
                    "rms_m": float(np.sqrt(np.mean(vals ** 2))),
                    "p95_abs_m": float(np.percentile(np.abs(vals), 95))}

        stats = {
            "tag": tag,
            "slope_gslc_on_gunw": slope,
            "correlation": r,
            "whole_scene": region_stats(np.ones_like(both)),
            "near_axis": region_stats(near_axis),
            "off_axis_within_25km": region_stats(~near_axis & (radius < 25_000)),
            "gslc_only_cells": int(gslc_only.sum()),
            "gslc_only_cells_near_axis": int((gslc_only & near_axis).sum()),
            "note": "GSLC-only cells have no GUNW check; the near-axis band "
                    "is where the two processings genuinely differ",
            "limits": COMPARE_LIMITS,
        }
        (RESULTS / f"gslc_gunw_compare_{tag}.json").write_text(
            _json.dumps(stats, indent=2) + "\n")
        rms = stats["whole_scene"]["rms_m"]
        print(f"{tag}: {int(both.sum())} common cells, slope {slope:+.3f}, "
              f"r {r:.3f}; RMS diff whole {100 * rms:.1f} cm, near-axis "
              f"{100 * stats['near_axis'].get('rms_m', float('nan')):.1f} cm "
              f"(p95 {100 * stats['near_axis'].get('p95_abs_m', float('nan')):.1f} cm), "
            f"off-axis {100 * stats['off_axis_within_25km'].get('rms_m', float('nan')):.2f} cm; "
              f"{stats['gslc_only_cells']} GSLC-only cells")

        # The limits are hard requirements only for the pairs that feed the
        # inversion. For other pairs a disagreement can itself be the finding
        # (a product unwrapping error), so it is reported, not fatal.
        violations = []
        lo, hi = COMPARE_LIMITS["slope"]
        if not (lo <= slope <= hi):
            violations.append(f"slope {slope:.3f} outside [{lo}, {hi}]")
        if r < COMPARE_LIMITS["r_min"]:
            violations.append(f"correlation {r:.3f} below {COMPARE_LIMITS['r_min']}")
        if rms > COMPARE_LIMITS["rms_all_max_m"]:
            violations.append(f"whole-scene RMS {rms:.3f} m over limit")
        near_rms = stats["near_axis"].get("rms_m")
        if near_rms is not None and near_rms > COMPARE_LIMITS["rms_near_axis_max_m"]:
            violations.append(f"near-axis RMS {near_rms:.3f} m over limit")
        if violations:
            if tag in INVERSION_PAIRS:
                raise ValueError(f"{tag}: " + "; ".join(violations))
            print(f"  {tag}: GUNW disagreement — " + "; ".join(violations))

        centres_e, centres_n = lattice_axes()
        extent = [centres_e[0] / 1e3, centres_e[-1] / 1e3,
                  centres_n[0] / 1e3, centres_n[-1] / 1e3]
        fig, axes = plt.subplots(1, 3, figsize=(16, 6))
        v = np.nanpercentile(np.abs(cells_u), 99)
        for ax, F, title in ((axes[0], cells_g, "GSLC 30 m reprocessing"),
                             (axes[1], cells_u, "GUNW product (80 m)"),
                             (axes[2], cells_g - cells_u, "difference")):
            im = ax.imshow(F * 100, cmap="RdBu_r", origin="lower",
                           vmin=-100 * v, vmax=100 * v, extent=extent)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel("UTM 37N easting (km)", fontsize=8)
            plt.colorbar(im, ax=ax, shrink=0.75, label="LOS (cm)")
        fig.suptitle(f"{tag}: unwrapped LOS, GSLC reprocessing vs GUNW "
                     f"(1.28-km cells: slope {slope:+.2f}, r {r:.3f}, "
                     f"RMS diff {100 * rms:.1f} cm)", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        out = RESULTS / f"gslc_gunw_compare_{tag}.png"
        fig.savefig(out, dpi=110)
        plt.close(fig)
        print(f"figure: {out}")


# Fixed geographic frame for the full-resolution wrapped-phase figures, so
# runs of different pairs and different processors can be laid side by side.
WRAP_FIG_EXTENT = (40.25, 12.05, 40.90, 12.53)  # west, south, east, north


def wrapped_figure(tags: list[str] | None = None) -> None:
    """Unfiltered 10-m wrapped phase per pair on a fixed lat/lon frame."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pyproj import Transformer

    from .params import RESULTS

    for tag in PAIRS:
        if tags and tag not in tags:
            continue
        ref = np.load(DATA_GSLC / f"{tag}_ref.npz")
        sec = np.load(DATA_GSLC / f"{tag}_sec.npz")
        z1, z2 = ref["z"], sec["z"]
        for z, npz in ((z1, ref), (z2, sec)):
            z[_epoch_invalid(npz, z.shape)] = 0
        ny, nx = z1.shape
        by = ny // 2
        phase = np.empty((by, nx), dtype=np.float32)
        power = np.empty((by, nx), dtype=np.float32)
        for b0 in range(0, by, 3000):
            b1 = min(b0 + 3000, by)
            ifg = (z1[b0 * 2:b1 * 2] * np.conj(z2[b0 * 2:b1 * 2])) \
                .reshape(b1 - b0, 2, nx).mean(axis=1)
            phase[b0:b1] = np.angle(ifg)
            power[b0:b1] = np.abs(ifg)
        x = ref["x"]
        y = ref["y"][:by * 2].reshape(-1, 2).mean(axis=1)

        west, south, east, north = WRAP_FIG_EXTENT
        step = 1.3e-4
        lons = np.arange(west, east, step)
        lats = np.arange(north, south, -step)
        to_utm = Transformer.from_crs(4326, int(ref["epsg"]), always_xy=True)
        out = np.full((len(lats), len(lons)), np.nan, dtype=np.float32)
        for i0 in range(0, len(lats), 400):
            i1 = min(i0 + 400, len(lats))
            LON, LAT = np.meshgrid(lons, lats[i0:i1])
            X, Y = to_utm.transform(LON, LAT)
            ix = np.round((X - x[0]) / 10.0).astype(int)
            iy = np.round((y[0] - Y) / 10.0).astype(int)
            ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < by)
            vals = np.full(ix.shape, np.nan, dtype=np.float32)
            vals[ok] = phase[iy[ok], ix[ok]]
            p = np.zeros(ix.shape, dtype=np.float32)
            p[ok] = power[iy[ok], ix[ok]]
            vals[p == 0] = np.nan
            out[i0:i1] = vals

        def acquisition_date(granule: str) -> str:
            for token in granule.split("_"):
                if len(token) > 8 and token[:8].isdigit() and "T" in token:
                    return f"{token[0:4]}-{token[4:6]}-{token[6:8]}"
            return "?"

        fig, ax = plt.subplots(figsize=(9.5, 7.5))
        im = ax.imshow(out, cmap="twilight_shifted", vmin=-np.pi, vmax=np.pi,
                       extent=[west, east, south, north],
                       aspect=1 / np.cos(np.radians(0.5 * (south + north))))
        granules = PAIRS[tag]
        ax.set_title(f"{tag}: {acquisition_date(granules[0])} → "
                     f"{acquisition_date(granules[1])}, GSLC interferogram, "
                     "10 m, unfiltered wrapped phase", fontsize=10)
        ax.set_xlabel("longitude (°E)")
        ax.set_ylabel("latitude (°N)")
        plt.colorbar(im, ax=ax, shrink=0.85, label="phase (rad)")
        fig.tight_layout()
        out_path = RESULTS / f"gslc_{tag}_wrapped10m.png"
        fig.savefig(out_path, dpi=110)
        plt.close(fig)
        print(f"figure: {out_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--crop", action="store_true")
    ap.add_argument("--download", action="store_true",
                    help="whole-granule download instead of ranged reads")
    ap.add_argument("--ifg", action="store_true")
    ap.add_argument("--unwrap", action="store_true")
    ap.add_argument("--unwrap-sensitivity", action="store_true",
                    help="compare SMOOTH/MCF with the saved DEFO/MCF "
                         "primary unwrap")
    ap.add_argument("--compare", action="store_true",
                    help="validate unwrapped GSLC LOS against the GUNW field")
    ap.add_argument("--wrapped", action="store_true",
                    help="full-resolution wrapped-phase figure per pair")
    ap.add_argument("--looks", type=int, default=LOOKS,
                    help="x looks (y looks = 2x; the grid is 10 m x, 5 m y)")
    ap.add_argument("--tags", nargs="*")
    args = ap.parse_args()
    if args.probe:
        probe(PAIRS["172_008_A"][0])
    if args.crop:
        crop(download=args.download, tags=args.tags)
    if args.ifg:
        make_ifg(tags=args.tags, looks=args.looks)
    if args.unwrap:
        unwrap(tags=args.tags)
    if args.unwrap_sensitivity:
        unwrap_sensitivity(tags=args.tags)
    if args.compare:
        compare(tags=args.tags)
    if args.wrapped:
        wrapped_figure(tags=args.tags)
    if not any((args.probe, args.crop, args.ifg, args.unwrap,
                args.unwrap_sensitivity, args.compare, args.wrapped)):
        ap.print_help()


if __name__ == "__main__":
    main()
