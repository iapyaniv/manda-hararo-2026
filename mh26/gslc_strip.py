"""Display-only northern strip of the track-165 co-event interferogram (Figure 1).

The analysis crop of the 165_082_D pair (`mh26.gslc`) is 90 km square around
the dyke and ends at 12.59 N; the GSLC frame continues north. This command
crops a strip north of that window from the same two GSLC granules, forms the
interferogram with the analysis settings (6 x 3 looks, Goldstein filter,
boxcar coherence), subtracts the GUNW ionosphere screen of the same pair over
the strip (streamed from the GUNW granule), removes the far-field plane of the
analysis crop (recovered from the saved products, see `analysis_plane`), and
saves the corrected LOS wrapped to one fringe:

    python -m mh26.gslc_strip        # -> data/gslc/165_082_D_north_display.npz

The strip is not unwrapped and no reported quantity reads it; `mh26.figure1`
draws it north of the analysis crop. Earthdata credentials in ~/.netrc.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from .fetch import COLLECTION as GUNW_COLLECTION, GRIDS as GUNW_GRIDS, detect_pol, read_pairs
from .gslc import (DATA_GSLC, PAIRS, UNWRAP_COHERENCE_MIN, _extract_crop, _iono_on_grid,
                   _remote_file, form_interferogram)
from .params import CENTER_E, CENTER_N, WAVELENGTH_L

TAG = "165_082_D"
STRIP_TAG = f"{TAG}_north"
# Product UTM metres (EPSG:32637): the analysis crop's x range; y from 3 km
# inside the analysis crop (which ends at CENTER_N + 45 km) to 62 km north of
# the centre, 2.5 km beyond the Figure 1 frame (12.72 N).
STRIP_WINDOW = (635_600.0, 725_600.0, CENTER_N + 42_000.0, CENTER_N + 62_000.0)
DISPLAY = DATA_GSLC / f"{STRIP_TAG}_display.npz"
FRINGE_M = WAVELENGTH_L / 2


def crop_strip() -> None:
    for role, granule in zip(("ref", "sec"), PAIRS[TAG]):
        out = DATA_GSLC / f"{STRIP_TAG}_{role}.npz"
        if out.exists():
            print(f"{STRIP_TAG} {role}: already extracted")
            continue
        print(f"{STRIP_TAG} {role}: {granule}")
        with _remote_file(granule) as f:
            _extract_crop(f, granule, out, window=STRIP_WINDOW)


def _remote_gunw(granule: str):
    import earthaccess
    import h5py

    earthaccess.login(strategy="netrc")
    hits = earthaccess.search_data(short_name=GUNW_COLLECTION, granule_name=granule)
    if not hits:
        raise SystemExit(f"GUNW granule not in archive: {granule}")
    url = [u for u in hits[0].data_links() if u.endswith(".h5")][0]
    fs = earthaccess.get_fsspec_https_session()
    fobj = fs.open(url, mode="rb", cache_type="blockcache", block_size=32 * 1024 * 1024)
    return h5py.File(fobj, "r", rdcc_nbytes=64 * 1024 * 1024)


def gunw_ionosphere(window: tuple[float, float, float, float], margin_m: float = 2_000.0) -> dict:
    """The GUNW ionosphere screen (radians, product posting) over `window` plus a margin."""
    granule = read_pairs()[TAG]
    x0, x1, y0, y1 = window
    with _remote_gunw(granule) as f:
        pol = detect_pol(f)
        base = f"{GUNW_GRIDS}/unwrappedInterferogram/{pol}"
        x = f[f"{base}/xCoordinates"][:]
        y = f[f"{base}/yCoordinates"][:]
        ix = np.where((x >= x0 - margin_m) & (x <= x1 + margin_m))[0]
        iy = np.where((y >= y0 - margin_m) & (y <= y1 + margin_m))[0]
        if len(ix) == 0 or len(iy) == 0:
            raise SystemExit("GUNW grid does not cover the strip")
        screen = f[f"{base}/ionospherePhaseScreen"][iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]
    return {"iono": np.asarray(screen, dtype=np.float32), "x": x[ix[0]:ix[-1] + 1],
            "y": y[iy[0]:iy[-1] + 1], "granule": granule}


def interpolate_screen(screen: dict, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    from scipy.interpolate import RegularGridInterpolator

    gy = screen["y"]
    order = np.argsort(gy)
    interpolator = RegularGridInterpolator(
        (gy[order], screen["x"]), np.nan_to_num(screen["iono"], nan=0.0)[order, :],
        bounds_error=False, fill_value=0.0)
    X, Y = np.meshgrid(x, y)
    return interpolator(np.column_stack([Y.ravel(), X.ravel()])).reshape(X.shape)


def analysis_plane() -> dict:
    """The far-field plane that `mh26.gslc` removed from the analysis crop,
    recovered from the saved products.

    The saved LOS is -(unwrapped - iono) * lambda / 4 pi - plane, and the
    unwrapped phase equals the saved filtered phase up to whole cycles, so
    -(phase - iono) - 4 pi LOS / lambda is the plane in radians, wrapped. Its
    slopes come from the mean phase step between neighbouring pixels and the
    offset from the mean phasor after the slopes are removed (known to one
    cycle, which a wrapped display cannot see). The concentration of the
    residual phasors (1 = exact) is returned as the check.
    """
    ifg = np.load(DATA_GSLC / f"{TAG}_ifg.npz")
    unw = np.load(DATA_GSLC / f"{TAG}_unw.npz")
    phase, coherence, x, y = ifg["phase"], ifg["coherence"], ifg["x"], ifg["y"]
    if unw["los"].shape != phase.shape or not np.array_equal(unw["x"], x):
        raise ValueError(f"{TAG}: saved interferogram and unwrapped products differ")
    iono = _iono_on_grid(TAG, x, y)
    residual = np.exp(1j * (-(phase.astype(np.float64) - iono)
                            - 4 * np.pi * unw["los"].astype(np.float64) / WAVELENGTH_L))
    good = (coherence > UNWRAP_COHERENCE_MIN) & np.isfinite(residual)
    residual = np.where(good, residual, 0.0)
    pair_x = residual[:, 1:] * np.conj(residual[:, :-1])
    pair_y = residual[1:, :] * np.conj(residual[:-1, :])
    step_x = float(np.angle(pair_x.sum()))
    step_y = float(np.angle(pair_y.sum()))
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    # radians per pixel -> LOS metres per 1e4 m, the plane's coordinate unit
    c1 = step_x * WAVELENGTH_L * 1e4 / (4 * np.pi * dx)
    c2 = step_y * WAVELENGTH_L * 1e4 / (4 * np.pi * dy)
    X, Y = np.meshgrid(x, y)
    ramp = np.exp(1j * 4 * np.pi * (c1 * (X - CENTER_E) + c2 * (Y - CENTER_N)) / (1e4 * WAVELENGTH_L))
    flat = residual * np.conj(ramp)
    c0 = float(np.angle(flat[good].mean())) * WAVELENGTH_L / (4 * np.pi)
    concentration = float(np.abs(flat[good].mean()))
    return {"c0_m": c0, "c1_m_per_10km": c1, "c2_m_per_10km": c2,
            "residual_concentration": concentration, "pixels": int(good.sum())}


def plane_on(x: np.ndarray, y: np.ndarray, plane: dict) -> np.ndarray:
    X, Y = np.meshgrid(x, y)
    return (plane["c0_m"] + plane["c1_m_per_10km"] * (X - CENTER_E) / 1e4
            + plane["c2_m_per_10km"] * (Y - CENTER_N) / 1e4)


def build() -> None:
    crop_strip()
    ref = np.load(DATA_GSLC / f"{STRIP_TAG}_ref.npz")
    sec = np.load(DATA_GSLC / f"{STRIP_TAG}_sec.npz")
    product = form_interferogram(ref, sec, STRIP_TAG)
    x, y = product["x"], product["y"]
    plane = analysis_plane()
    print(f"analysis plane: {json.dumps(plane)}")
    if plane["residual_concentration"] < 0.99:
        raise ValueError("the recovered plane does not reproduce the saved products")
    screen = gunw_ionosphere(STRIP_WINDOW)
    iono = interpolate_screen(screen, x, y)
    los = -(product["phase"].astype(np.float64) - iono) * WAVELENGTH_L / (4 * np.pi)
    los = np.mod(los - plane_on(x, y, plane), FRINGE_M)
    np.savez_compressed(
        DISPLAY, los=los.astype(np.float32), coherence=product["coherence"], x=x, y=y,
        epsg=product["epsg"], looks=product["looks"],
        granules=np.array(list(PAIRS[TAG]) + [screen["granule"]]),
        window=np.array(STRIP_WINDOW), plane=np.array(json.dumps(plane)),
        note=np.array("LOS m wrapped to one fringe, +toward satellite; GUNW iono "
                      "subtracted; analysis far-field plane removed; display only"))
    print(f"saved {DISPLAY.name} {los.shape}; y {y.min():.0f}-{y.max():.0f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--plane-only", action="store_true",
                        help="print the recovered analysis plane and stop")
    args = parser.parse_args()
    if args.plane_only:
        print(json.dumps(analysis_plane(), indent=1))
        return
    build()


if __name__ == "__main__":
    main()
