"""The September 2005 dyke in Envisat: the track-49 pair processed here.

Wright et al. (2006, Fig. 2a) show the descending track-49 interferogram of
6 May to 28 October 2005. The ESA archive holds that pair (frame 3357, IS2,
`ASA_IMS_1PNESA20050506_071516_000000152037_00049_16636_0000` and
`ASA_IMS_1PNESA20051028_071514_000000152042_00049_19141_0000`); it was
processed with GMTSAR 6.2 (`p2p_processing.csh ENVI_SLC`: DORIS precise
orbits, SRTM 3-arcsecond terrain with the EGM2008 geoid added, 200-m
Gaussian filter, Goldstein filter, decimation 2, no unwrapping), the cosine
and sine of the filtered phase geocoded separately with `proj_ra2ll.csh`,
and this module turns them and the geocoded coherence into the display
product of Figure 1 panel (b), `data/envisat-2005/t049_20050506_20051028.npz`.
The phase is taken after geocoding: GMTSAR projects grids through a block
median, which on wrapped radians mixes values across the +-pi boundary (its
own `phasefilt_ll.grd` is wrong by up to half a cycle at such nodes).

The product counts the way Ian Hamling's track-300 grids do, so panels (b)
and (c) share one colour convention: `range_change_cm` is the range change
wrapped at one fringe (lambda/2 = 2.8118 cm), increasing with motion AWAY
from the satellite. The sense of the processor's phase is not assumed: the
geocoded phase is compared with the forward model of the September 2005
opening (I. Hamling's revised 10-km model, mh26.hamling) in the pair's own
look geometry, and the sign whose wrapped phase steps agree with the data
is taken; the run fails if the two signs are not clearly separated. The
look geometry (heading, ground-to-satellite azimuth and incidence at the
scene centre) comes from the reference orbit (the GMTSAR LED file) and is
stored with the product for the panel's arrows.

The same step-agreement test, applied to the track-300 grids, is how their
sense was established (October 2008 and September 2006 grids: agreement
0.75-0.86 for z increasing with range increase, 0.0-0.4 for the opposite).

Run from the repository root:

  python -m mh26.envisat_2005 --gmtsar <work dir> [--pair 20050506_20051028]

writes the product, `results/envisat_2005.json` and `results/envisat_2005.md`.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

import numpy as np
from pyproj import Transformer
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from scipy.io import netcdf_file

from . import hamling
from .params import LAME_LAMBDA, RESULTS, SHEAR_MODULUS

REPO = Path(__file__).resolve().parents[1]
PRODUCT = REPO / "data" / "envisat-2005" / "t049_20050506_20051028.npz"
FRINGE_CM = 2.8118            # Envisat lambda/2 in centimetres (lambda 5.6236 cm)
COHERENCE_MIN = 0.15          # below this the geocoded phase is blanked
# The product covers Figure 1's context frame with a margin, on a 0.0005-degree
# grid (about 55 m) made by averaging the unit phasors of the processor's grid.
PRODUCT_EXTENT = {"lon": (40.20, 41.14), "lat": (11.84, 12.92)}
PRODUCT_STEP = 0.0005
STEP_WINDOW = (0.4, 2.6)      # radians: model steps used in the sense test
# The catalogue identifiers of the two scenes; ESA delivers slices generated
# on demand whose names differ in the duration field (18 s instead of 15 s),
# so the files are found by their date, not by these names.
CATALOGUE_SCENES = {
    "20050506": "ASA_IMS_1PNESA20050506_071516_000000152037_00049_16636_0000",
    "20051028": "ASA_IMS_1PNESA20051028_071514_000000152042_00049_19141_0000",
}
TO_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32637", always_xy=True)
WGS84_A, WGS84_F = 6378137.0, 1 / 298.257223563


def read_gmt_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A COARDS GMT grid as (z rows north first, lon ascending, lat descending)."""
    grid = netcdf_file(str(path), "r", mmap=False)
    v = grid.variables
    names = [n for n in ("lon", "x") if n in v], [n for n in ("lat", "y") if n in v]
    if not names[0] or not names[1]:
        raise ValueError(f"{path}: no lon/lat or x/y coordinate variables")
    lon = np.array(v[names[0][0]][:], dtype=np.float64)
    lat = np.array(v[names[1][0]][:], dtype=np.float64)
    lon = np.where(lon < -180.0, lon + 360.0, lon)     # GMTSAR geocodes to -360..0 here
    z = np.array(v["z"][:], dtype=np.float64)
    grid.close()
    if lat[0] < lat[-1]:
        lat, z = lat[::-1], z[::-1]
    return z, lon, lat


def scene_stem(gmtsar_dir: Path, date: str) -> str:
    """The delivered product of `date` (yyyymmdd) among the raw inputs."""
    stems = sorted(p.stem for p in (gmtsar_dir / "raw").glob(f"ASA_IMS_1PNESA{date}_*.N1"))
    if len(stems) != 1:
        raise FileNotFoundError(f"{gmtsar_dir / 'raw'}: expected one ASA_IMS product of {date}, "
                                f"found {stems}")
    return stems[0]


def read_prm(path: Path) -> dict:
    values = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = (s.strip() for s in line.split("=", 1))
            try:
                values[key] = float(value)
            except ValueError:
                values[key] = value
    return values


def read_led(path: Path) -> np.ndarray:
    """State vectors (seconds of day, x, y, z, vx, vy, vz) from a GMTSAR LED."""
    rows = np.loadtxt(path, skiprows=1)
    return np.column_stack([rows[:, 2], rows[:, 3:9]])


def prm_mid_seconds(prm: dict) -> float:
    """Seconds of day at the middle of the scene; the PRM clocks are yyyyddd.fraction."""
    mid = (prm["SC_clock_start"] + prm["SC_clock_stop"]) / 2
    return 86400.0 * (mid - np.floor(mid))


def ecef_to_enu_matrix(lon_deg: float, lat_deg: float) -> np.ndarray:
    lon, lat = np.radians(lon_deg), np.radians(lat_deg)
    return np.array([[-np.sin(lon), np.cos(lon), 0.0],
                     [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)],
                     [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]])


def geodetic_to_ecef(lon_deg: float, lat_deg: float, h: float = 0.0) -> np.ndarray:
    lon, lat = np.radians(lon_deg), np.radians(lat_deg)
    e2 = WGS84_F * (2 - WGS84_F)
    n = WGS84_A / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.array([(n + h) * np.cos(lat) * np.cos(lon), (n + h) * np.cos(lat) * np.sin(lon),
                     (n * (1 - e2) + h) * np.sin(lat)])


def look_geometry(prm: dict, led: np.ndarray, lon0: float, lat0: float) -> dict:
    """Heading, ground-to-satellite azimuth and incidence at (lon0, lat0),
    from the orbit interpolated at the scene's mid time."""
    t_mid = prm_mid_seconds(prm)
    t = led[:, 0]
    if t_mid < t.min() or t_mid > t.max():
        raise ValueError(f"scene mid time {t_mid:.1f} s outside the LED span "
                         f"{t.min():.1f}-{t.max():.1f} s")
    deg = min(5, len(t) - 1)
    pos = np.array([np.polyval(np.polyfit(t - t_mid, led[:, k], deg), 0.0) for k in (1, 2, 3)])
    vel = np.array([np.polyval(np.polyfit(t - t_mid, led[:, k], deg), 0.0) for k in (4, 5, 6)])
    rot = ecef_to_enu_matrix(lon0, lat0)
    ground = geodetic_to_ecef(lon0, lat0)
    u = rot @ (pos - ground)
    u /= np.linalg.norm(u)
    v = rot @ vel
    return {"heading_deg": float(np.degrees(np.arctan2(v[0], v[1])) % 360),
            "los_azimuth_deg": float(np.degrees(np.arctan2(u[0], u[1])) % 360),
            "incidence_deg": float(np.degrees(np.arccos(u[2]))),
            "unit_east_north_up": [float(a) for a in u]}


def model_los_cm(lon: np.ndarray, lat: np.ndarray, unit_enu: np.ndarray,
                 ok: np.ndarray) -> tuple[np.ndarray, dict]:
    """LOS displacement (cm, toward the satellite positive) of the revised
    September 2005 model at the grid nodes flagged `ok`."""
    from pyrocko.modelling import okada_ext

    delivery = hamling.load_openings()
    opening, meta = hamling.load_updated_sept05(delivery)
    edges = delivery["depth_edges_km"][: opening.shape[1] + 1] * 1000.0
    sources, dislocations = [], []
    for c in range(opening.shape[0]):
        for d in range(opening.shape[1]):
            w = float(opening[c, d])
            if w <= 0.0:
                continue
            top, bottom = edges[d], edges[d + 1]
            half_len = 500.0 * float(delivery["columns_width_km"][c])
            sources.append([delivery["columns_n_m"][c], delivery["columns_e_m"][c],
                            0.5 * (top + bottom), delivery["columns_strike_deg"][c], 90.0,
                            -half_len, half_len, -0.5 * (bottom - top), 0.5 * (bottom - top)])
            dislocations.append([0.0, 0.0, w])
    lon2, lat2 = np.meshgrid(lon, lat)
    east, north = TO_UTM.transform(lon2[ok], lat2[ok])
    rec = np.zeros((ok.sum(), 3))
    rec[:, 0], rec[:, 1] = north, east
    res = okada_ext.okada(np.array(sources), np.array(dislocations), rec, LAME_LAMBDA,
                          SHEAR_MODULUS, nthreads=0, rotate_sdn=0, stack_sources=1)
    los = np.full(ok.shape, np.nan)
    los[ok] = 100.0 * (res[:, 1] * unit_enu[0] + res[:, 0] * unit_enu[1] - res[:, 2] * unit_enu[2])
    return los, {"patches": len(sources), "max_opening_m": float(opening.max()),
                 "model_file": meta["file"], "model_sha256": meta["sha256"]}


def step_agreement(phase_grid: np.ndarray, phase_model: np.ndarray, ok: np.ndarray,
                   lag: int) -> tuple[float, int]:
    """Mean cos(grid step - model step) over down and right neighbours at
    `lag` nodes, where the model step lies inside STEP_WINDOW."""
    values, count = [], 0
    for axis in (0, 1):
        n = phase_grid.shape[axis]
        idx = np.minimum(np.arange(n) + lag, n - 1)
        dg = np.angle(np.exp(1j * (np.take(phase_grid, idx, axis=axis) - phase_grid)))
        dm = np.take(phase_model, idx, axis=axis) - phase_model
        use = ok & np.take(ok, idx, axis=axis) & (np.abs(dm) > STEP_WINDOW[0]) \
            & (np.abs(dm) < STEP_WINDOW[1])
        if axis == 0:
            use[n - lag:, :] = False
        else:
            use[:, n - lag:] = False
        values.append(np.cos(dg[use] - dm[use]))
        count += int(use.sum())
    values = np.concatenate(values)
    return (float(values.mean()) if values.size else float("nan")), count


def product_grid(range_phase: np.ndarray, corr: np.ndarray, ok: np.ndarray, lon: np.ndarray,
                 lat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The wrapped range change (cm) and coherence on the product grid: unit
    phasors and coherence averaged from the processor's grid, blank where
    fewer than half of the source nodes are valid."""
    dlon, dlat = abs(lon[1] - lon[0]), abs(lat[1] - lat[0])
    src_transform = from_origin(lon.min() - dlon / 2, lat.max() + dlat / 2, dlon, dlat)
    (x0, x1), (y0, y1) = PRODUCT_EXTENT["lon"], PRODUCT_EXTENT["lat"]
    cols, rows = int(round((x1 - x0) / PRODUCT_STEP)), int(round((y1 - y0) / PRODUCT_STEP))
    dst_transform = from_origin(x0, y1, PRODUCT_STEP, PRODUCT_STEP)
    lon_out = x0 + (np.arange(cols) + 0.5) * PRODUCT_STEP
    lat_out = y1 - (np.arange(rows) + 0.5) * PRODUCT_STEP

    def average(field: np.ndarray) -> np.ndarray:
        destination = np.full((rows, cols), np.nan, dtype=np.float32)
        reproject(source=np.ascontiguousarray(field, dtype=np.float32), destination=destination,
                  src_transform=src_transform, src_crs="EPSG:4326", src_nodata=np.nan,
                  dst_transform=dst_transform, dst_crs="EPSG:4326", dst_nodata=np.nan,
                  resampling=Resampling.average)
        return destination.astype(np.float64)

    weight = average(ok.astype(np.float64))
    real = average(np.where(ok, np.cos(range_phase), 0.0))
    imag = average(np.where(ok, np.sin(range_phase), 0.0))
    keep = np.isfinite(weight) & (weight > 0.5)
    wrapped = np.mod(np.arctan2(imag, real), 2 * np.pi) / (2 * np.pi) * FRINGE_CM
    coherence = average(np.where(np.isfinite(corr), corr, np.nan))
    return np.where(keep, wrapped, np.nan), coherence, lon_out, lat_out


def build(gmtsar_dir: Path, pair: str, out: Path = PRODUCT) -> dict:
    first, second = pair.split("_")
    # GMTSAR names the pair directory by year and day of year.
    doy = [dt.datetime.strptime(d, "%Y%m%d").strftime("%Y%j") for d in (first, second)]
    intf = gmtsar_dir / "intf" / f"{doy[0]}_{doy[1]}"
    if not intf.is_dir():
        candidates = sorted(p.name for p in (gmtsar_dir / "intf").glob("*"))
        raise FileNotFoundError(f"{intf} missing; interferogram directories: {candidates}")
    for name in ("phasefilt_cos_ll.grd", "phasefilt_sin_ll.grd"):
        if not (intf / name).is_file():
            raise FileNotFoundError(
                f"{intf / name} missing: geocode the cosine and sine of the filtered phase "
                "(gmt grdmath phasefilt.grd COS = phasefilt_cos.grd, SIN likewise, then "
                "proj_ra2ll.csh ../../topo/trans.dat on each); the phase must be taken after "
                "geocoding")
    cos, lon, lat = read_gmt_grid(intf / "phasefilt_cos_ll.grd")
    sin, lon_s, lat_s = read_gmt_grid(intf / "phasefilt_sin_ll.grd")
    corr, lon_c, lat_c = read_gmt_grid(intf / "corr_ll.grd")
    same = (cos.shape == sin.shape == corr.shape and np.allclose(lon, lon_s) and np.allclose(lon, lon_c)
            and np.allclose(lat, lat_s) and np.allclose(lat, lat_c))
    if not same:
        raise ValueError("phasefilt_cos_ll, phasefilt_sin_ll and corr_ll are on different grids")
    phase = np.arctan2(sin, cos)
    ok = np.isfinite(cos) & np.isfinite(sin) & np.isfinite(corr) & (corr >= COHERENCE_MIN)

    scenes = {date: scene_stem(gmtsar_dir, date) for date in (first, second)}
    prm = read_prm(gmtsar_dir / "raw" / f"{scenes[first]}.PRM")
    led = read_led(gmtsar_dir / "raw" / f"{scenes[first]}.LED")
    wavelength_cm = 100.0 * float(prm["radar_wavelength"])
    if abs(wavelength_cm / 2 - FRINGE_CM) > 0.001:
        raise ValueError(f"PRM wavelength {wavelength_cm:.4f} cm is not Envisat's")
    lon0, lat0 = float(np.mean(lon[np.any(ok, axis=0)])), float(np.mean(lat[np.any(ok, axis=1)]))
    view = look_geometry(prm, led, lon0, lat0)

    los, model_meta = model_los_cm(lon, lat, np.array(view["unit_east_north_up"]), ok)
    scores = {}
    for sign in (+1, -1):
        phase_model = sign * 2 * np.pi * los / FRINGE_CM
        scores[sign] = {f"lag{lag}": step_agreement(phase, np.where(ok, phase_model, 0.0), ok, lag)
                        for lag in (1, 2, 4)}
    mean_score = {s: float(np.nanmean([v[0] for v in scores[s].values()])) for s in scores}
    toward = max(mean_score, key=mean_score.get)
    if not (mean_score[toward] > 0.5 and mean_score[-toward] < 0.3):
        raise ValueError(f"the sense of the processed phase is not clear: step agreement "
                         f"{mean_score[+1]:+.3f} for phase increasing toward the satellite, "
                         f"{mean_score[-1]:+.3f} for the opposite")
    # Hamling's grids and panel (c) count range INCREASE; `toward` is the sign
    # of the processor's phase per unit of motion toward the satellite.
    range_phase = -toward * phase
    range_change_cm, coherence_out, lon_out, lat_out = product_grid(range_phase, corr, ok, lon, lat)

    out = out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    # Half precision keeps the file at 8 MB: 0.002 cm steps on the wrapped
    # range change and 0.001 on the coherence, far below what is drawn.
    np.savez_compressed(
        out, range_change_cm=range_change_cm.astype(np.float16),
        coherence=coherence_out.astype(np.float16),
        lon=lon_out, lat=lat_out, fringe_cm=FRINGE_CM, coherence_min=COHERENCE_MIN,
        heading_deg=view["heading_deg"], los_azimuth_deg=view["los_azimuth_deg"],
        incidence_deg=view["incidence_deg"], unit_east_north_up=view["unit_east_north_up"],
        scenes=np.array([scenes[first], scenes[second]]), dates=np.array([first, second]),
        catalogue_scenes=np.array([CATALOGUE_SCENES[first], CATALOGUE_SCENES[second]]),
        track=49, sense="range change (cm) wrapped at one fringe, increasing with motion "
                        "away from the satellite, as in the Hamling track-300 grids")
    report = {
        "product": str(out.relative_to(REPO)) if out.is_relative_to(REPO) else str(out), "pair": pair,
        "scenes": [scenes[first], scenes[second]],
        "catalogue_scenes": [CATALOGUE_SCENES[first], CATALOGUE_SCENES[second]],
        "gmtsar_dir": str(gmtsar_dir), "intf_dir": intf.name, "grid": {"rows": int(phase.shape[0]), "cols": int(phase.shape[1]),
                                                "lon": [float(lon.min()), float(lon.max())],
                                                "lat": [float(lat.min()), float(lat.max())],
                                                "spacing_deg": [float(abs(lon[1] - lon[0])),
                                                                float(abs(lat[1] - lat[0]))]},
        "product_grid": {"rows": int(range_change_cm.shape[0]), "cols": int(range_change_cm.shape[1]),
                         "step_deg": PRODUCT_STEP, "lon": list(PRODUCT_EXTENT["lon"]),
                         "lat": list(PRODUCT_EXTENT["lat"]),
                         "valid_fraction": float(np.isfinite(range_change_cm).mean())},
        "valid_fraction": float(ok.mean()), "coherence_min": COHERENCE_MIN,
        "coherence_median_valid": float(np.median(corr[ok])),
        "wavelength_cm": wavelength_cm, "fringe_cm": FRINGE_CM, "view": view,
        "model": model_meta, "model_los_cm_range": [float(np.nanmin(los)), float(np.nanmax(los))],
        "step_agreement": {("toward" if s > 0 else "away"): {k: {"score": v[0], "steps": v[1]}
                                                           for k, v in scores[s].items()}
                           for s in scores},
        "processor_phase_increases_with_motion": "toward" if toward > 0 else "away",
        "sense": "range change increasing with motion away from the satellite",
    }
    (RESULTS / "envisat_2005.json").write_text(json.dumps(report, indent=2) + "\n")
    (RESULTS / "envisat_2005.md").write_text(report_text(report))
    return report


def report_text(r: dict) -> str:
    v, g, pg = r["view"], r["grid"], r["product_grid"]
    toward, away = r["step_agreement"]["toward"], r["step_agreement"]["away"]
    lines = [
        "# The September 2005 dyke in Envisat track 49 (this processing)",
        "",
        f"Pair {r['pair']} (descending track 49, frame 3357, IS2): `{r['scenes'][0]}` and "
        f"`{r['scenes'][1]}` (ESA's on-demand slices of the catalogue scenes "
        f"`{r['catalogue_scenes'][0]}` and `{r['catalogue_scenes'][1]}`), the pair of "
        "Wright et al. (2006) Fig. 2a. Processed with GMTSAR 6.2 "
        "(`p2p_processing.csh ENVI_SLC`, DORIS precise orbits, SRTM 3-arcsecond terrain on the "
        "ellipsoid, 200-m filter, Goldstein filter, decimation 2, no unwrapping), the cosine and "
        "sine of the filtered phase geocoded separately and the phase taken after geocoding; "
        f"product `{r['product']}`.",
        "",
        f"Processor grid: {g['rows']} x {g['cols']} nodes at {g['spacing_deg'][0]:.5f} x "
        f"{g['spacing_deg'][1]:.5f} degrees, {g['lon'][0]:.3f}-{g['lon'][1]:.3f} E, "
        f"{g['lat'][0]:.3f}-{g['lat'][1]:.3f} N; {100 * r['valid_fraction']:.1f}% of nodes kept "
        f"(coherence >= {r['coherence_min']}, median kept coherence {r['coherence_median_valid']:.2f}). "
        f"Product grid: {pg['rows']} x {pg['cols']} nodes at {pg['step_deg']} degrees over "
        f"{pg['lon'][0]}-{pg['lon'][1]} E, {pg['lat'][0]}-{pg['lat'][1]} N, unit phasors and "
        f"coherence averaged from the processor grid, {100 * pg['valid_fraction']:.1f}% valid.",
        "",
        f"Look geometry at the scene centre from the reference orbit: heading {v['heading_deg']:.1f} deg, "
        f"ground-to-satellite azimuth {v['los_azimuth_deg']:.1f} deg, incidence {v['incidence_deg']:.1f} deg.",
        "",
        "Sense check against the forward model of I. Hamling's revised September 2005 opening "
        f"({r['model']['patches']} patches, maximum {r['model']['max_opening_m']:.1f} m) in this geometry "
        f"(model LOS {r['model_los_cm_range'][0]:+.0f} to {r['model_los_cm_range'][1]:+.0f} cm): mean cosine "
        "of the wrapped phase-step difference between data and model, for model steps of 0.4-2.6 rad, "
        "at lags of 1, 2 and 4 nodes:",
        "",
        "| processor phase increases with motion | lag 1 | lag 2 | lag 4 |",
        "|---|---|---|---|",
        f"| toward the satellite | {toward['lag1']['score']:+.3f} ({toward['lag1']['steps']}) | "
        f"{toward['lag2']['score']:+.3f} ({toward['lag2']['steps']}) | "
        f"{toward['lag4']['score']:+.3f} ({toward['lag4']['steps']}) |",
        f"| away from the satellite | {away['lag1']['score']:+.3f} ({away['lag1']['steps']}) | "
        f"{away['lag2']['score']:+.3f} ({away['lag2']['steps']}) | "
        f"{away['lag4']['score']:+.3f} ({away['lag4']['steps']}) |",
        "",
        f"The processor's phase increases with motion {r['processor_phase_increases_with_motion']} "
        "from the satellite; the product stores the range change wrapped at one fringe "
        f"({r['fringe_cm']} cm), increasing with motion away, the sense of the Hamling track-300 "
        "grids (established by the same test on the October 2008 and September 2006 grids).",
        "",
        "The agreement is a check of the interferogram's identity and sign, not a fit: the model "
        "was derived from other data and the near-axis fringes are unresolved at this spacing.",
        "",
    ]
    return "\n".join(lines)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gmtsar", type=Path, required=True,
                        help="GMTSAR working directory (with raw/ and intf/)")
    parser.add_argument("--pair", default="20050506_20051028")
    parser.add_argument("--out", type=Path, default=PRODUCT)
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    result = build(args.gmtsar.resolve(), args.pair, args.out)
    print(json.dumps({k: result[k] for k in ("product", "valid_fraction", "view",
                                              "processor_phase_increases_with_motion")}, indent=2))
