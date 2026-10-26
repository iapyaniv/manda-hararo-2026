"""The GRL location figure (Figure 1), in the panel layout agreed with the coauthors.

Panels
  (a) the July 2026 intrusion in the NISAR descending interferogram of 7 to
      19 July (165_082_D, 30-m GSLC chain): the corrected unwrapped LOS field
      the paper models (ionosphere screen subtracted, far-field plane
      removed), screened at the supplement's coherence threshold, smoothed by
      one pixel in the complex domain and rewrapped at one L-band fringe. The
      analysis crop ends at 12.59 N; north of it the panel shows the
      display-only strip of mh26.gslc_strip (same granules and settings, the
      same corrections, not unwrapped), when that product exists; with the
      model dyke plane
      (the part with opening above one third of the peak drawn thick), the
      candidate graben fault, the 2005 dyke (model plane, thick where it
      opened), the one located July 2026 earthquake (two agency locations
      with their stated errors), the volcanic centres at their mapped
      positions and the segment names;
  (b) the September 2005 dyke in the descending track-49 pair of 6 May to
      28 October 2005 (the pair of Wright et al. 2006, Fig. 2a) processed
      here (mh26.envisat_2005), drawn like (c) when that product exists;
      the delivered track-300 set holds no 2005 interferogram (its "sept"
      grid is the September 2006 dyke), so without the product the panel is
      a placeholder with the 2005 model plane;
  (c) the October 2008 dyke in Ian Hamling's wrapped Envisat track-300
      interferogram (Hamling et al. 2010; acquisition dates to be confirmed),
      multilooked as a complex phasor;
  (d) Sentinel-1 descending cumulative LOS displacement 23 October 2014 to
      11 July 2026 (P. Espín Bedón's product; see
      data/sentinel1-cumulative/README) wrapped at one Sentinel-1 fringe, as
      in B. Ireland's figure, so that the local Ado'Ale uplift shows against
      the broad signal; pixels whose neighbourhood phase is inconsistent are
      omitted; with the time-series sites of (e) and (f) and the Sentinel-1
      look direction;
  (e, f) the Ado'Ale series and the graben-minus-flank series of B. Ireland's
      Sentinel-1 time-series figure, replotted from values digitized from it
      (data/ireland-series; DRAFT: to be replaced by the series themselves).

Everything is drawn in geographic coordinates with one Copernicus GLO-30
hillshade (tile N12 E040 from the repository and N11 E040, fetched, for the
southern strip). Nothing here is an analysis quantity: the figure only
displays saved products and the delivered models.

Run from the repository root:

  python -m mh26.figure1 [--dem-south <N11 tile>]          results/figure1.{png,pdf}
  python -m mh26.figure1 --notes                            results/figure1_draft.{png,pdf},
                                                            with the notes to the coauthors
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from cmcrameri import cm as cmc  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LightSource, Normalize  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from pyproj import Transformer  # noqa: E402
from rasterio.transform import from_origin  # noqa: E402
from rasterio.warp import Resampling, reproject  # noqa: E402
from scipy.io import netcdf_file  # noqa: E402
from scipy.ndimage import gaussian_filter, uniform_filter  # noqa: E402

from .hamling import (  # noqa: E402
    EPOCH_TAGS, OPENING_CMAP, PRIMARY_EXCLUDED, load_openings, load_updated_sept05)
from .params import REPO, RESULTS, WAVELENGTH_L  # noqa: E402
from .summary_figures import DEM_PATH, geometry_lines, load_mapped_gslc_fit  # noqa: E402

UTM_CRS = "EPSG:32637"
LL_CRS = "EPSG:4326"
TO_LL = Transformer.from_crs(32637, 4326, always_xy=True)

NISAR_UNW = REPO / "data" / "gslc" / "165_082_D_unw.npz"
# Display-only strip north of the analysis crop (python -m mh26.gslc_strip).
NISAR_NORTH = REPO / "data" / "gslc" / "165_082_D_north_display.npz"
NISAR_CROP_NORTH_LAT = 12.59
NISAR_LABEL = "NISAR, 7 – 19 Jul 2026"
S1_PRE = (REPO / "data" / "sentinel1-cumulative" / "Pre"
          / "pre_079D_07694_131313_20141023_20260711.cum.geo.tif")
NISAR_FRINGE_M = WAVELENGTH_L / 2
NISAR_FRINGE_CM = 100 * NISAR_FRINGE_M
# Sentinel-1 C-band (5.5465 cm): the cumulative map is wrapped at half the
# wavelength, the interval that reproduces B. Ireland's wrapped panel.
S1_FRINGE_MM = 1000 * 0.055465 / 2
# Local phase consistency below which the wrapped Sentinel-1 map is not
# drawn: the mean resultant length, over 5 x 5 pixels, of the phase steps
# between neighbouring pixels. A smooth ramp of any slope keeps a constant
# step and passes; the decorrelated ground west of the rift and the
# unwrapping-noise cells on the axis do not. (The concentration of the
# phasors themselves would blank steep but smooth fringes, such as the
# centre of the Dabbahu bullseye.)
S1_CONSISTENCY_MIN = 0.6
# Ian Hamling's wrapped Envisat track-300 interferograms of the earlier dykes
# (one GMT grid per dyke, wrapped LOS in centimetres over one Envisat fringe).
ENVISAT_DIR = REPO / "data" / "ihamling" / "opening_files"
ENVISAT_FRINGE_CM = 2.8118
# The delivered grids are the 2006 to 2010 dykes ("sept" is September 2006).
ENVISAT_PANELS = {"c": ("oct", "Oct 2008 dike\nEnvisat track 300")}
# Envisat ascending, right-looking: heading about 348 deg, so the
# ground-to-satellite azimuth is about 258 deg (the lobes of the October 2008
# grid match the ascending forward model, not the descending one). The grids
# count range increase: motion away from the satellite raises z.
ENVISAT_VIEW = {"heading_deg": 348.0, "los_deg": 258.0}
# The September 2005 dyke: the track-49 pair of Wright et al. (2006) processed
# here (mh26.envisat_2005), counted like the track-300 grids. Panel (b) is a
# placeholder until this product exists.
ENVISAT_2005 = REPO / "data" / "envisat-2005" / "t049_20050506_20051028.npz"
# Text S1 screens the GSLC observations at this coherence.
COHERENCE_MIN = 0.25
# The thick part of the dyke line: columns whose maximum opening reaches this
# fraction of the peak (17.5 km of the 28 km model plane for the mapped fit).
OPENED_FRACTION = 1 / 3

# Second terrain tile for the strip south of 12 N (same source as DEM_PATH).
DEM_SOUTH_PATH = REPO / "data" / "dem" / "Copernicus_DSM_COG_10_N11_00_E040_00_DEM.tif"
DEM_SOUTH_URL = (
    "https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_N11_00_E040_00_DEM/"
    "Copernicus_DSM_COG_10_N11_00_E040_00_DEM.tif"
)
# The six site series digitized from B. Ireland's figure; see
# data/ireland-series/README.md. Draft stand-in for the series themselves.
SERIES_CSV = REPO / "data" / "ireland-series" / "ireland_079D_sites_digitized.csv"

# Display frames (degrees). Frame A is the NISAR panel; frame B is shared by
# the three context maps so that they compare directly.
FRAME_A = {"lon": (40.24, 41.08), "lat": (11.76, 12.72), "step": 0.0004}
FRAME_B = {"lon": (40.28, 41.06), "lat": (11.92, 12.84), "step": 0.0005}

# Volcanic centres at their mapped positions, not the model sources: Dabbahu
# and Manda Hararo from the Global Volcanism Program; Gabho as located by
# Wright et al. (2006) (quoted by Grandin et al. 2010); the Ado'Ale Volcanic
# Complex at about 12 deg 22 min N (Ebinger et al. 2008), on the rift axis.
VOLCANOES = {
    "Gabho": (40.53, 12.69),
    "Dabbahu": (40.48, 12.595),
    "Ado'Ale": (40.60, 12.37),
    "Manda Hararo": (40.82, 12.17),
}
# The one located July 2026 earthquake at its IGSSA location (compiled by
# D. Keir). Errors in km (latitude, longitude), drawn as plain error bars: no
# joint confidence level is implied. The ISC location of the same event
# (11.79 N, 40.64 E, errors 28 and 16 km) lies on the southern edge of the
# frame and is left out.
EARTHQUAKE = {"lon": 40.853, "lat": 12.042, "err_lat_km": 3.8, "err_lon_km": 3.6}
# B. Ireland's time-series sites (positions read from the figure and
# confirmed against the cumulative map). The graben series is the graben
# site minus the flank site 8.5 km to the north-east.
SITES = {
    "e": {"name": "Ado'Ale", "lon": 40.698, "lat": 12.330, "series": "Ado'Ale"},
    "f": {"name": "graben", "lon": 40.584, "lat": 12.393, "series": "Graben subsidence"},
}
FLANK_REFERENCE = (40.649, 12.436)

# Viewing geometry for the arrows. NISAR is left-looking: on the descending
# pass 165 (heading about 195 deg) the saved ground-to-satellite unit vector
# has azimuth 285 deg and the dyke's west flank is the positive one.
# Sentinel-1 is right-looking: the LiCSAR 079D unit vectors (ground to
# satellite) have azimuth 101 deg, so the heading is about 191 deg.
NISAR_VIEW = {"heading_deg": 195.0, "los_deg": 285.0}
S1_VIEW = {"heading_deg": 191.0, "los_deg": 101.0}

HALO = [path_effects.Stroke(linewidth=2.2, foreground="white"), path_effects.Normal()]
THIN_HALO = [path_effects.Stroke(linewidth=1.6, foreground="white"), path_effects.Normal()]

KM_PER_DEG_LAT = 111.32
FS_TICK, FS_LABEL, FS_TITLE, FS_NOTE = 6.5, 6.2, 6.4, 6.0
FS_MODEL = 6.6   # the three model-geometry labels share one size and weight


def km_per_deg_lon(lat: float) -> float:
    return KM_PER_DEG_LAT * math.cos(math.radians(lat))


# ---------------------------------------------------------------------------
# Raster helpers


def frame_grid(frame: dict) -> tuple[object, tuple[int, int], list[float]]:
    (lon0, lon1), (lat0, lat1), step = frame["lon"], frame["lat"], frame["step"]
    width = int(round((lon1 - lon0) / step))
    height = int(round((lat1 - lat0) / step))
    transform = from_origin(lon0, lat1, step, step)
    return transform, (height, width), [lon0, lon1, lat0, lat1]


def warp(source, src_transform, src_crs, dst_transform, dst_shape, resampling,
         src_nodata=np.nan, destination=None):
    init = destination is None
    if init:
        destination = np.full(dst_shape, np.nan, dtype=np.float32)
    reproject(
        source=np.ascontiguousarray(source, dtype=np.float32),
        destination=destination,
        src_transform=src_transform,
        src_crs=src_crs,
        src_nodata=src_nodata,
        dst_transform=dst_transform,
        dst_crs=LL_CRS,
        dst_nodata=np.nan,
        resampling=resampling,
        init_dest_nodata=init,
    )
    return destination


def terrain_shade(frame: dict, dem_paths: list[Path]) -> np.ndarray:
    transform, shape, _ = frame_grid(frame)
    terrain = None
    for path in dem_paths:
        with rasterio.open(path) as dataset:
            source = dataset.read(1, masked=True).filled(np.nan)
            terrain = warp(source, dataset.transform, dataset.crs, transform, shape,
                           Resampling.average, destination=terrain)
    if not np.isfinite(terrain).any():
        raise ValueError("no terrain in the frame")
    filled = np.where(np.isfinite(terrain), terrain, np.nanmedian(terrain))
    filled = gaussian_filter(filled, sigma=0.7)
    lat_mid = 0.5 * sum(frame["lat"])
    dx = frame["step"] * km_per_deg_lon(lat_mid) * 1000
    dy = frame["step"] * KM_PER_DEG_LAT * 1000
    shade = LightSource(azdeg=315, altdeg=42).hillshade(
        filled, vert_exag=1.3, dx=dx, dy=dy, fraction=1.1)
    return 0.70 + 0.30 * shade


def blend(field: np.ndarray, shade: np.ndarray, cmap, norm, alpha: float) -> np.ndarray:
    background = np.repeat(shade[:, :, None], 3, axis=2)
    colours = cmap(norm(field))[..., :3]
    colours = np.clip(colours * (0.80 + 0.28 * shade[:, :, None]), 0, 1)
    out = background.copy()
    valid = np.isfinite(field)
    out[valid] = (1 - alpha) * background[valid] + alpha * colours[valid]
    return out


def utm_transform(x: np.ndarray, y: np.ndarray) -> object:
    # x ascends, y descends; arrays hold cell centres.
    dx = float(x[1] - x[0])
    dy = float(y[0] - y[1])
    return from_origin(float(x[0]) - dx / 2, float(y[0]) + dy / 2, dx, dy)


def nisar_rewrapped_cm(frame: dict, coherence_min: float) -> tuple[np.ndarray, bool]:
    """The corrected LOS field (metres, positive toward the satellite),
    screened by coherence, wrapped to one fringe and smoothed by one pixel
    and resampled as a complex phasor, so that integer-cycle steps in the
    unwrapped field cannot leak into fractional fringes. The analysis crop is
    drawn wherever it has a value; the display strip fills the frame north of
    it. Returns the field and whether the strip was drawn."""
    wrapped = _rewrapped_cm(NISAR_UNW, frame, coherence_min)
    if not NISAR_NORTH.exists():
        print(f"{NISAR_NORTH.name} missing: (a) ends at the analysis crop "
              f"({NISAR_CROP_NORTH_LAT}°N); run python -m mh26.gslc_strip")
        return wrapped, False
    north = _rewrapped_cm(NISAR_NORTH, frame, coherence_min)
    # Switch to the strip 1 km inside the analysis crop, not at its last row:
    # the smoothing and resampling thin out at the crop edge and would draw a
    # hairline there (the two fields agree over the 3-km overlap).
    crop = np.load(NISAR_UNW, allow_pickle=False)
    _, lat_top = TO_LL.transform(float(np.mean(crop["x"])), float(crop["y"].max()) - 1000.0)
    _, (height, _), _ = frame_grid(frame)
    lat_rows = frame["lat"][1] - (np.arange(height) + 0.5) * frame["step"]
    use_strip = (lat_rows > lat_top)[:, None] & np.isfinite(north)
    return np.where(use_strip, north, wrapped), True


def _rewrapped_cm(path: Path, frame: dict, coherence_min: float) -> np.ndarray:
    data = np.load(path, allow_pickle=False)
    transform, shape, _ = frame_grid(frame)
    src_transform = utm_transform(data["x"], data["y"])
    los = np.asarray(data["los"], dtype=np.float64)
    coherence = np.asarray(data["coherence"], dtype=np.float64)
    valid = np.isfinite(los) & (coherence >= coherence_min)
    phase = 2 * np.pi * np.where(valid, los, 0.0) / NISAR_FRINGE_M
    weight = gaussian_filter(valid.astype(np.float64), 1.0)
    real = gaussian_filter(np.where(valid, np.cos(phase), 0.0), 1.0)
    imag = gaussian_filter(np.where(valid, np.sin(phase), 0.0), 1.0)
    keep = valid & (weight > 0.5)
    real = np.where(keep, real / np.maximum(weight, 1e-6), np.nan)
    imag = np.where(keep, imag / np.maximum(weight, 1e-6), np.nan)
    real_w = warp(real, src_transform, UTM_CRS, transform, shape, Resampling.average)
    imag_w = warp(imag, src_transform, UTM_CRS, transform, shape, Resampling.average)
    wrapped = np.mod(np.arctan2(imag_w, real_w), 2 * np.pi) / (2 * np.pi) * NISAR_FRINGE_CM
    wrapped[~np.isfinite(real_w) | ~np.isfinite(imag_w)] = np.nan
    return wrapped


def phase_step_consistency(unit: np.ndarray, valid: np.ndarray, size: int = 5) -> np.ndarray:
    """Mean resultant length, over size x size pixels, of the phase steps to
    the right and downward neighbours (the smaller of the two directions);
    zero where fewer than half the window's steps exist."""
    out = np.full(unit.shape, np.inf)
    for axis in (0, 1):
        step = np.zeros_like(unit)
        pair = np.zeros(unit.shape, dtype=bool)
        lead = (slice(1, None), slice(None)) if axis == 0 else (slice(None), slice(1, None))
        lag = (slice(None, -1), slice(None)) if axis == 0 else (slice(None), slice(None, -1))
        pair[lag] = valid[lead] & valid[lag]
        step[lag] = np.where(pair[lag], unit[lead] * np.conj(unit[lag]), 0.0)
        count = uniform_filter(pair.astype(np.float64), size)
        length = np.hypot(uniform_filter(step.real, size), uniform_filter(step.imag, size))
        length = np.where(count >= 0.5, length / np.maximum(count, 1e-6), 0.0)
        out = np.minimum(out, length)
    return out


def s1_wrapped_cm(frame: dict) -> np.ndarray:
    """The delivered cumulative map wrapped at one Sentinel-1 fringe, smoothed
    by one pixel and resampled as a complex phasor. The south-western corner
    (west of the rift) holds about -100 mm with 25 to 35 mm pixel scatter,
    which wraps to noise: pixels whose 5 x 5 neighbourhood has inconsistent
    phase steps (phase_step_consistency below S1_CONSISTENCY_MIN) are left
    blank, as are the map's own gaps."""
    transform, shape, _ = frame_grid(frame)
    with rasterio.open(S1_PRE) as dataset:
        source = dataset.read(1, masked=True).filled(np.nan)
        src_transform, src_crs = dataset.transform, dataset.crs
    valid = np.isfinite(source)
    phase = 2 * np.pi * np.where(valid, source, 0.0) / S1_FRINGE_MM
    unit = np.where(valid, np.exp(1j * phase), 0.0)
    unit_real, unit_imag = unit.real.copy(), unit.imag.copy()
    keep = valid & (phase_step_consistency(unit, valid) >= S1_CONSISTENCY_MIN)
    weight = gaussian_filter(keep.astype(np.float64), 1.0)
    real = gaussian_filter(np.where(keep, unit_real, 0.0), 1.0)
    imag = gaussian_filter(np.where(keep, unit_imag, 0.0), 1.0)
    real = np.where(keep, real / np.maximum(weight, 1e-6), np.nan)
    imag = np.where(keep, imag / np.maximum(weight, 1e-6), np.nan)
    real_w = warp(real, src_transform, src_crs, transform, shape, Resampling.average)
    imag_w = warp(imag, src_transform, src_crs, transform, shape, Resampling.average)
    wrapped = np.mod(np.arctan2(imag_w, real_w), 2 * np.pi) / (2 * np.pi) * S1_FRINGE_MM / 10
    wrapped[~np.isfinite(real_w) | ~np.isfinite(imag_w)] = np.nan
    return wrapped


def read_grd(path: Path) -> tuple[np.ndarray, object]:
    """A GMT grid in either classic layout, as (rows north first, transform)."""
    grid = netcdf_file(str(path), "r", mmap=False)
    v = grid.variables
    if "x_range" in v:
        x0, x1 = (float(a) for a in v["x_range"][:])
        y0, y1 = (float(a) for a in v["y_range"][:])
        nx, ny = (int(a) for a in v["dimension"][:])
        z = np.array(v["z"][:], dtype=np.float64).reshape(ny, nx)
        transform = from_origin(x0, y1, (x1 - x0) / nx, (y1 - y0) / ny)
    else:
        x = np.array(v["x"][:], dtype=np.float64)
        y = np.array(v["y"][:], dtype=np.float64)
        z = np.array(v["z"][:], dtype=np.float64)
        if y[0] < y[-1]:
            z = z[::-1]
        dx, dy = abs(x[1] - x[0]), abs(y[1] - y[0])
        transform = from_origin(x.min() - dx / 2, y.max() + dy / 2, dx, dy)
    grid.close()
    return z, transform


def wrapped_cm_on_frame(z_cm: np.ndarray, src_transform, src_crs, frame: dict,
                        looks: int) -> np.ndarray:
    """A wrapped field in centimetres of one Envisat fringe, multilooked as a
    complex phasor over `looks` x `looks` cells and resampled the same way
    onto the frame."""
    ok = np.isfinite(z_cm)
    phase = 2 * np.pi * np.where(ok, z_cm, 0.0) / ENVISAT_FRINGE_CM
    weight = uniform_filter(ok.astype(np.float64), looks)
    real = uniform_filter(np.where(ok, np.cos(phase), 0.0), looks)
    imag = uniform_filter(np.where(ok, np.sin(phase), 0.0), looks)
    keep = ok & (weight > 0.5)
    real = np.where(keep, real / np.maximum(weight, 1e-6), np.nan)
    imag = np.where(keep, imag / np.maximum(weight, 1e-6), np.nan)
    transform, shape, _ = frame_grid(frame)
    real_w = warp(real, src_transform, src_crs, transform, shape, Resampling.average)
    imag_w = warp(imag, src_transform, src_crs, transform, shape, Resampling.average)
    wrapped = np.mod(np.arctan2(imag_w, real_w), 2 * np.pi) / (2 * np.pi) * ENVISAT_FRINGE_CM
    wrapped[~np.isfinite(real_w) | ~np.isfinite(imag_w)] = np.nan
    return wrapped


def envisat_wrapped_cm(frame: dict, name: str, looks: int = 3) -> np.ndarray:
    """One of the track-300 interferograms on the frame, in centimetres of
    one fringe."""
    z, src_transform = read_grd(ENVISAT_DIR / f"t300_{name}_wrp_dyke.grd")
    if not np.isfinite(z).any() or np.nanmax(z) == np.nanmin(z):
        raise ValueError(f"track-300 grid {name!r} holds no interferogram")
    return wrapped_cm_on_frame(z, src_transform, UTM_CRS, frame, looks)


def envisat_2005_wrapped_cm(frame: dict, looks: int = 3) -> tuple[np.ndarray, dict]:
    """The track-49 interferogram of the September 2005 dyke on the frame, in
    centimetres of one fringe, with the look geometry stored in the product."""
    with np.load(ENVISAT_2005, allow_pickle=False) as product:
        z = product["range_change_cm"].astype(np.float64)
        lon, lat = product["lon"], product["lat"]
        view = {"heading_deg": float(product["heading_deg"]),
                "los_deg": float(product["los_azimuth_deg"])}
    dx, dy = abs(lon[1] - lon[0]), abs(lat[1] - lat[0])
    src_transform = from_origin(lon.min() - dx / 2, lat.max() + dy / 2, dx, dy)
    return wrapped_cm_on_frame(z, src_transform, "EPSG:4326", frame, looks), view


def hamling_profiles() -> dict:
    """Column centres (lon, lat) along the Hamling plane with the depth-mean
    (0 to 10 km) opening of the updated September 2005 model and of the
    cumulative 2006 to 2010 sequence (thirteen canonical files)."""
    delivery = load_openings()
    lon, lat = TO_LL.transform(delivery["columns_e_m"], delivery["columns_n_m"])
    cum_tags = [t for t in EPOCH_TAGS[1:] if t not in PRIMARY_EXCLUDED]
    cumulative = sum(delivery["openings_m"][t] for t in cum_tags)
    updated, _ = load_updated_sept05(delivery)
    return {
        "lon": np.asarray(lon),
        "lat": np.asarray(lat),
        "sept05": updated[:, :10].mean(axis=1),
        "cumulative": cumulative[:, :10].mean(axis=1),
        "n_cumulative": len(cum_tags),
    }


def dyke_2026_lines() -> dict:
    """The model plane and candidate fault traces in lon/lat, plus the part of
    the plane whose columns open to at least OPENED_FRACTION of the peak."""
    _, geometry, magnitudes = load_mapped_gslc_fit()
    out = {}
    for name, (e, n) in geometry_lines(geometry).items():
        lon, lat = TO_LL.transform(e, n)
        out[name] = (np.asarray(lon), np.asarray(lat))
    kinds = np.asarray(geometry["kinds"])
    opening = np.asarray(magnitudes)[kinds == "dike"].reshape(tuple(geometry["dike_shape"]))
    # dike_shape is (depth rows, along-strike columns); C order is the saved order.
    column_max = opening.max(axis=0)
    edges = np.asarray(geometry["dike_s_edges_m"], dtype=float)
    opened = np.where(column_max >= OPENED_FRACTION * column_max.max())[0]
    s = np.linspace(edges[opened[0]], edges[opened[-1] + 1], 50)
    e_c, n_c = geometry["trace_centre"]
    strike = math.radians(geometry["strike"])
    lon, lat = TO_LL.transform(e_c + s * math.sin(strike), n_c + s * math.cos(strike))
    out["dike_open"] = (np.asarray(lon), np.asarray(lat))
    out["opened_length_km"] = float((edges[opened[-1] + 1] - edges[opened[0]]) / 1000)
    out["plane_length_km"] = float((edges[-1] - edges[0]) / 1000)
    return out


def site_series(csv_path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Decimal year and LOS (mm) per site name from the digitized series file."""
    if not csv_path.exists():
        raise FileNotFoundError(
            f"{csv_path}: the series behind panels (e) and (f) are missing; see "
            "data/ireland-series/README.md")
    series: dict[str, list[tuple[float, float]]] = {}
    with open(csv_path, newline="") as handle:
        for row in csv.DictReader(handle):
            day = dt.date.fromisoformat(row["epoch"])
            start = dt.date(day.year, 1, 1)
            length = (dt.date(day.year + 1, 1, 1) - start).days
            year = day.year + (day - start).days / length
            series.setdefault(row["site"], []).append((year, float(row["los_mm"])))
    return {site: tuple(np.asarray(v).T) for site, v in series.items()}


# ---------------------------------------------------------------------------
# Drawing helpers


def mm_axes(fig, x_mm, y_mm, w_mm, h_mm, **kwargs):
    width_mm, height_mm = (v * 25.4 for v in fig.get_size_inches())
    ax = fig.add_axes([x_mm / width_mm, y_mm / height_mm, w_mm / width_mm, h_mm / height_mm],
                      **kwargs)
    ax._size_mm = (w_mm, h_mm)
    return ax


def map_axes(ax, frame: dict, *, lon_labels: bool, lat_labels: bool, tick_step=0.2):
    (lon0, lon1), (lat0, lat1) = frame["lon"], frame["lat"]
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(lat0, lat1)
    ax.set_aspect(1 / math.cos(math.radians(0.5 * (lat0 + lat1))))
    lon_ticks = np.arange(math.ceil(lon0 / tick_step) * tick_step, lon1 + 1e-9, tick_step)
    lat_ticks = np.arange(math.ceil(lat0 / tick_step) * tick_step, lat1 + 1e-9, tick_step)
    ax.set_xticks(lon_ticks)
    ax.set_yticks(lat_ticks)
    ax.set_xticklabels([f"{v:.1f}°E" for v in lon_ticks] if lon_labels else [])
    ax.set_yticklabels([f"{v:.1f}°N" for v in lat_ticks] if lat_labels else [])
    ax.tick_params(labelsize=FS_TICK, length=2.2, width=0.6, pad=1.5, direction="out")
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)


def panel_letter(ax, letter: str, *, x=0.015, y=0.985):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.4},
            zorder=20)


def panel_title(ax, text: str, *, fontsize=FS_TITLE):
    ax.text(0.985, 0.985, text, transform=ax.transAxes, ha="right", va="top",
            fontsize=fontsize, linespacing=1.15,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.4},
            zorder=20)


def scale_bar(ax, frame: dict, length_km: float, *, x0=0.05, y0=0.045):
    (lon0, lon1), (lat0, lat1) = frame["lon"], frame["lat"]
    lat_mid = 0.5 * (lat0 + lat1)
    lon_a = lon0 + x0 * (lon1 - lon0)
    lat_a = lat0 + y0 * (lat1 - lat0)
    lon_b = lon_a + length_km / km_per_deg_lon(lat_mid)
    line = ax.plot([lon_a, lon_b], [lat_a, lat_a], color="black", lw=2.2, solid_capstyle="butt",
                   zorder=15)[0]
    line.set_path_effects(THIN_HALO)
    txt = ax.text(0.5 * (lon_a + lon_b), lat_a + 0.012 * (lat1 - lat0), f"{length_km:.0f} km",
                  ha="center", va="bottom", fontsize=FS_LABEL, fontweight="bold", zorder=15)
    txt.set_path_effects(THIN_HALO)


def inset_colorbar(ax, mappable, label: str, ticks, *, x=0.52, w=0.44, tick_labels=None,
                   extend="neither", bottom_mm=1.2):
    """A horizontal colour bar inside the panel, on a white box that also holds
    the tick numbers below it and the label above it (sizes in millimetres)."""
    _, h_mm = ax._size_mm
    bar_mm, below_mm, above_mm = 1.3, 3.4, 3.4
    y_bar = (bottom_mm + below_mm) / h_mm
    h = bar_mm / h_mm
    pad_x = 0.02
    ax.add_patch(FancyBboxPatch((x - pad_x, bottom_mm / h_mm), w + 2 * pad_x,
                                (below_mm + bar_mm + above_mm) / h_mm,
                                boxstyle="round,pad=0.004", transform=ax.transAxes,
                                facecolor="white", edgecolor="none", alpha=0.88, zorder=16))
    cax = ax.inset_axes([x, y_bar, w, h])
    cbar = ax.figure.colorbar(mappable, cax=cax, orientation="horizontal", ticks=ticks,
                              extend=extend)
    cbar.outline.set_linewidth(0.5)
    cax.tick_params(labelsize=FS_NOTE, length=1.8, width=0.5, pad=1.0)
    if tick_labels is not None:
        cax.set_xticklabels(tick_labels)
    cax.set_zorder(17)
    ax.text(x + 0.5 * w, y_bar + h + 0.8 / h_mm, label, transform=ax.transAxes, ha="center",
            va="bottom", fontsize=FS_LABEL, zorder=18)
    return cbar


def haloed_text(ax, x, y, text, **kwargs):
    kwargs.setdefault("fontsize", FS_LABEL)
    kwargs.setdefault("zorder", 14)
    txt = ax.text(x, y, text, **kwargs)
    txt.set_path_effects(HALO)
    return txt


def draw_2026_geometry(ax, lines: dict, *, lw_dike=1.6, lw_fault=1.0, labels: bool):
    plane = ax.plot(*lines["dike"], color="#c81e1e", lw=0.55 * lw_dike, zorder=12)[0]
    plane.set_path_effects(THIN_HALO)
    opened = ax.plot(*lines["dike_open"], color="#c81e1e", lw=lw_dike, zorder=12)[0]
    opened.set_path_effects(HALO)
    if "fault" in lines:
        fault = ax.plot(*lines["fault"], color="#111111", lw=lw_fault, ls=(0, (3.5, 2)),
                        zorder=12)[0]
        fault.set_path_effects(HALO)
    if labels:
        lon, lat = lines["dike"]
        haloed_text(ax, lon[-1] - 0.012, lat[-1] - 0.004, "2026 dike model",
                    color="#c81e1e", fontsize=FS_MODEL, ha="right", va="top")
        if "fault" in lines:
            flon, flat = lines["fault"]
            haloed_text(ax, flon[0] + 0.03, flat[0] + 0.012, "2026 fault model",
                        color="#111111", fontsize=FS_MODEL, ha="left", va="bottom")


def draw_plane_trace(ax, profiles: dict, *, lw=0.8, label: bool):
    """The 2005 model plane: thin over its full length, thick where the
    September 2005 model opened (the southern 14 km did not)."""
    line = ax.plot(profiles["lon"], profiles["lat"], color="#3a3a3a", lw=0.6 * lw, zorder=11)[0]
    line.set_path_effects(THIN_HALO)
    points = np.column_stack([profiles["lon"], profiles["lat"]])
    segments = np.stack([points[:-1], points[1:]], axis=1)
    opened = (profiles["sept05"][:-1] > 0) & (profiles["sept05"][1:] > 0)
    ax.add_collection(LineCollection(segments[opened], colors="white", linewidths=lw + 1.4,
                                     zorder=11, capstyle="round"))
    ax.add_collection(LineCollection(segments[opened], colors="#3a3a3a", linewidths=lw,
                                     zorder=11, capstyle="round"))
    if label:
        haloed_text(ax, profiles["lon"][-1] + 0.015, profiles["lat"][-1] - 0.03,
                    "2005 dike model", color="#3a3a3a", fontsize=FS_MODEL, ha="left",
                    va="top")


def draw_opening_ribbon(ax, profiles: dict, values: np.ndarray, norm, *, lw=4.4):
    points = np.column_stack([profiles["lon"], profiles["lat"]])
    segments = np.stack([points[:-1], points[1:]], axis=1)
    seg_values = 0.5 * (values[:-1] + values[1:])
    keep = seg_values > 0
    outline = LineCollection(segments[keep], colors="#303030", linewidths=lw + 1.3,
                             zorder=10, capstyle="round")
    ax.add_collection(outline)
    ribbon = LineCollection(segments[keep], cmap=OPENING_CMAP, norm=norm, linewidths=lw,
                            zorder=11, capstyle="round")
    ribbon.set_array(seg_values[keep])
    ax.add_collection(ribbon)
    return ribbon


def draw_volcanoes(ax, *, names: bool, size=38):
    offsets = {"Gabho": (0.02, 0.0, "left", "center"),
               "Dabbahu": (-0.02, 0.0, "right", "center"),
               "Ado'Ale": (-0.02, 0.0, "right", "center"),
               "Manda Hararo": (0.02, 0.0, "left", "center")}
    for name, (lon, lat) in VOLCANOES.items():
        ax.scatter([lon], [lat], marker="^", s=size, facecolor="white", edgecolor="black",
                   linewidth=0.8, zorder=13)
        if names:
            dx, dy, ha, va = offsets[name]
            haloed_text(ax, lon + dx, lat + dy, name, fontsize=7.0, fontweight="bold",
                        ha=ha, va=va)


def draw_earthquake(ax):
    q = EARTHQUAKE
    lat_err = q["err_lat_km"] / KM_PER_DEG_LAT
    lon_err = q["err_lon_km"] / km_per_deg_lon(q["lat"])
    bars = ax.errorbar([q["lon"]], [q["lat"]], xerr=[[lon_err], [lon_err]],
                       yerr=[[lat_err], [lat_err]], fmt="none", ecolor="#202020",
                       elinewidth=0.8, capsize=1.8, capthick=0.8, zorder=13)
    for part in bars[2]:
        part.set_path_effects(THIN_HALO)
    ax.scatter([q["lon"]], [q["lat"]], marker="*", s=75, facecolor="#ffd400",
               edgecolor="black", linewidth=0.7, zorder=14)
    haloed_text(ax, q["lon"] + 0.02, q["lat"] - 0.028, "17 Jul 2026\nearthquake (IGSSA)",
                ha="left", va="top", linespacing=1.1)


def tucked_label(quadrant: str, shafts, half_w: float, half_h: float,
                 clear_mm: float) -> tuple[float, float]:
    """Centre, in mm from the arrow crossing, of a label box with the half
    extents ``half_w`` and ``half_h`` in the quadrant named like a compass
    corner ("NE", "NW", "SE" or "SW"). The box starts at the crossing and
    each shaft (unit direction and half-length, head included) that it sits
    beside pushes it out along the shaft's normal until the box clears the
    shaft axis by ``clear_mm``, so the label stays as near the crossing as
    the clearance allows."""
    sx, sy = (1.0 if "E" in quadrant else -1.0), (1.0 if "N" in quadrant else -1.0)
    lx, ly = sx * half_w, sy * half_h
    for _ in range(4):                  # the pushes interact only through the shafts' skew
        for ux, uy, reach in shafts:
            nx, ny = (-uy, ux) if -uy * sx + ux * sy >= 0.0 else (uy, -ux)
            corners = [(lx + i * half_w, ly + j * half_h) for i in (-1, 1) for j in (-1, 1)]
            along = [cx * ux + cy * uy for cx, cy in corners]
            if max(along) < -reach - clear_mm or min(along) > reach + clear_mm:
                continue
            short = clear_mm - min(cx * nx + cy * ny for cx, cy in corners)
            if short > 0.0:
                lx, ly = lx + short * nx, ly + short * ny
    return lx, ly


def draw_look_arrows(ax, view: dict, *, corner: str, x: float, y: float,
                     flight_mm=4.5, los_mm=3.2, pad_mm=0.6, los_label=None):
    """Flight direction and ground-to-satellite look direction: two arrows
    through one centre on a white box sized from the measured arrows and
    labels so that nothing crosses its edge. ``x``, ``y`` (axes fraction)
    place the box corner named by ``corner`` ("NE", "NW", "SE" or "SW"); the
    arrow half-lengths are in millimetres. Each label sits beyond its arrow
    head, unless ``los_label`` names the quadrant between the arrows into
    which the LOS label is tucked instead: the box is then narrower."""
    w_mm, h_mm = ax._size_mm
    px_per_mm = ax.figure.dpi / 25.4
    renderer = ax.figure.canvas.get_renderer()
    head_mm = 0.7                       # clearance round the shaft ends and heads
    items = []
    for label, colour, angle, half_mm, quadrant in (
            ("flight", "#222222", view["heading_deg"], flight_mm, None),
            ("LOS", "#c81e1e", view["los_deg"], los_mm, los_label)):
        ux, uy = math.sin(math.radians(angle)), math.cos(math.radians(angle))
        text = ax.text(0.0, 0.0, label, transform=ax.transAxes, ha="center", va="center",
                       fontsize=FS_NOTE, color=colour, zorder=18)
        extent = text.get_window_extent(renderer)
        items.append((text, ux, uy, half_mm, colour, quadrant,
                      0.5 * extent.width / px_per_mm, 0.5 * extent.height / px_per_mm))
    shafts = [(ux, uy, half_mm + head_mm) for _, ux, uy, half_mm, *_ in items]
    arrows, labels, xs, ys = [], [], [], []
    for text, ux, uy, half_mm, colour, quadrant, half_w, half_h in items:
        if quadrant is None:
            # The label centre sits beyond the head by the label's own extent
            # along the arrow, so its near edge clears the head by 0.6 mm.
            reach = half_mm + head_mm + 0.6 + abs(ux) * half_w + abs(uy) * half_h
            lx, ly = reach * ux, reach * uy
        else:
            # 0.9 mm from the shaft axis leaves 0.5 mm beside the head barbs.
            lx, ly = tucked_label(quadrant, shafts, half_w, half_h, clear_mm=0.9)
        arrows.append((ux, uy, half_mm, colour))
        labels.append((text, lx, ly))
        for s in (half_mm + head_mm, -(half_mm + head_mm)):
            xs += [s * ux - head_mm, s * ux + head_mm]
            ys += [s * uy - head_mm, s * uy + head_mm]
        xs += [lx - half_w, lx + half_w]
        ys += [ly - half_h, ly + half_h]
    x0, x1 = min(xs) - pad_mm, max(xs) + pad_mm
    y0, y1 = min(ys) - pad_mm, max(ys) + pad_mm
    cx = x - (x1 if "E" in corner else x0) / w_mm
    cy = y - (y1 if "N" in corner else y0) / h_mm
    ax.add_patch(FancyBboxPatch((cx + x0 / w_mm, cy + y0 / h_mm), (x1 - x0) / w_mm,
                                (y1 - y0) / h_mm, boxstyle="round,pad=0.004",
                                transform=ax.transAxes, facecolor="white", edgecolor="none",
                                alpha=0.86, zorder=16))
    for ux, uy, half_mm, colour in arrows:
        dx, dy = half_mm * ux / w_mm, half_mm * uy / h_mm
        ax.annotate("", xy=(cx + dx, cy + dy), xytext=(cx - dx, cy - dy),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops={"arrowstyle": "-|>", "lw": 1.1, "color": colour,
                                "shrinkA": 0, "shrinkB": 0}, zorder=18)
    for text, lx, ly in labels:
        text.set_position((cx + lx / w_mm, cy + ly / h_mm))


def series_axes(ax, year, los, *, title: str, ylabel: str, xlabel: bool, title_left=False):
    ax.scatter(year, los, s=2.2, color="black", linewidths=0, zorder=5)
    ax.set_xlim(2014.7, 2026.7)
    ax.set_xticks([2016, 2020, 2024])
    ax.set_xticklabels(["2016", "2020", "2024"] if xlabel else [])
    ax.tick_params(labelsize=FS_TICK, length=2.0, width=0.5, pad=1.5)
    ax.set_ylabel(ylabel, fontsize=FS_LABEL, labelpad=1.5)
    if xlabel:
        ax.set_xlabel("year", fontsize=FS_LABEL, labelpad=1.0)
    ax.text(0.22 if title_left else 0.97, 0.95, title, transform=ax.transAxes,
            ha="left" if title_left else "right", va="top", fontsize=FS_TITLE, fontweight="bold")
    ax.grid(True, lw=0.3, color="#cccccc", zorder=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)


# ---------------------------------------------------------------------------
# The figure


def build(dem_south: Path, out_stem: Path, coherence_min: float, notes: bool) -> list[Path]:
    plt.rcParams.update({"font.size": 7, "axes.linewidth": 0.7, "pdf.fonttype": 42})
    dem_paths = [DEM_PATH, dem_south]

    shade_a = terrain_shade(FRAME_A, dem_paths)
    shade_b = terrain_shade(FRAME_B, dem_paths)
    _, _, extent_a = frame_grid(FRAME_A)
    _, _, extent_b = frame_grid(FRAME_B)
    lines_2026 = dyke_2026_lines()
    profiles = hamling_profiles()
    series = site_series(SERIES_CSV)
    cyclic = cmc.romaO

    # Layout in millimetres: (a) at the left, the context block at the right.
    left = 11.0
    a_w = 81.0
    a_h = a_w * (FRAME_A["lat"][1] - FRAME_A["lat"][0]) / (
        (FRAME_A["lon"][1] - FRAME_A["lon"][0]) * math.cos(math.radians(12.24)))
    b_w = 36.0
    b_h = b_w * (FRAME_B["lat"][1] - FRAME_B["lat"][0]) / (
        (FRAME_B["lon"][1] - FRAME_B["lon"][0]) * math.cos(math.radians(12.38)))
    y_base = 20.5 if notes else 9.5          # the notes block needs 11 mm below the panels
    top = y_base + a_h
    col1, col2 = 103.0, 149.0
    fig_h = top + 2.0
    fig = plt.figure(figsize=(190 / 25.4, fig_h / 25.4))

    # --- (a) NISAR 2026 ---------------------------------------------------
    ax_a = mm_axes(fig, left, y_base, a_w, a_h)
    wrapped, north_drawn = nisar_rewrapped_cm(FRAME_A, coherence_min)
    norm_a = Normalize(0, NISAR_FRINGE_CM)
    ax_a.imshow(blend(wrapped, shade_a, cyclic, norm_a, alpha=0.90), extent=extent_a,
                origin="upper", interpolation="nearest", zorder=1)
    map_axes(ax_a, FRAME_A, lon_labels=True, lat_labels=True)
    draw_plane_trace(ax_a, profiles, label=True)
    draw_2026_geometry(ax_a, lines_2026, labels=True)
    draw_volcanoes(ax_a, names=True)
    draw_earthquake(ax_a)
    for name, lon, lat in (("Dabbahu MS", 40.30, 12.46), ("Hararo MS", 40.40, 12.02)):
        haloed_text(ax_a, lon, lat, name, fontsize=7.4, fontstyle="italic", ha="left",
                    va="center", color="#202020")
    draw_look_arrows(ax_a, NISAR_VIEW, corner="NE", x=0.975, y=0.945, los_label="NE")
    scale_bar(ax_a, FRAME_A, 20)
    inset_colorbar(ax_a, plt.cm.ScalarMappable(norm=norm_a, cmap=cyclic),
                   "LOS toward satellite, one fringe (cm)",
                   ticks=[0, NISAR_FRINGE_CM / 2, NISAR_FRINGE_CM],
                   tick_labels=["0", f"{NISAR_FRINGE_CM / 2:.1f}", f"{NISAR_FRINGE_CM:.1f}"],
                   x=0.50, w=0.46)
    panel_letter(ax_a, "a")
    panel_title(ax_a, NISAR_LABEL)

    # --- (b) the September 2005 dyke in the track-49 pair --------------------
    norm_env = Normalize(0, ENVISAT_FRINGE_CM)
    ax_b = mm_axes(fig, col1, top - b_h, b_w, b_h)
    b_drawn = ENVISAT_2005.is_file()
    if b_drawn:
        field_b, view_b = envisat_2005_wrapped_cm(FRAME_B)
        ax_b.imshow(blend(field_b, shade_b, cyclic, norm_env, alpha=0.90), extent=extent_b,
                    origin="upper", interpolation="nearest", zorder=1)
    else:
        ax_b.imshow(np.repeat(shade_b[:, :, None], 3, axis=2), extent=extent_b, origin="upper",
                    interpolation="nearest", zorder=1)
    map_axes(ax_b, FRAME_B, lon_labels=True, lat_labels=True)
    draw_plane_trace(ax_b, profiles, lw=0.8 if b_drawn else 1.2, label=False)
    draw_2026_geometry(ax_b, lines_2026, lw_dike=1.0, labels=False)
    draw_volcanoes(ax_b, names=False, size=16)
    if b_drawn:
        draw_look_arrows(ax_b, view_b, corner="NE", x=0.975, y=0.81,
                         flight_mm=3.4, los_mm=2.4, los_label="NW")
        inset_colorbar(ax_b, plt.cm.ScalarMappable(norm=norm_env, cmap=cyclic),
                       "one fringe (cm)", ticks=[0, ENVISAT_FRINGE_CM],
                       tick_labels=["0", f"{ENVISAT_FRINGE_CM:.1f}"], x=0.44, w=0.50)
        title_b = "Sept 2005 dike\nEnvisat track 49"
    else:
        haloed_text(ax_b, 0.96, 0.30, "Envisat interferogram\nof the Sept 2005 dike\n(Wright et al., 2006)\nto be supplied",
                    transform=ax_b.transAxes, ha="right", va="center", fontsize=FS_LABEL,
                    linespacing=1.25, color="#202020")
        haloed_text(ax_b, profiles["lon"][0] + 0.02, profiles["lat"][0] - 0.01,
                    "2005 dike model", color="#3a3a3a", fontsize=FS_MODEL, ha="left",
                    va="top")
        title_b = "Sept 2005 dike"
    scale_bar(ax_b, FRAME_B, 20)
    panel_letter(ax_b, "b")
    panel_title(ax_b, title_b, fontsize=FS_NOTE)

    # --- (c) Envisat track-300 interferogram of the October 2008 dyke ---------
    grid_name, title = ENVISAT_PANELS["c"]
    ax_c = mm_axes(fig, col2, top - b_h, b_w, b_h)
    field = envisat_wrapped_cm(FRAME_B, grid_name)
    ax_c.imshow(blend(field, shade_b, cyclic, norm_env, alpha=0.90), extent=extent_b,
                origin="upper", interpolation="nearest", zorder=1)
    map_axes(ax_c, FRAME_B, lon_labels=True, lat_labels=False)
    draw_2026_geometry(ax_c, lines_2026, lw_dike=1.0, labels=False)
    draw_volcanoes(ax_c, names=False, size=16)
    draw_look_arrows(ax_c, ENVISAT_VIEW, corner="NE", x=0.975, y=0.81,
                     flight_mm=3.4, los_mm=2.4, los_label="NE")
    scale_bar(ax_c, FRAME_B, 20)
    inset_colorbar(ax_c, plt.cm.ScalarMappable(norm=norm_env, cmap=cyclic),
                   "one fringe (cm)", ticks=[0, ENVISAT_FRINGE_CM],
                   tick_labels=["0", f"{ENVISAT_FRINGE_CM:.1f}"], x=0.44, w=0.50)
    panel_letter(ax_c, "c")
    panel_title(ax_c, title, fontsize=FS_NOTE)

    # --- (d) Sentinel-1 cumulative 2014 to 2026, wrapped ---------------------
    ax_d = mm_axes(fig, col1, y_base, b_w, b_h)
    cum_wrapped = s1_wrapped_cm(FRAME_B)
    norm_d = Normalize(0, S1_FRINGE_MM / 10)
    ax_d.imshow(blend(cum_wrapped, shade_b, cyclic, norm_d, alpha=0.90), extent=extent_b,
                origin="upper", interpolation="nearest", zorder=1)
    map_axes(ax_d, FRAME_B, lon_labels=True, lat_labels=True)
    draw_plane_trace(ax_d, profiles, lw=0.6, label=False)
    draw_2026_geometry(ax_d, lines_2026, lw_dike=1.0, labels=False)
    draw_volcanoes(ax_d, names=False, size=16)
    for letter, site in SITES.items():
        ax_d.scatter([site["lon"]], [site["lat"]], s=22, facecolor="white", edgecolor="black",
                     linewidth=0.8, zorder=13)
        haloed_text(ax_d, site["lon"] + 0.018, site["lat"] + 0.004, letter, fontsize=6.8,
                    fontweight="bold", ha="left", va="center")
    ax_d.scatter([FLANK_REFERENCE[0]], [FLANK_REFERENCE[1]], s=16, facecolor="none",
                 edgecolor="black", linewidth=0.8, zorder=13)
    haloed_text(ax_d, FLANK_REFERENCE[0] + 0.016, FLANK_REFERENCE[1] + 0.004, "f′",
                fontsize=6.8, fontweight="bold", ha="left", va="center")
    draw_look_arrows(ax_d, S1_VIEW, corner="SW", x=0.025, y=0.25, flight_mm=3.4, los_mm=2.4,
                     los_label="NW")
    scale_bar(ax_d, FRAME_B, 20)
    inset_colorbar(ax_d, plt.cm.ScalarMappable(norm=norm_d, cmap=cyclic),
                   "one fringe (cm)", ticks=[0, S1_FRINGE_MM / 10],
                   tick_labels=["0", f"{S1_FRINGE_MM / 10:.1f}"], x=0.44, w=0.50)
    panel_letter(ax_d, "d")
    panel_title(ax_d, "Sentinel-1 descending\n23 Oct 2014 – 11 Jul 2026\ncumulative, wrapped",
                fontsize=FS_NOTE)

    # --- (e), (f) the two series ---------------------------------------------
    gap = 7.0
    h_each = (b_h - gap) / 2
    ax_e = mm_axes(fig, col2 + 6.0, y_base + h_each + gap, b_w - 6.0, h_each)
    year, los = series[SITES["e"]["series"]]
    series_axes(ax_e, year, los, title="Ado'Ale", ylabel="LOS (mm)", xlabel=False,
                title_left=True)
    panel_letter(ax_e, "e")
    ax_f = mm_axes(fig, col2 + 6.0, y_base, b_w - 6.0, h_each)
    year, los = series[SITES["f"]["series"]]
    series_axes(ax_f, year, los, title="graben − flank", ylabel="f − f′ (mm)",
                xlabel=True)
    panel_letter(ax_f, "f")

    if notes:
        draft_notes(fig, fig_h, left, coherence_min, north_drawn, b_drawn)

    outputs = []
    for suffix, kwargs in ((".png", {"dpi": 400}), (".pdf", {"dpi": 400})):
        path = out_stem.with_suffix(suffix)
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def draft_notes(fig, fig_h: float, left: float, coherence_min: float, north_drawn: bool,
                b_drawn: bool) -> None:
    """The notes to the coauthors under the draft version of the figure."""
    north_note = (f"north of {NISAR_CROP_NORTH_LAT}°N the panel shows the display strip of "
                  "the same pair (same corrections, not unwrapped)" if north_drawn else
                  f"the analysis crop ends at {NISAR_CROP_NORTH_LAT}°N; the display strip "
                  "north of it is not built")
    b_note = ("(b) the 6 May – 28 Oct 2005 pair of track 49 (the pair of Wright et al., 2006, Fig. 2a) "
              "processed here; one colour cycle is 2.8 cm of range change, as in (c)." if b_drawn else
              "(b) the delivered track-300 grids are the 2006–2010 dikes (the \"sept\" grid is the "
              "September 2006 dike); the 2005 interferogram is awaited.")
    fig.text(left / 190.0, 0.8 / fig_h,
             f"Draft. (a) samples below coherence {coherence_min:.2f} are omitted; {north_note}.\n"
             f"{b_note}\n"
             "(c) Hamling et al. (2010), acquisition dates to be confirmed; one colour cycle is "
             "2.8 cm of range change.\n"
             "(d) the cumulative map wrapped at one Sentinel-1 fringe (2.8 cm), as in B. Ireland "
             "(2026, unpublished figure); blank where the neighbourhood phase is inconsistent.\n"
             "(e), (f): digitized from B. Ireland (2026, unpublished figure); relative to 23 Oct 2014 "
             "in the original processing reference.",
             fontsize=FS_NOTE, ha="left", va="bottom", color="#333333", linespacing=1.25)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dem-south", type=Path, default=DEM_SOUTH_PATH,
                        help="Copernicus GLO-30 tile N11 E040 (ignored data, fetched from "
                             f"{DEM_SOUTH_URL})")
    parser.add_argument("--notes", action="store_true",
                        help="add the notes to the coauthors under the panels (the draft version, "
                             "written to results/figure1_draft by default)")
    parser.add_argument("--out", type=Path,
                        help="output stem (PNG and PDF are written); results/figure1 by default")
    parser.add_argument("--coherence-min", type=float, default=COHERENCE_MIN,
                        help="NISAR samples below this coherence are not displayed")
    args = parser.parse_args()
    out = args.out or RESULTS / ("figure1_draft" if args.notes else "figure1")
    for path in build(args.dem_south, out, args.coherence_min, args.notes):
        print(path)


if __name__ == "__main__":
    main()
