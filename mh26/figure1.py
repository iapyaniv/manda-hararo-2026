"""Draft of the GRL location figure (Figure 1), in the panel layout agreed with the coauthors.

Panels
  (a) the July 2026 intrusion in the NISAR descending interferogram of 7 to
      19 July (165_082_D, 30-m GSLC chain, a 90 km processing crop that ends
      at 12.59 N): the corrected unwrapped LOS field the paper models
      (ionosphere screen subtracted, far-field plane removed), screened at the
      supplement's coherence threshold, smoothed by one pixel in the complex
      domain and rewrapped at one L-band fringe; with the model dyke plane
      (the part with opening above one third of the peak drawn thick), the
      candidate graben fault, the 2005 dyke (model plane, thick where it
      opened), the one located July 2026 earthquake (two agency locations
      with their stated errors), the volcanic centres at their mapped
      positions and the segment names;
  (b) the September 2005 dyke and (c) the October 2008 dyke in Ian Hamling's
      wrapped Envisat track-300 interferograms (Hamling et al. 2009, 2010;
      acquisition dates to be confirmed), multilooked as complex phasors;
  (d) Sentinel-1 descending cumulative LOS displacement 23 October 2014 to
      11 July 2026 (P. Espín Bedón's product, re-referenced; see
      data/sentinel1-cumulative/README), with the time-series sites of (e)
      and (f) and the Sentinel-1 look direction;
  (e, f) the Ado'Ale series and the graben-minus-flank series of B. Ireland's
      Sentinel-1 time-series figure, replotted from values digitized from it
      (data/ireland-series; DRAFT: to be replaced by the series themselves).

Everything is drawn in geographic coordinates with one Copernicus GLO-30
hillshade (tile N12 E040 from the repository and N11 E040, fetched, for the
southern strip). Nothing here is an analysis quantity: the figure only
displays saved products and the delivered models.

Run from the repository root:

  python -m mh26.figure1 [--dem-south <N11 tile>]
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
from matplotlib.colors import LightSource, Normalize, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
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
NISAR_LABEL = "NISAR, 7 – 19 Jul 2026"
S1_PRE = (REPO / "data" / "sentinel1-cumulative" / "Pre"
          / "pre_079D_07694_131313_20141023_20260711.cum.geo.tif")
NISAR_FRINGE_M = WAVELENGTH_L / 2
NISAR_FRINGE_CM = 100 * NISAR_FRINGE_M
# Ian Hamling's wrapped Envisat track-300 interferograms of the earlier dykes
# (one GMT grid per dyke, wrapped LOS in centimetres over one Envisat fringe).
ENVISAT_DIR = REPO / "data" / "ihamling" / "opening_files"
ENVISAT_FRINGE_CM = 2.8118
ENVISAT_PANELS = {"b": ("sept", "Sept 2005 dyke\nEnvisat track 300"),
                  "c": ("oct", "Oct 2008 dyke\nEnvisat track 300")}
# Envisat ascending, right-looking: heading about 348 deg, so the
# ground-to-satellite azimuth is about 258 deg.
ENVISAT_VIEW = {"heading_deg": 348.0, "los_deg": 258.0}
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
# B. Ireland's six site series as digitized from his figure; see
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
# The one located July 2026 earthquake, as located by two agencies
# (compiled by D. Keir). Errors in km (latitude, longitude), drawn as plain
# error bars: no joint confidence level is implied.
EARTHQUAKES = {
    "IGSSA": {"lon": 40.853, "lat": 12.042, "err_lat_km": 3.8, "err_lon_km": 3.6},
    "ISC": {"lon": 40.6394, "lat": 11.7869, "err_lat_km": 27.8, "err_lon_km": 15.7},
}
# B. Ireland's time-series sites (positions read from his figure and
# confirmed against the cumulative map). His graben series is the graben
# site minus the flank site 8.5 km to the north-east.
SITES = {
    "e": {"name": "Ado'Ale", "lon": 40.698, "lat": 12.330, "series": "Ado'Ale"},
    "f": {"name": "graben", "lon": 40.584, "lat": 12.393, "series": "Graben subsidence"},
}
FLANK_REFERENCE = (40.649, 12.436)
S1_REFERENCE_AREA = ((40.93682, 40.94682), (12.63819, 12.64819))

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


def nisar_rewrapped_cm(frame: dict, coherence_min: float) -> np.ndarray:
    """The corrected unwrapped LOS field (metres, positive toward the
    satellite), screened by coherence, then wrapped to one fringe and
    smoothed by one pixel and resampled as a complex phasor, so that
    integer-cycle steps in the unwrapped field cannot leak into fractional
    fringes."""
    data = np.load(NISAR_UNW, allow_pickle=False)
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


def s1_cumulative_mm(frame: dict) -> np.ndarray:
    """The delivered cumulative map, as delivered: its own gaps stay gaps and
    nothing else is masked. The south-western corner (west of the rift) holds
    about -100 mm with 25 to 35 mm pixel scatter in this track and +350 mm in
    the ascending one; whether that is plate motion relative to the reference
    area or ground to be masked is the processors' call."""
    transform, shape, _ = frame_grid(frame)
    with rasterio.open(S1_PRE) as dataset:
        source = dataset.read(1, masked=True).filled(np.nan)
        return warp(source, dataset.transform, dataset.crs, transform, shape,
                    Resampling.bilinear)


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


def envisat_wrapped_cm(frame: dict, name: str, looks: int = 3) -> np.ndarray:
    """One of the track-300 interferograms, multilooked as a complex phasor
    over `looks` x `looks` cells and resampled the same way, in centimetres of
    one fringe."""
    z, src_transform = read_grd(ENVISAT_DIR / f"t300_{name}_wrp_dyke.grd")
    ok = np.isfinite(z)
    if not ok.any() or np.nanmax(z) == np.nanmin(z):
        raise ValueError(f"track-300 grid {name!r} holds no interferogram")
    phase = 2 * np.pi * np.where(ok, z, 0.0) / ENVISAT_FRINGE_CM
    weight = uniform_filter(ok.astype(np.float64), looks)
    real = uniform_filter(np.where(ok, np.cos(phase), 0.0), looks)
    imag = uniform_filter(np.where(ok, np.sin(phase), 0.0), looks)
    keep = ok & (weight > 0.5)
    real = np.where(keep, real / np.maximum(weight, 1e-6), np.nan)
    imag = np.where(keep, imag / np.maximum(weight, 1e-6), np.nan)
    transform, shape, _ = frame_grid(frame)
    real_w = warp(real, src_transform, UTM_CRS, transform, shape, Resampling.average)
    imag_w = warp(imag, src_transform, UTM_CRS, transform, shape, Resampling.average)
    wrapped = np.mod(np.arctan2(imag_w, real_w), 2 * np.pi) / (2 * np.pi) * ENVISAT_FRINGE_CM
    wrapped[~np.isfinite(real_w) | ~np.isfinite(imag_w)] = np.nan
    return wrapped


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
        haloed_text(ax, lon[-1] - 0.012, lat[-1] - 0.004, "2026 dyke (model)",
                    color="#c81e1e", fontsize=6.8, fontweight="bold", ha="right", va="top")
        if "fault" in lines:
            flon, flat = lines["fault"]
            haloed_text(ax, flon[0] + 0.03, flat[0] + 0.012, "candidate fault",
                        color="#111111", fontsize=FS_LABEL, ha="left", va="bottom")


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
                    "2005 dyke (model)", color="#3a3a3a", fontsize=FS_LABEL, ha="left",
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


def draw_earthquakes(ax, frame: dict):
    lat0 = frame["lat"][0]
    for agency, q in EARTHQUAKES.items():
        lat_err = q["err_lat_km"] / KM_PER_DEG_LAT
        lon_err = q["err_lon_km"] / km_per_deg_lon(q["lat"])
        bars = ax.errorbar([q["lon"]], [q["lat"]], xerr=[[lon_err], [lon_err]],
                           yerr=[[lat_err], [lat_err]], fmt="none", ecolor="#202020",
                           elinewidth=0.8, capsize=1.8, capthick=0.8, zorder=13)
        for part in bars[2]:
            part.set_path_effects(THIN_HALO)
        ax.scatter([q["lon"]], [q["lat"]], marker="*", s=75, facecolor="#ffd400",
                   edgecolor="black", linewidth=0.7, zorder=14)
        if agency == "ISC":
            haloed_text(ax, q["lon"] - 0.02, q["lat"] + 0.03, "ISC (same event)",
                        ha="right", va="bottom", zorder=19)
            south = q["lat"] - lat_err
            if south < lat0:
                haloed_text(ax, q["lon"] - 0.02, lat0 + 0.004,
                            f"error bar continues to {south:.2f}°N", fontsize=FS_NOTE,
                            ha="right", va="bottom", zorder=19)
        else:
            haloed_text(ax, q["lon"] + 0.02, q["lat"] - 0.028, "17 Jul 2026\nearthquake (IGSSA)",
                        ha="left", va="top", linespacing=1.1)


def draw_look_arrows(ax, view: dict, *, x, y, w, h, flight_mm=6.0, los_mm=4.2):
    """Flight direction and ground-to-satellite LOS azimuth in a white box."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004",
                                transform=ax.transAxes, facecolor="white", edgecolor="none",
                                alpha=0.86, zorder=16))
    w_mm, h_mm = ax._size_mm
    cx, cy = x + 0.5 * w, y + 0.5 * h
    for angle, length_mm, label, colour in ((view["heading_deg"], flight_mm, "flight", "#222222"),
                                            (view["los_deg"], los_mm, "LOS", "#c81e1e")):
        dx = length_mm * math.sin(math.radians(angle)) / w_mm
        dy = length_mm * math.cos(math.radians(angle)) / h_mm
        ax.annotate("", xy=(cx + dx, cy + dy), xytext=(cx - dx, cy - dy),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops={"arrowstyle": "-|>", "lw": 1.1, "color": colour,
                                "shrinkA": 0, "shrinkB": 0}, zorder=18)
        ax.text(cx + 1.35 * dx, cy + 1.35 * dy, label, transform=ax.transAxes, ha="center",
                va="center", fontsize=FS_NOTE, color=colour, zorder=18)


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


def build(dem_south: Path, out_stem: Path, coherence_min: float) -> list[Path]:
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
    y_base = 13.0
    top = y_base + a_h
    col1, col2 = 103.0, 149.0
    fig_h = top + 2.0
    fig = plt.figure(figsize=(190 / 25.4, fig_h / 25.4))

    # --- (a) NISAR 2026 ---------------------------------------------------
    ax_a = mm_axes(fig, left, y_base, a_w, a_h)
    wrapped = nisar_rewrapped_cm(FRAME_A, coherence_min)
    norm_a = Normalize(0, NISAR_FRINGE_CM)
    ax_a.imshow(blend(wrapped, shade_a, cyclic, norm_a, alpha=0.90), extent=extent_a,
                origin="upper", interpolation="nearest", zorder=1)
    map_axes(ax_a, FRAME_A, lon_labels=True, lat_labels=True)
    draw_plane_trace(ax_a, profiles, label=True)
    draw_2026_geometry(ax_a, lines_2026, labels=True)
    draw_volcanoes(ax_a, names=True)
    draw_earthquakes(ax_a, FRAME_A)
    for name, lon, lat in (("Dabbahu MS", 40.30, 12.46), ("Hararo MS", 40.40, 12.02)):
        haloed_text(ax_a, lon, lat, name, fontsize=7.4, fontstyle="italic", ha="left",
                    va="center", color="#202020")
    draw_look_arrows(ax_a, NISAR_VIEW, x=0.77, y=0.79, w=0.20, h=0.15)
    scale_bar(ax_a, FRAME_A, 20)
    inset_colorbar(ax_a, plt.cm.ScalarMappable(norm=norm_a, cmap=cyclic),
                   "LOS toward satellite, one fringe (cm)",
                   ticks=[0, NISAR_FRINGE_CM / 2, NISAR_FRINGE_CM],
                   tick_labels=["0", f"{NISAR_FRINGE_CM / 2:.1f}", f"{NISAR_FRINGE_CM:.1f}"],
                   x=0.50, w=0.46)
    panel_letter(ax_a, "a")
    panel_title(ax_a, NISAR_LABEL)

    # --- (b), (c) Envisat track-300 interferograms of two earlier dykes -------
    norm_env = Normalize(0, ENVISAT_FRINGE_CM)
    for letter, ax_x, lat_labels in (("b", col1, True), ("c", col2, False)):
        grid_name, title = ENVISAT_PANELS[letter]
        ax = mm_axes(fig, ax_x, top - b_h, b_w, b_h)
        field = envisat_wrapped_cm(FRAME_B, grid_name)
        ax.imshow(blend(field, shade_b, cyclic, norm_env, alpha=0.90), extent=extent_b,
                  origin="upper", interpolation="nearest", zorder=1)
        map_axes(ax, FRAME_B, lon_labels=True, lat_labels=lat_labels)
        draw_2026_geometry(ax, lines_2026, lw_dike=1.0, labels=False)
        draw_volcanoes(ax, names=False, size=16)
        if letter == "b":
            draw_look_arrows(ax, ENVISAT_VIEW, x=0.67, y=0.60, w=0.30, h=0.19,
                             flight_mm=4.6, los_mm=3.2)
        scale_bar(ax, FRAME_B, 20)
        inset_colorbar(ax, plt.cm.ScalarMappable(norm=norm_env, cmap=cyclic),
                       "one fringe (cm)", ticks=[0, ENVISAT_FRINGE_CM],
                       tick_labels=["0", f"{ENVISAT_FRINGE_CM:.1f}"], x=0.44, w=0.50)
        panel_letter(ax, letter)
        panel_title(ax, title, fontsize=FS_NOTE)

    # --- (d) Sentinel-1 cumulative 2014 to 2026 -----------------------------
    ax_d = mm_axes(fig, col1, y_base, b_w, b_h)
    cum_mm = s1_cumulative_mm(FRAME_B)
    norm_d = TwoSlopeNorm(vmin=-100, vcenter=0, vmax=300)
    ax_d.imshow(blend(cum_mm, shade_b, cmc.vik, norm_d, alpha=0.88), extent=extent_b,
                origin="upper", interpolation="nearest", zorder=1)
    map_axes(ax_d, FRAME_B, lon_labels=True, lat_labels=True)
    draw_plane_trace(ax_d, profiles, lw=0.6, label=False)
    draw_2026_geometry(ax_d, lines_2026, lw_dike=1.0, labels=False)
    draw_volcanoes(ax_d, names=False, size=16)
    (ref_lon, ref_lat) = S1_REFERENCE_AREA
    ax_d.add_patch(Rectangle((ref_lon[0], ref_lat[0]), ref_lon[1] - ref_lon[0],
                             ref_lat[1] - ref_lat[0], facecolor="none", edgecolor="black",
                             lw=0.8, zorder=12))
    haloed_text(ax_d, ref_lon[0] - 0.012, ref_lat[0] + 0.005, "ref", fontsize=FS_NOTE,
                ha="right", va="center")
    for letter, site in SITES.items():
        ax_d.scatter([site["lon"]], [site["lat"]], s=22, facecolor="white", edgecolor="black",
                     linewidth=0.8, zorder=13)
        haloed_text(ax_d, site["lon"] + 0.018, site["lat"] + 0.004, letter, fontsize=6.8,
                    fontweight="bold", ha="left", va="center")
    ax_d.scatter([FLANK_REFERENCE[0]], [FLANK_REFERENCE[1]], s=16, facecolor="none",
                 edgecolor="black", linewidth=0.8, zorder=13)
    haloed_text(ax_d, FLANK_REFERENCE[0] + 0.016, FLANK_REFERENCE[1] + 0.004, "f′",
                fontsize=6.8, fontweight="bold", ha="left", va="center")
    draw_look_arrows(ax_d, S1_VIEW, x=0.03, y=0.26, w=0.30, h=0.19, flight_mm=4.6, los_mm=3.2)
    scale_bar(ax_d, FRAME_B, 20)
    inset_colorbar(ax_d, plt.cm.ScalarMappable(norm=norm_d, cmap=cmc.vik),
                   "LOS (mm)", ticks=[-100, 0, 150, 300], x=0.44, w=0.50,
                   extend="min")
    panel_letter(ax_d, "d")
    panel_title(ax_d, "Sentinel-1 descending\n23 Oct 2014 – 11 Jul 2026", fontsize=FS_NOTE)

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

    fig.text(left / 190.0, 0.8 / fig_h,
             "Draft. (a) shows the 90 km NISAR processing crop, which ends at 12.59°N; "
             f"samples below coherence {coherence_min:.2f} are omitted.\n"
             "(b), (c): Envisat track-300 interferograms of Hamling et al. (2009, 2010), "
             "acquisition dates to be confirmed; one colour cycle is 2.8 cm of range change.\n"
             "(e), (f): digitized from B. Ireland's figure, pending his series; relative to "
             "23 Oct 2014 in the original processing reference, not the reference area of (d).",
             fontsize=FS_NOTE, ha="left", va="bottom", color="#333333", linespacing=1.25)

    outputs = []
    for suffix, kwargs in ((".png", {"dpi": 400}), (".pdf", {"dpi": 400})):
        path = out_stem.with_suffix(suffix)
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dem-south", type=Path, default=DEM_SOUTH_PATH,
                        help="Copernicus GLO-30 tile N11 E040 (ignored data, fetched from "
                             f"{DEM_SOUTH_URL})")
    parser.add_argument("--out", type=Path, default=RESULTS / "figure1_draft",
                        help="output stem (PNG and PDF are written)")
    parser.add_argument("--coherence-min", type=float, default=COHERENCE_MIN,
                        help="NISAR samples below this coherence are not displayed")
    args = parser.parse_args()
    for path in build(args.dem_south, args.out, args.coherence_min):
        print(path)


if __name__ == "__main__":
    main()
