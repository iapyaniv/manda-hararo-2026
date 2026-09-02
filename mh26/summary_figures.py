"""Publication-style summary figures from the saved analysis results.

Eight figures in the presentation of Xu et al. (2026,
doi:10.1016/j.geog.2026.06.002). Every figure reads saved results and data
files; nothing is refitted.

1. ``summary_los_model_residual``: rewrapped observed, model, and residual
   LOS maps with rift-perpendicular profiles for all four co-event
   geometries, from the trace-supported distributed fit to the reprocessed
   30-m (GSLC) inputs (``results/distributed_fit_mapped_gslc.json``).
2. ``summary_decomposition_enu``: the free East/North/Up solve on 1.28-km
   cells with formal uncertainties, geometry coverage, and four-geometry
   closure (recomputed and checked against
   ``results/decomposition_joint_gslc.json``).
3. ``summary_source_geometry``: the saved dike-opening and east-fault slip
   distributions in physical geometry (3-D view, plan view, cross-section,
   and the unfolded fault plane).
4. ``summary_phase_gradient``: unwrapped LOS and the signed
   rift-perpendicular wrapped-phase difference (after Sandwell & Price,
   1998, doi:10.1029/1998JB900008; cf. Xu et al., 2020,
   doi:10.1126/science.abd1690) with the mapped ridges and cross-rift
   profiles.
5. ``summary_fault_stress``: dike-induced normal and Coulomb stress change
   on the trace-supported east-fault receivers, with map-space fields at a
   fixed depth and a GSLC/GUNW input comparison on identical receivers.
6. ``summary_decomposition_dike_frame``: the four decay-normalized LOS
   inputs (top row) and the same free solve rotated into the dike frame —
   dike-perpendicular, dike-parallel (positive near-north, masked above
   5 cm formal σ), and up — with the rotated formal uncertainties
   (``mh26.decompose3d.dike_frame``).
7. ``summary_los_model_residual_wright2006``: the observed/model/residual
   view for the fit constrained to the Wright et al. (2006) 2005 dike
   plane (``results/distributed_fit_geometry_<hash>_gslc.json``). The
   segmented plane has no single straight frame, so its map trace and
   profile frame come from the saved patch arrays.
8. ``summary_uperp_wright2006``: the opening cross-section on that plane
   (distance along the 2005 plane versus depth, after Hamling et al.,
   2010, Fig. 2) from the same fit — the publication-grade version of
   ``results/uperp_2026_wright2006.png``.

Inputs: the gslc data archive (or a rebuild, see ``data/README.md``), the
saved results listed above, and one Copernicus GLO-30 terrain tile used
only for the map background and profiles:

    https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N12_00_E040_00_DEM/Copernicus_DSM_COG_10_N12_00_E040_00_DEM.tif

Place the tile at ``data/dem/Copernicus_DSM_COG_10_N12_00_E040_00_DEM.tif``.
Its SHA-256 hash is pinned below and verified at run time.

Run after ``mh26.distributed --mapped --gslc``, ``mh26.distributed
--geometry results/geometry_wright2006.json --gslc``, ``mh26.decompose3d
--gslc``, ``mh26.faults``, and ``mh26.stress``:

    python -m mh26.summary_figures              # all eight figures
    python -m mh26.summary_figures --only los   # one figure

Writes ``results/summary_*.png`` (300 dpi) with matching vector PDF files,
and, on a full run, ``results/summary_figures.json`` with the input and
output identity hashes. Raw data files are covered by the archive checksum
manifests in ``data/`` and are not re-hashed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import string
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from cmcrameri import cm as cmc
from matplotlib.cm import ScalarMappable
from matplotlib.colors import (
    BoundaryNorm,
    LightSource,
    ListedColormap,
    Normalize,
    TwoSlopeNorm,
)
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MultipleLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from pyrocko.modelling import okada_ext
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.warp import reproject
from rasterio.warp import transform as warp_transform
import rasterio
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

from .decompose3d import dike_frame, lattice_axes, load_cells, solve_cells_3d
from .distributed import (
    NISAR_LABEL,
    S1_GEOMETRY,
    _strike_vectors,
    geometry_from_json,
    geometry_sha256,
)
from .gslc import DATA_GSLC, UNWRAP_COHERENCE_MIN
from .gunw import load_pair, los_unit_neu
from .params import (
    CENTER_E,
    CENTER_N,
    EPSG_UTM37N,
    INVERSION_PAIRS,
    LAME_LAMBDA,
    REPO,
    RESULTS,
    S1_TRACK_PAIRS,
    SHEAR_MODULUS,
    WAVELENGTH_L,
)
from .s1_tracks import (
    MIN_COEVENT_ABS_R,
    coevent_pixel_stats,
    load_frame,
    load_pair_los,
    reference,
)
from .stress import (
    DIKE_STANDOFF_M,
    FRICTION,
    MAIN_DIP,
    PROFILE_DZ_M,
    RECEIVER_DIPS,
    fault_patch_receivers,
    fault_profile_receivers,
    resolve_normal_slip,
    selftest,
    side_stats,
    stress_at,
)


# Terrain tile (map background and profile topography only; no analysis
# quantity depends on it).
DEM_PATH = REPO / "data" / "dem" / "Copernicus_DSM_COG_10_N12_00_E040_00_DEM.tif"
DEM_URL = (
    "https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_N12_00_E040_00_DEM/"
    "Copernicus_DSM_COG_10_N12_00_E040_00_DEM.tif"
)
DEM_SHA256 = "90ad0f4b67afff380496e9e35a0dd74ea4dcb915ab0e64655b28b30cf2da6d24"

# Common display frame (UTM 37N metres) and model-evaluation posting. The
# 200-m posting is display only; the inversion support stays at the
# ~1.3-km blocks reported with the fit.
X_MIN, X_MAX = 660_000.0, 701_000.0
Y_MIN, Y_MAX = 1_329_000.0, 1_368_000.0
DISPLAY_STEP_M = 200.0
UTM_CRS = f"EPSG:{EPSG_UTM37N}"

# LOS figure display.
DISPLAY_CYCLE_M = 0.10
RESIDUAL_LIMIT_M = 0.10

# Decomposition figure display.
DISPLACEMENT_LIMIT_CM = 60.0
SIGMA_LIMIT_CM = 8.0
CLOSURE_LIMIT_SIGMA = 2.0
# Display stipple only: cells of the weakly observed horizontal component
# (north, or dike-parallel) above this formal sigma are dotted so a clipped
# near-axis cell cannot read as signal. The solve is unchanged.
WEAK_SIGMA_STIPPLE_CM = 6.0

# Source figure display.
DIKE_LIMIT_M = 2.6
FAULT_LIMIT_M = 0.65
# batlowW_r starts at paper-white, so low fault slip vanishes on a white
# page; the truncated ramp keeps zero slip a visible pale tint.
FAULT_CMAP = ListedColormap(
    cmc.batlowW_r(np.linspace(0.14, 1.0, 256)), name="batlowW_r_floor"
)

# Phase-gradient figure display.
COHERENCE_MIN_GRADIENT = 0.30
COHERENCE_MIN_LOS = 0.20
LOS_LIMIT_CM = 70.0
PHASE_INCREMENT_LIMIT_RAD = 0.30
SECTION_P_LIMIT_M = 5_000.0
ALIAS_FRACTION_CUT = 0.25
VEIL_ALPHA = 0.55
PLAN_X_MIN, PLAN_X_MAX = 672_000.0, 695_000.0
PLAN_Y_MIN, PLAN_Y_MAX = 1_335_000.0, 1_365_000.0

# Stress figure display.
STRESS_LIMIT_MPA = 20.0
DIFFERENCE_LIMIT_MPA = 2.0
MAP_STRESS_LIMIT_MPA = 5.0
MAP_DEPTH_M = 1_000.0
MAP_BOX_ASPECT = 0.68

RC_PARAMS = {
    "font.family": "DejaVu Sans",
    "font.size": 7.5,
    "axes.titlesize": 8.2,
    "axes.labelsize": 7,
    "figure.titlesize": 11.5,
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_dem() -> Path:
    if not DEM_PATH.is_file():
        raise SystemExit(
            f"missing terrain tile: {DEM_PATH}\n"
            f"Download it from {DEM_URL}\n"
            "See data/README.md."
        )
    actual = sha256_file(DEM_PATH)
    if actual != DEM_SHA256:
        raise SystemExit(f"terrain tile SHA-256 {actual} does not match the pinned hash")
    return DEM_PATH


WRIGHT_FIT_GSLC = "distributed_fit_geometry_4a7b2386c8b5_gslc.json"


def load_gslc_fit(name: str) -> tuple[dict, dict, np.ndarray]:
    payload = json.loads((RESULTS / name).read_text())
    if payload.get("nisar_source") != "gslc":
        raise ValueError(
            f"{name}: nisar_source is {payload.get('nisar_source')!r}; the "
            "summary figures draw GSLC observations, so the fit must use "
            "the same chain")
    geometry = geometry_from_json(payload["geometry"])
    if payload.get("geometry_sha256") != geometry_sha256(geometry):
        raise ValueError(f"{name} geometry hash mismatch")
    magnitudes = np.asarray(payload["dislocation_magnitudes_m"], dtype=float)
    if magnitudes.shape != geometry["kinds"].shape:
        raise ValueError("saved dislocations do not match the saved geometry")
    return payload, geometry, magnitudes


def load_mapped_gslc_fit() -> tuple[dict, dict, np.ndarray]:
    return load_gslc_fit("distributed_fit_mapped_gslc.json")


def load_wright_gslc_fit() -> tuple[dict, dict, np.ndarray]:
    # The file name pins the geometry-contract hash: a revised contract
    # changes the hash, so this figure cannot silently rebuild from a
    # stale fit of the old plane.
    payload, geometry, magnitudes = load_gslc_fit(WRIGHT_FIT_GSLC)
    if "wright2006_provenance" not in geometry:
        raise ValueError(f"{WRIGHT_FIT_GSLC} was not fitted on the Wright plane")
    return payload, geometry, magnitudes


# ---------------------------------------------------------------------------
# Shared map style: display grid, terrain, graticule, and source geometry.


def display_grid() -> tuple[np.ndarray, np.ndarray, object]:
    width = int(round((X_MAX - X_MIN) / DISPLAY_STEP_M))
    height = int(round((Y_MAX - Y_MIN) / DISPLAY_STEP_M))
    xs = X_MIN + DISPLAY_STEP_M * (np.arange(width) + 0.5)
    ys = Y_MAX - DISPLAY_STEP_M * (np.arange(height) + 0.5)
    transform = from_origin(X_MIN, Y_MAX, DISPLAY_STEP_M, DISPLAY_STEP_M)
    return xs, ys, transform


def load_terrain(destination_transform: object, shape: tuple[int, int]) -> np.ndarray:
    with rasterio.open(require_dem()) as dataset:
        source = dataset.read(1, masked=True).filled(np.nan)
        terrain = np.full(shape, np.nan, dtype=np.float64)
        reproject(
            source=source,
            destination=terrain,
            src_transform=dataset.transform,
            src_crs=dataset.crs,
            src_nodata=dataset.nodata,
            dst_transform=destination_transform,
            dst_crs=UTM_CRS,
            dst_nodata=np.nan,
            resampling=Resampling.average,
            init_dest_nodata=True,
        )
    if np.isfinite(terrain).mean() < 0.99:
        raise ValueError("terrain tile does not cover the display frame")
    return gaussian_filter(terrain, sigma=0.65)


def hillshade(terrain: np.ndarray) -> np.ndarray:
    shade = LightSource(azdeg=315, altdeg=42).hillshade(
        terrain,
        vert_exag=1.35,
        dx=DISPLAY_STEP_M,
        dy=DISPLAY_STEP_M,
        fraction=1.1,
    )
    return 0.72 + 0.28 * shade


def blended_map(
    field: np.ndarray,
    shade: np.ndarray,
    cmap: object,
    norm: object,
    *,
    alpha: float,
) -> np.ndarray:
    background = np.repeat(shade[:, :, None], 3, axis=2)
    colors = cmap(norm(field))[..., :3]
    colors = np.clip(colors * (0.78 + 0.30 * shade[:, :, None]), 0, 1)
    result = background.copy()
    valid = np.isfinite(field)
    result[valid] = (1 - alpha) * background[valid] + alpha * colors[valid]
    return result


def _require_straight_frame(geometry: dict) -> None:
    # The Wright-2006 contract's scalar strike/trace_centre are only the
    # chord of a segmented plane, and its grafted fault-axis metadata refers
    # to a different frame; rebuilding positions from them here would place
    # sources tens of kilometres off their saved patches.
    if "wright2006_provenance" in geometry:
        raise ValueError(
            "geometry uses the segmented Wright-2006 plane; this module's "
            "straight-frame reconstruction does not apply — use the patch "
            "arrays")


def geometry_lines(geometry: dict) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    _require_straight_frame(geometry)
    e_c, n_c = geometry["trace_centre"]
    strike = np.radians(geometry["strike"])
    along_e = np.sin(strike)
    along_n = np.cos(strike)
    dike_s = np.linspace(
        geometry["dike_s_edges_m"][0], geometry["dike_s_edges_m"][-1], 200
    )
    lines = {
        "dike": (e_c + dike_s * along_e, n_c + dike_s * along_n),
    }
    if "east" in geometry["enabled_fault_sides"]:
        s = np.asarray(geometry["fault_s_centres_m"], dtype=float)
        p = np.asarray(geometry["fault_column_offsets_m"]["east"], dtype=float)
        lines["fault"] = (
            e_c + s * along_e - p * np.cos(strike),
            n_c + s * along_n + p * np.sin(strike),
        )
    return lines


def draw_geometry(ax: object, lines: dict, *, legend: bool = False) -> None:
    halo = [
        path_effects.Stroke(linewidth=2.1, foreground="white"),
        path_effects.Normal(),
    ]
    dike = ax.plot(
        *lines["dike"], color="#b21f2d", lw=1.05, label="dike", zorder=6
    )[0]
    dike.set_path_effects(halo)
    if "fault" in lines:
        fault = ax.plot(
            *lines["fault"],
            color="#202020",
            lw=0.9,
            ls=(0, (4, 2)),
            label="east ridge",
            zorder=6,
        )[0]
        fault.set_path_effects(halo)
    if legend:
        legend_box = ax.legend(
            loc="lower right",
            fontsize=6.0,
            frameon=True,
            framealpha=0.9,
            borderpad=0.35,
            handlelength=2.4,
        )
        # above the haloed geometry lines, so they cannot cross its text
        legend_box.set_zorder(10)


def _utm_to_lon(easting_m: np.ndarray, northing_m: float) -> np.ndarray:
    values = np.ravel(np.asarray(easting_m, dtype=float))
    if values.size == 0:
        return np.asarray(easting_m, dtype=float)
    lon, _ = warp_transform(
        UTM_CRS, "EPSG:4326", values.tolist(), [northing_m] * values.size
    )
    return np.asarray(lon, dtype=float).reshape(np.shape(easting_m))


def _lon_to_utm(lon: np.ndarray, northing_m: float) -> np.ndarray:
    values = np.ravel(np.asarray(lon, dtype=float))
    if values.size == 0:
        return np.asarray(lon, dtype=float)
    _, lat = warp_transform(
        UTM_CRS, "EPSG:4326", [0.5 * (X_MIN + X_MAX)], [northing_m]
    )
    easting, _ = warp_transform(
        "EPSG:4326", UTM_CRS, values.tolist(), [lat[0]] * values.size
    )
    return np.asarray(easting, dtype=float).reshape(np.shape(lon))


def _utm_to_lat(northing_m: np.ndarray, easting_m: float) -> np.ndarray:
    values = np.ravel(np.asarray(northing_m, dtype=float))
    if values.size == 0:
        return np.asarray(northing_m, dtype=float)
    _, lat = warp_transform(
        UTM_CRS, "EPSG:4326", [easting_m] * values.size, values.tolist()
    )
    return np.asarray(lat, dtype=float).reshape(np.shape(northing_m))


def _lat_to_utm(lat: np.ndarray, easting_m: float) -> np.ndarray:
    values = np.ravel(np.asarray(lat, dtype=float))
    if values.size == 0:
        return np.asarray(lat, dtype=float)
    lon, _ = warp_transform(
        UTM_CRS, "EPSG:4326", [easting_m], [0.5 * (Y_MIN + Y_MAX)]
    )
    _, northing = warp_transform(
        "EPSG:4326", UTM_CRS, [lon[0]] * values.size, values.tolist()
    )
    return np.asarray(northing, dtype=float).reshape(np.shape(lat))


def add_lonlat_axes(
    ax: object,
    *,
    top: bool = False,
    right: bool = False,
    step_deg: float = 0.1,
    fontsize: float = 6.0,
) -> None:
    """Add geographic degree ticks on the outer edges of a UTM map panel.

    The tick mapping fixes the frame's top northing (for longitude) or right
    easting (for latitude); across this frame the resulting placement error
    is far below one tick spacing.
    """
    if top:
        y_ref = max(ax.get_ylim())
        secondary = ax.secondary_xaxis(
            "top",
            functions=(
                lambda x: _utm_to_lon(np.asarray(x), y_ref),
                lambda lon: _lon_to_utm(lon, y_ref),
            ),
        )
        secondary.xaxis.set_major_locator(MultipleLocator(step_deg))
        secondary.xaxis.set_major_formatter(
            FuncFormatter(lambda value, _: f"{value:.1f}°E")
        )
        secondary.tick_params(labelsize=fontsize, length=2.0, width=0.6, pad=1.2)
        secondary.spines["top"].set_linewidth(0.65)
    if right:
        x_ref = max(ax.get_xlim())
        secondary = ax.secondary_yaxis(
            "right",
            functions=(
                lambda y: _utm_to_lat(np.asarray(y), x_ref),
                lambda lat: _lat_to_utm(lat, x_ref),
            ),
        )
        secondary.yaxis.set_major_locator(MultipleLocator(step_deg))
        secondary.yaxis.set_major_formatter(
            FuncFormatter(lambda value, _: f"{value:.1f}°N")
        )
        secondary.tick_params(labelsize=fontsize, length=2.0, width=0.6, pad=1.2)
        secondary.spines["right"].set_linewidth(0.65)


def add_scale_and_north(ax: object) -> None:
    x0, y0 = X_MIN + 2_000, Y_MIN + 2_000
    ax.plot([x0, x0 + 10_000], [y0, y0], color="black", lw=2.4, zorder=8)
    ax.plot([x0, x0], [y0 - 300, y0 + 300], color="black", lw=1, zorder=8)
    ax.plot(
        [x0 + 10_000, x0 + 10_000],
        [y0 - 300, y0 + 300],
        color="black",
        lw=1,
        zorder=8,
    )
    ax.text(x0 + 5_000, y0 + 700, "10 km", ha="center", va="bottom", fontsize=6.5)
    ax.annotate(
        "N",
        xy=(X_MIN + 2_000, Y_MAX - 1_600),
        xytext=(X_MIN + 2_000, Y_MAX - 5_400),
        ha="center",
        va="center",
        fontsize=7,
        fontweight="bold",
        arrowprops={"arrowstyle": "-|>", "lw": 1.1, "color": "black"},
        zorder=8,
    )


def add_panel_label(ax: object, index: int, *, right: bool = False) -> None:
    ax.text(
        0.975 if right else 0.018,
        0.975,
        f"({string.ascii_lowercase[index]})",
        transform=ax.transAxes,
        ha="right" if right else "left",
        va="top",
        fontsize=7.2,
        fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.3},
        zorder=10,
    )


# ---------------------------------------------------------------------------
# Strike-frame geometry shared by the source, gradient, and stress figures.


def local_basis(geometry: dict) -> tuple[np.ndarray, np.ndarray]:
    _require_straight_frame(geometry)
    strike = np.radians(geometry["strike"])
    along_ne = np.array([np.cos(strike), np.sin(strike)])
    east_side_ne = np.array([np.sin(strike), -np.cos(strike)])
    return along_ne, east_side_ne


def patch_corners_local(geometry: dict) -> np.ndarray:
    """Patch corners in the (s, p, z) strike frame, kilometres.

    ``s`` is along the saved strike (positive SSE), ``p`` is perpendicular
    (positive ENE, the east side), and ``z`` is depth. The frame matches
    ``mh26.distributed._strike_vectors``.
    """
    along_ne, east_side_ne = local_basis(geometry)
    centre_e, centre_n = geometry["trace_centre"]
    origin_ne = np.array([centre_n, centre_e])
    polygons = []
    for patch in geometry["patches"]:
        north, east, depth, strike_deg, dip_deg, al1, al2, aw1, aw2 = patch
        strike = np.radians(strike_deg)
        dip = np.radians(dip_deg)
        patch_along = np.array([np.cos(strike), np.sin(strike), 0.0])
        patch_down_dip = np.array(
            [-np.sin(strike) * np.cos(dip), np.cos(strike) * np.cos(dip), np.sin(dip)]
        )
        centre = np.array([north, east, depth])
        corners_nez = np.array(
            [
                centre + al1 * patch_along + aw1 * patch_down_dip,
                centre + al2 * patch_along + aw1 * patch_down_dip,
                centre + al2 * patch_along + aw2 * patch_down_dip,
                centre + al1 * patch_along + aw2 * patch_down_dip,
            ]
        )
        relative_ne = corners_nez[:, :2] - origin_ne
        s = relative_ne @ along_ne
        p = relative_ne @ east_side_ne
        polygons.append(np.column_stack([s, p, corners_nez[:, 2]]) / 1_000)
    return np.asarray(polygons)


def peak_fault_section_s_m(geometry: dict, magnitudes: np.ndarray) -> float:
    """Along-strike position (m) of the east-fault patch with peak slip.

    The source and gradient figures derive their cross-section from this one
    value, so a refit cannot move one figure and not the other.
    """
    polygons = patch_corners_local(geometry)
    fault = geometry["kinds"] == "fault_east"
    centres_s_km = polygons[fault][:, :, 0].mean(axis=1)
    peak = int(np.argmax(np.asarray(magnitudes, dtype=float)[fault]))
    return float(centres_s_km[peak] * 1_000)


def signed_cross_rift_increment(
    phase: np.ndarray,
    coherence: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    strike_deg: float,
) -> np.ndarray:
    """Signed wrapped phase change across one native pixel, positive ENE."""
    valid = np.isfinite(phase) & (coherence > COHERENCE_MIN_GRADIENT)
    east = np.full(phase.shape, np.nan, dtype=np.float64)
    north = np.full(phase.shape, np.nan, dtype=np.float64)

    east_pairs = valid[:, :-1] & valid[:, 1:]
    east_diff = np.angle(np.exp(1j * (phase[:, 1:] - phase[:, :-1])))
    east[:, :-1][east_pairs] = east_diff[east_pairs]

    north_pairs = valid[:-1, :] & valid[1:, :]
    row_diff = np.angle(np.exp(1j * (phase[1:, :] - phase[:-1, :])))
    y_direction = np.sign(float(np.median(np.diff(y))))
    north[:-1, :][north_pairs] = y_direction * row_diff[north_pairs]

    dx = float(np.median(np.diff(x)))
    dy = float(abs(np.median(np.diff(y))))
    if not (np.isclose(dx, 30.0, atol=0.5) and np.isclose(dy, 30.0, atol=0.5)):
        raise ValueError(
            f"signed phase display expects 30-m pixels; got {dx:.2f} by {dy:.2f} m"
        )

    strike = np.radians(strike_deg)
    perp_east = -np.cos(strike)
    perp_north = np.sin(strike)
    return perp_east * east + perp_north * north


def aliasing_flag(
    los: np.ndarray,
    y: np.ndarray,
    strike_deg: float,
    alias_limit_m: float,
) -> np.ndarray:
    """True where the unwrapped cross-rift LOS rate exceeds ``alias_limit_m``
    per pixel, so the wrapped difference can exceed π and alias.

    Uses the same forward-difference stencil and ENE projection as the
    signed phase increment, so the flag aligns with the signed map. Pixels
    without a computable rate stay unflagged.
    """
    valid = np.isfinite(los)
    east = np.full(los.shape, np.nan, dtype=np.float64)
    north = np.full(los.shape, np.nan, dtype=np.float64)
    east_pairs = valid[:, :-1] & valid[:, 1:]
    east[:, :-1][east_pairs] = (los[:, 1:] - los[:, :-1])[east_pairs]
    north_pairs = valid[:-1, :] & valid[1:, :]
    y_direction = np.sign(float(np.median(np.diff(y))))
    north[:-1, :][north_pairs] = (
        y_direction * (los[1:, :] - los[:-1, :])
    )[north_pairs]
    strike = np.radians(strike_deg)
    rate = -np.cos(strike) * east + np.sin(strike) * north
    return np.abs(rate) > alias_limit_m


def source_transform(
    field: np.ndarray, xs: np.ndarray, ys: np.ndarray
) -> tuple[np.ndarray, object]:
    values = np.asarray(field, dtype=np.float64)
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    if x[0] > x[-1]:
        values = values[:, ::-1]
        x = x[::-1]
    if y[0] < y[-1]:
        values = values[::-1]
        y = y[::-1]
    dx = float(np.median(np.diff(x)))
    dy = float(abs(np.median(np.diff(y))))
    return values, from_origin(x[0] - dx / 2, y[0] + dy / 2, dx, dy)


def reproject_field(
    field: np.ndarray,
    xs: np.ndarray,
    ys: np.ndarray,
    source_crs: str,
    destination_transform: object,
    destination_shape: tuple[int, int],
    *,
    resampling: Resampling,
) -> np.ndarray:
    source, transform = source_transform(field, xs, ys)
    destination = np.full(destination_shape, np.nan, dtype=np.float64)
    reproject(
        source=source,
        destination=destination,
        src_transform=transform,
        src_crs=source_crs,
        src_nodata=np.nan,
        dst_transform=destination_transform,
        dst_crs=UTM_CRS,
        dst_nodata=np.nan,
        resampling=resampling,
        init_dest_nodata=True,
    )
    return destination


# ---------------------------------------------------------------------------
# Figure 1: observed, model, and residual LOS with profiles (four geometries).


LOOKS = tuple(
    [
        {
            "key": tag,
            "kind": "nisar",
            "tag": tag,
            "label": NISAR_LABEL[tag],
            "payload_name": NISAR_LABEL[tag].replace("NISAR", "NISAR/GSLC"),
        }
        for tag in INVERSION_PAIRS
    ]
    + [
        {
            "key": f"{frame}:{pair}",
            "kind": "s1",
            "frame": frame,
            "pair": pair,
            "label": S1_GEOMETRY[frame],
            "payload_name": S1_GEOMETRY[frame],
        }
        for frame, pair, role in S1_TRACK_PAIRS
        if role == "co-event"
    ]
)


def load_nisar_look(
    look: dict,
    factor: float,
    transform: object,
    shape: tuple[int, int],
    points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    source = np.load(DATA_GSLC / f"{look['tag']}_unw.npz")
    los = np.where(
        source["coherence"] > UNWRAP_COHERENCE_MIN,
        source["los"].astype(np.float64) / factor,
        np.nan,
    )
    observation = reproject_field(
        los,
        source["x"],
        source["y"],
        UTM_CRS,
        transform,
        shape,
        resampling=Resampling.average,
    )
    geometry_cube = load_pair(look["tag"])
    units = los_unit_neu(geometry_cube, points)
    return observation, units.reshape(*shape, 3)


def load_s1_look(
    look: dict,
    model: np.ndarray,
    factor: float,
    transform: object,
    shape: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    frame = load_frame(look["frame"])
    referenced = reference(frame, load_pair_los(frame, look["pair"]))
    statistics = coevent_pixel_stats(frame, referenced, model)
    # The fit refuses frames below this gate; re-check it here so an input
    # change cannot silently flip a sign the fit would reject.
    if abs(statistics["correlation"]) < MIN_COEVENT_ABS_R:
        raise ValueError(
            f"{look['label']}: co-event correlation "
            f"{statistics['correlation']:.3f} is below the sign gate "
            f"{MIN_COEVENT_ABS_R}"
        )
    aligned = statistics["applied_observation_sign"] * referenced / factor
    observation = reproject_field(
        aligned,
        frame["lons"],
        frame["lats"],
        "EPSG:4326",
        transform,
        shape,
        resampling=Resampling.average,
    )
    components = [
        reproject_field(
            np.where(frame["vec_ok"], frame[name], np.nan),
            frame["lons"],
            frame["lats"],
            "EPSG:4326",
            transform,
            shape,
            resampling=Resampling.bilinear,
        )
        for name in ("n", "e", "u")
    ]
    units = np.stack(components, axis=-1)
    norm = np.linalg.norm(units, axis=-1)
    good = np.isfinite(norm) & (norm > 0.9)
    units[good] /= norm[good, None]
    units[~good] = np.nan
    return observation, units


def model_surface_neu(
    geometry: dict, magnitudes: np.ndarray, points: np.ndarray
) -> np.ndarray:
    disl = geometry["disl"] * magnitudes[:, None]
    receivers = np.zeros((len(points), 3), dtype=np.float64)
    receivers[:, :2] = points
    result = okada_ext.okada(
        geometry["patches"],
        disl,
        receivers,
        LAME_LAMBDA,
        SHEAR_MODULUS,
        nthreads=0,
        rotate_sdn=0,
        stack_sources=1,
    )
    return np.column_stack([result[:, 0], result[:, 1], -result[:, 2]])


def profile_values(
    field: np.ndarray,
    xs: np.ndarray,
    ys: np.ndarray,
    frame: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """Median across a 1-km swath of seven rift-perpendicular lines."""
    interpolator = RegularGridInterpolator(
        (ys[::-1], xs), field[::-1], bounds_error=False, fill_value=np.nan
    )
    distance = np.linspace(-18_000.0, 18_000.0, 361)
    strike = np.radians(frame["strike_deg"])
    e_c, n_c = frame["origin_e"], frame["origin_n"]
    samples = []
    for along_offset in np.linspace(-500.0, 500.0, 7):
        east = e_c + along_offset * np.sin(strike) - distance * np.cos(strike)
        north = n_c + along_offset * np.cos(strike) + distance * np.sin(strike)
        samples.append(interpolator(np.column_stack([north, east])))
    stacked = np.stack(samples)
    output = np.full(stacked.shape[1], np.nan)
    valid = np.isfinite(stacked).any(axis=0)
    output[valid] = np.nanmedian(stacked[:, valid], axis=0)
    return distance / 1e3, output


def profile_line(frame: dict) -> tuple[np.ndarray, np.ndarray]:
    distance = np.array([-18_000.0, 18_000.0])
    strike = np.radians(frame["strike_deg"])
    e_c, n_c = frame["origin_e"], frame["origin_n"]
    return (
        e_c - distance * np.cos(strike),
        n_c + distance * np.sin(strike),
    )


def draw_profile_line(ax: object, frame: dict, *, labels: bool = False) -> None:
    east, north = profile_line(frame)
    line = ax.plot(
        east, north, color="#303030", lw=0.65, ls=(0, (3, 2)), zorder=5
    )[0]
    line.set_path_effects(
        [path_effects.Stroke(linewidth=1.4, foreground="white"), path_effects.Normal()]
    )
    if labels:
        ax.text(east[0], north[0], "A", fontsize=6.5, fontweight="bold",
                ha="right", va="top", zorder=8)
        ax.text(east[1], north[1], "A′", fontsize=6.5, fontweight="bold",
                ha="left", va="bottom", zorder=8)


def straight_profile_frame(geometry: dict) -> dict:
    _require_straight_frame(geometry)
    centre_e, centre_n = geometry["trace_centre"]
    return {"origin_e": float(centre_e), "origin_n": float(centre_n),
            "strike_deg": float(geometry["strike"])}


def patch_top_trace(geometry: dict, kind: str) -> tuple[np.ndarray, np.ndarray]:
    """Surface trace of one source sheet from its saved shallowest patch
    edges, as NaN-separated (easting, northing) segments.

    Segmented geometries have no usable scalar chord frame (see
    ``_require_straight_frame``), so the trace must come from the patches.
    """
    kinds = np.asarray(geometry["kinds"])
    patches = np.asarray(geometry["patches"], dtype=float)[kinds == kind]
    if patches.size == 0:
        raise ValueError(f"geometry has no {kind!r} patches")
    strike = np.radians(patches[:, 3])
    dip = np.radians(patches[:, 4])
    # Up-dip edge midpoints, in the patch frame of patch_corners_local.
    edge_n = patches[:, 0] - patches[:, 7] * np.sin(strike) * np.cos(dip)
    edge_e = patches[:, 1] + patches[:, 7] * np.cos(strike) * np.cos(dip)
    edge_z = patches[:, 2] + patches[:, 7] * np.sin(dip)
    top = edge_z < edge_z.min() + 1.0
    along_n, along_e = np.cos(strike), np.sin(strike)
    mean_n, mean_e = along_n[top].mean(), along_e[top].mean()
    order = np.argsort(edge_n[top] * mean_n + edge_e[top] * mean_e)
    east, north = [], []
    for index in np.flatnonzero(top)[order]:
        for along in (patches[index, 5], patches[index, 6]):
            east.append(edge_e[index] + along * along_e[index])
            north.append(edge_n[index] + along * along_n[index])
        east.append(np.nan)
        north.append(np.nan)
    return np.asarray(east), np.asarray(north)


def patch_profile_frame(geometry: dict, magnitudes: np.ndarray) -> dict:
    """Profile frame from the saved patches: origin at the opening-weighted
    dike centre, azimuth the opening-weighted patch strike.

    The 180° flip keeps the along axis within 90° of the saved chord, so
    positive profile distance stays ENE (the map figures' reading order)
    whichever sense the patch strikes are stored in.
    """
    kinds = np.asarray(geometry["kinds"])
    dike = kinds == "dike"
    patches = np.asarray(geometry["patches"], dtype=float)[dike]
    if not np.allclose(patches[:, 4], 90.0):
        raise ValueError("the profile origin assumes a vertical dike")
    weights = (np.clip(np.asarray(magnitudes, dtype=float)[dike], 0.0, None)
               * np.asarray(geometry["areas"], dtype=float)[dike])
    total = weights.sum()
    if total <= 0.0:
        raise ValueError("no positive dike opening to place the profile on")
    weights = weights / total
    strikes = np.radians(patches[:, 3])
    mean = np.arctan2(weights @ np.sin(strikes), weights @ np.cos(strikes))
    if np.cos(mean - np.radians(float(geometry["strike"]))) < 0.0:
        mean += np.pi
    return {"origin_e": float(weights @ patches[:, 1]),
            "origin_n": float(weights @ patches[:, 0]),
            "strike_deg": float(np.degrees(mean) % 360.0)}


def patch_fault_offset_km(geometry: dict, frame: dict) -> float | None:
    """Signed profile-axis distance (km, ENE positive) where the east-fault
    trace passes closest to the profile origin."""
    if "east" not in geometry["enabled_fault_sides"]:
        return None
    east, north = patch_top_trace(geometry, "fault_east")
    finite = np.isfinite(east)
    relative_e = east[finite] - frame["origin_e"]
    relative_n = north[finite] - frame["origin_n"]
    strike = np.radians(frame["strike_deg"])
    along = relative_e * np.sin(strike) + relative_n * np.cos(strike)
    perpendicular = -relative_e * np.cos(strike) + relative_n * np.sin(strike)
    return float(perpendicular[np.argmin(np.abs(along))] / 1e3)


def draw_look_direction(ax: object, unit: np.ndarray, kind: str) -> None:
    horizontal = np.array([unit[1], unit[0]], dtype=float)
    horizontal /= np.linalg.norm(horizontal)
    # Flight direction from the horizontal ground-to-satellite vector and the
    # mission's fixed look side: Sentinel-1 is right-looking, NISAR is
    # left-looking, so the heading is the look azimuth rotated 90 degrees.
    if kind == "s1":
        flight = np.array([horizontal[1], -horizontal[0]])
    else:
        flight = np.array([-horizontal[1], horizontal[0]])
    halo = [
        path_effects.Stroke(linewidth=2.6, foreground="white"),
        path_effects.Normal(),
    ]
    centre = np.array([0.86, 0.85])
    flight_delta = 0.19 * flight
    flight_arrow = ax.annotate(
        "",
        xy=centre + flight_delta / 2,
        xytext=centre - flight_delta / 2,
        xycoords="axes fraction",
        arrowprops={
            "arrowstyle": "-|>",
            "lw": 0.9,
            "color": "#111111",
            "mutation_scale": 6.5,
        },
        zorder=9,
    )
    flight_arrow.arrow_patch.set_path_effects(halo)
    look_arrow = ax.annotate(
        "",
        xy=centre + 0.115 * horizontal,
        xytext=centre + 0.02 * horizontal,
        xycoords="axes fraction",
        arrowprops={
            "arrowstyle": "-|>",
            "lw": 1.25,
            "color": "#111111",
            "mutation_scale": 8.5,
        },
        zorder=9,
    )
    look_arrow.arrow_patch.set_path_effects(halo)
    ax.text(
        centre[0],
        centre[1] + 0.105,
        "+LOS",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.2,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.8},
        zorder=9,
    )


def draw_los_map_axes(
    ax: object,
    row: int,
    column: int,
    terrain: np.ndarray,
    xs: np.ndarray,
    ys: np.ndarray,
) -> None:
    levels = np.arange(
        np.ceil(np.nanmin(terrain) / 100) * 100,
        np.floor(np.nanmax(terrain) / 100) * 100 + 1,
        200,
    )
    ax.contour(
        xs, ys, terrain, levels=levels,
        colors="#303030", linewidths=0.28, alpha=0.24, zorder=4,
    )
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_aspect("equal")
    ax.set_xticks([665_000, 675_000, 685_000, 695_000])
    ax.set_yticks([1_335_000, 1_345_000, 1_355_000, 1_365_000])
    if row == len(LOOKS) - 1:
        ax.set_xticklabels(["665", "675", "685", "695"])
        if column == 1:
            ax.set_xlabel("UTM 37N easting (km)", fontsize=7)
    else:
        ax.set_xticklabels([])
    if column == 0:
        ax.set_yticklabels(["1335", "1345", "1355", "1365"])
        ax.set_ylabel("UTM 37N northing (km)", fontsize=7)
    else:
        ax.set_yticklabels([])
    ax.tick_params(length=2.5, width=0.6, labelsize=6.5, pad=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.65)
        spine.set_color("#303030")


def los_model_residual_figure(
    payload: dict,
    geometry: dict,
    magnitudes: np.ndarray,
    *,
    lines: dict,
    profile_frame: dict,
    fault_offset_km: float | None,
    stem: str,
    suptitle: str,
    subtitle: str,
) -> list[Path]:
    model_okada = np.asarray(
        json.loads((RESULTS / "okada_fit.json").read_text())["values"][:7],
        dtype=float,
    )
    factors = {row["name"]: row["decay_factor"] for row in payload["datasets"]}

    xs, ys, transform = display_grid()
    shape = (len(ys), len(xs))
    X, Y = np.meshgrid(xs, ys)
    points = np.column_stack([Y.ravel(), X.ravel()])
    terrain = load_terrain(transform, shape)
    shade = hillshade(terrain)

    observations, look_units, unit_medians = {}, {}, {}
    for look in LOOKS:
        if look["kind"] == "nisar":
            observation, units = load_nisar_look(
                look, factors[look["payload_name"]], transform, shape, points
            )
        else:
            observation, units = load_s1_look(
                look, model_okada, factors[look["payload_name"]], transform, shape
            )
        observations[look["key"]] = observation
        look_units[look["key"]] = units
        unit_medians[look["key"]] = np.nanmedian(units.reshape(-1, 3), axis=0)
        print(f"loaded {look['label']}")

    total_neu = model_surface_neu(geometry, magnitudes, points).reshape(*shape, 3)
    planes = {row["name"]: row["model_farfield_plane_m"]
              for row in payload["datasets"]}
    models = {}
    for look in LOOKS:
        prediction = np.einsum("ijk,ijk->ij", total_neu, look_units[look["key"]])
        prediction += payload["dataset_offsets_m"][look["payload_name"]]
        coef = planes[look["payload_name"]]
        # The observed fields are far-field plane referenced; show the model
        # through the estimator's matching projection so the residual panels
        # equal the fitted residuals.
        prediction -= (coef[0] + coef[1] * (X - CENTER_E)
                       + coef[2] * (Y - CENTER_N))
        models[look["key"]] = prediction

    plt.rcParams.update(RC_PARAMS)
    cyclic_cmap = cmc.romaO
    residual_cmap = cmc.vik
    cycle_norm = Normalize(0, DISPLAY_CYCLE_M)
    residual_norm = TwoSlopeNorm(
        vmin=-RESIDUAL_LIMIT_M, vcenter=0, vmax=RESIDUAL_LIMIT_M
    )

    fig = plt.figure(figsize=(7.5, 8.0))
    # The fourth, empty column keeps the profile axis text and the residual
    # column's outer °N labels apart.
    grid = fig.add_gridspec(
        len(LOOKS),
        5,
        width_ratios=[1, 1, 1, 0.40, 1.18],
        left=0.075,
        right=0.945,
        bottom=0.16,
        top=0.91,
        hspace=0.10,
        wspace=0.10,
    )
    column_titles = ("Observed", "Distributed model", "Residual", "Profile AA′")
    profile_limit_cm = 85
    panel_index = 0
    first_row_map_axes = []
    profile_axes = []

    for row, look in enumerate(LOOKS):
        observation = observations[look["key"]]
        model_field = models[look["key"]]
        residual = observation - model_field
        map_fields = (
            np.mod(observation, DISPLAY_CYCLE_M),
            np.mod(model_field, DISPLAY_CYCLE_M),
            residual,
        )
        for column, field in enumerate(map_fields):
            ax = fig.add_subplot(grid[row, column])
            if row == 0:
                first_row_map_axes.append(ax)
            if column < 2:
                rgb = blended_map(field, shade, cyclic_cmap, cycle_norm, alpha=0.84)
            else:
                rgb = blended_map(field, shade, residual_cmap, residual_norm, alpha=0.79)
            ax.imshow(
                rgb,
                extent=[X_MIN, X_MAX, Y_MIN, Y_MAX],
                origin="upper",
                interpolation="nearest",
                zorder=1,
            )
            draw_los_map_axes(ax, row, column, terrain, xs, ys)
            draw_profile_line(ax, profile_frame, labels=(row == 0 and column == 0))
            draw_geometry(ax, lines, legend=(row == 0 and column == 1))
            if row == 0 and column == 1:
                for text in ax.get_legend().get_texts():
                    text.set_fontsize(6.0)
            if column == 0:
                draw_look_direction(ax, unit_medians[look["key"]], look["kind"])
            if row == 0 or column == 2:
                add_lonlat_axes(ax, top=(row == 0), right=(column == 2))
            if row == 0:
                ax.set_title(column_titles[column], fontweight="bold", pad=13)
            ax.text(
                0.018,
                0.975,
                f"({string.ascii_lowercase[panel_index]})",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=7.3,
                fontweight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                      "pad": 1.4},
                zorder=9,
            )
            if column == 0:
                ax.text(
                    0.018,
                    0.91,
                    look["label"],
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=6.6,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                          "pad": 1.3},
                    zorder=9,
                )
            if row == 0 and column == 0:
                add_scale_and_north(ax)
            if column == 2:
                rms_cm = 100 * payload["rms_m"]["distributed_with_faults"][
                    look["payload_name"]
                ]
                ax.text(
                    0.975,
                    0.025,
                    f"fit RMS {rms_cm:.1f} cm",
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=6.4,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84,
                          "pad": 1.4},
                    zorder=9,
                )
            panel_index += 1

        profile_ax = fig.add_subplot(grid[row, 4])
        profile_axes.append(profile_ax)
        distance, observation_profile = profile_values(
            observation, xs, ys, profile_frame)
        _, model_profile = profile_values(model_field, xs, ys, profile_frame)
        _, terrain_profile = profile_values(terrain, xs, ys, profile_frame)
        topography_ax = profile_ax.twinx()
        topography_ax.fill_between(
            distance, 0, terrain_profile,
            color="#d7d7d7", alpha=0.62, linewidth=0, zorder=0,
        )
        topography_ax.set_ylim(0, max(1_200, np.nanmax(terrain_profile) * 1.08))
        topography_ax.set_yticks([500, 1_000])
        topography_ax.tick_params(length=2, width=0.5, labelsize=6.0, pad=1)
        topography_ax.set_ylabel("elev. (m)", fontsize=6.2, color="#666666", labelpad=1)
        profile_ax.plot(distance, 100 * observation_profile, color="#202020",
                        lw=1.05, label="observed", zorder=4)
        profile_ax.plot(distance, 100 * model_profile, color="#c43c39",
                        lw=1.35, label="model", zorder=5)
        profile_ax.axvline(0, color="#b21f2d", lw=0.8, zorder=3)
        if fault_offset_km is not None:
            profile_ax.axvline(fault_offset_km, color="#202020", lw=0.75,
                               ls=(0, (4, 2)), zorder=3)
        profile_ax.axhline(0, color="#777777", lw=0.45, zorder=2)
        profile_ax.set_xlim(-18, 18)
        profile_ax.set_ylim(-profile_limit_cm, profile_limit_cm)
        profile_ax.set_yticks([-80, -40, 0, 40, 80])
        profile_ax.tick_params(length=2.5, width=0.6, labelsize=6.5, pad=1.5)
        profile_ax.set_ylabel("LOS (cm)", fontsize=6.4, labelpad=1)
        if row == len(LOOKS) - 1:
            profile_ax.set_xlabel("distance across dike, WSW → ENE (km)", fontsize=7)
        else:
            profile_ax.set_xticklabels([])
        if row == 0:
            profile_ax.set_title(column_titles[3], fontweight="bold", pad=13)
            profile_ax.legend(
                loc="upper left", ncols=2, fontsize=6.1, frameon=False,
                handlelength=2.2, columnspacing=1.1,
            )
        profile_ax.text(
            0.975, 0.975, f"({string.ascii_lowercase[panel_index]})",
            transform=profile_ax.transAxes, ha="right", va="top",
            fontsize=7.3, fontweight="bold", zorder=9,
        )
        profile_ax.grid(color="#b0b0b0", alpha=0.25, lw=0.45)
        for spine in profile_ax.spines.values():
            spine.set_linewidth(0.65)
        for spine in topography_ax.spines.values():
            spine.set_linewidth(0.65)
        profile_ax.set_box_aspect((Y_MAX - Y_MIN) / (X_MAX - X_MIN) / 1.18)
        panel_index += 1

    fig.canvas.draw()
    reference_positions = [ax.get_position() for ax in first_row_map_axes]
    colorbar_y = 0.075
    colorbar_height = 0.016
    cycle_axis = fig.add_axes(
        [
            reference_positions[0].x0,
            colorbar_y,
            reference_positions[1].x1 - reference_positions[0].x0,
            colorbar_height,
        ]
    )
    cycle_bar = fig.colorbar(
        ScalarMappable(norm=cycle_norm, cmap=cyclic_cmap),
        cax=cycle_axis,
        orientation="horizontal",
    )
    cycle_bar.set_ticks([0, DISPLAY_CYCLE_M / 2, DISPLAY_CYCLE_M])
    cycle_bar.set_ticklabels(["0", "5", "10"])
    cycle_bar.set_label(
        "displacement colour cycle\n(one cycle = 10 cm toward the satellite)",
        fontsize=7,
    )
    cycle_bar.ax.tick_params(labelsize=6.5, length=2)

    residual_axis = fig.add_axes(
        [
            reference_positions[2].x0,
            colorbar_y,
            reference_positions[2].width,
            colorbar_height,
        ]
    )
    residual_bar = fig.colorbar(
        ScalarMappable(norm=residual_norm, cmap=residual_cmap),
        cax=residual_axis,
        orientation="horizontal",
    )
    residual_bar.set_ticks([-RESIDUAL_LIMIT_M, 0, RESIDUAL_LIMIT_M])
    residual_bar.set_ticklabels(["−10", "0", "+10"])
    residual_bar.set_label("observed − model (cm; clipped)", fontsize=7)
    residual_bar.ax.tick_params(labelsize=6.5, length=2)

    profile_position = profile_axes[0].get_position()
    fig.text(
        profile_position.x0,
        0.083,
        "Profiles: + toward satellite\nmedian across a 1-km swath\n"
        "200-m display grid · gray = no data\nRMS at ≈1.3-km inversion blocks",
        ha="left",
        va="center",
        fontsize=6.0,
        color="#4a4a4a",
    )

    fig.suptitle(
        suptitle,
        x=0.5,
        y=0.992,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.944,
        subtitle,
        ha="center",
        va="bottom",
        fontsize=7.2,
        color="#3a3a3a",
    )

    png = RESULTS / f"{stem}.png"
    pdf = RESULTS / f"{stem}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def build_los_figure() -> list[Path]:
    payload, geometry, magnitudes = load_mapped_gslc_fit()
    fault_offset_km = None
    if "east" in geometry["enabled_fault_sides"]:
        offsets = np.asarray(geometry["fault_column_offsets_m"]["east"])
        s_centres = np.asarray(geometry["fault_s_centres_m"])
        fault_offset_km = float(offsets[np.argmin(np.abs(s_centres))]) / 1e3
    return los_model_residual_figure(
        payload,
        geometry,
        magnitudes,
        lines=geometry_lines(geometry),
        profile_frame=straight_profile_frame(geometry),
        fault_offset_km=fault_offset_km,
        stem="summary_los_model_residual",
        suptitle="Manda Hararo 2026 — four-geometry distributed-source fit",
        subtitle="reprocessed 30-m (GSLC) inputs · 18 July reference · "
                 "distributed opening + trace-supported east fault",
    )


def build_wright_los_figure() -> list[Path]:
    payload, geometry, magnitudes = load_wright_gslc_fit()
    profile_frame = patch_profile_frame(geometry, magnitudes)
    lines = {"dike": patch_top_trace(geometry, "dike")}
    if "east" in geometry["enabled_fault_sides"]:
        lines["fault"] = patch_top_trace(geometry, "fault_east")
    return los_model_residual_figure(
        payload,
        geometry,
        magnitudes,
        lines=lines,
        profile_frame=profile_frame,
        fault_offset_km=patch_fault_offset_km(geometry, profile_frame),
        stem="summary_los_model_residual_wright2006",
        suptitle="Manda Hararo 2026 — fit on the 2005 dike plane",
        subtitle="reprocessed 30-m (GSLC) inputs · 18 July reference · "
                 "SE-extended Wright et al. (2006) plane + "
                 "trace-supported east fault",
    )


def build_wright_uperp_figure() -> list[Path]:
    payload, geometry, magnitudes = load_wright_gslc_fit()
    provenance = geometry["wright2006_provenance"]
    nz, ns = geometry["dike_shape"]
    opening = np.asarray(magnitudes, dtype=float)[:nz * ns].reshape(nz, ns)

    s_edges = np.asarray(geometry["dike_s_edges_m"], dtype=float) / 1e3
    z_edges = np.asarray(geometry["dike_z_edges_m"], dtype=float) / 1e3
    s_centres = 0.5 * (s_edges[:-1] + s_edges[1:])
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    # Boundary patch values repeat at the edges so the contours cover the
    # saved model domain without extending beyond it (as in
    # wright2006.uperp_panel).
    s_plot = np.concatenate(([s_edges[0]], s_centres, [s_edges[-1]]))
    z_plot = np.concatenate(([z_edges[0]], z_centres, [z_edges[-1]]))
    opening_plot = np.pad(opening, 1, mode="edge")

    plt.rcParams.update(RC_PARAMS)
    dike_norm = Normalize(0, DIKE_LIMIT_M)
    extend = "max" if float(opening.max()) > DIKE_LIMIT_M else "neither"
    fig = plt.figure(figsize=(7.5, 3.0))
    ax = fig.add_axes([0.075, 0.175, 0.815, 0.545])
    filled = ax.contourf(
        s_plot,
        z_plot,
        opening_plot,
        levels=np.linspace(0.0, DIKE_LIMIT_M, 27),
        cmap=cmc.lajolla_r,
        extend=extend,
        antialiased=True,
    )
    filled.set_rasterized(True)
    line_levels = [level for level in (0.5, 1.0, 1.5, 2.0)
                   if level < float(opening.max())]
    if line_levels:
        ax.contour(s_plot, z_plot, opening_plot, levels=line_levels,
                   colors="0.15", linewidths=0.5)

    n_ext = provenance["extension_se"]["n_patches"]
    if n_ext:
        original_se = s_edges[ns - n_ext]
        ax.axvspan(original_se, s_edges[-1], color="0.65", alpha=0.12,
                   zorder=2)
        ax.axvline(original_se, color="0.25", lw=0.85, ls=(0, (4, 3)),
                   zorder=4)
        ax.text(
            original_se + 0.35,
            0.5,
            "2005 SE limit",
            transform=ax.get_xaxis_transform(),
            rotation=90,
            ha="left",
            va="center",
            fontsize=6.4,
            color="0.20",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 1.2},
            zorder=5,
        )
        ax.text(
            0.5 * (original_se + s_edges[-1]),
            0.04,
            f"{s_edges[-1] - original_se:.0f} km SE\nextension",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=6.2,
            color="0.20",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78,
                  "pad": 1.2},
            zorder=5,
        )

    ax.set_xlim(s_edges[0], s_edges[-1])
    ax.set_ylim(z_edges[0], z_edges[-1])
    ax.invert_yaxis()
    ax.set_xlabel("distance along the 2005 dike plane, NW → SE (km)")
    ax.set_ylabel("depth (km)")
    ax.tick_params(length=2.5, width=0.6, labelsize=6.8, pad=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.65)
    ax.set_title("July 2026", loc="left", fontsize=9.0, fontweight="bold",
                 pad=5)
    ax.set_title(
        f"0–{z_edges[-1]:g} km model depth · contours every 0.5 m",
        loc="right",
        fontsize=6.8,
        color="0.30",
        pad=6,
    )
    ax.text(
        0.985,
        0.94,
        f"max {payload['dike']['max_opening_m']:.2f} m · "
        f"{payload['dike']['dike_volume_km3']:.3f} km³",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.4,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
              "pad": 1.4},
        zorder=6,
    )

    colorbar_axis = fig.add_axes([0.905, 0.175, 0.016, 0.545])
    colorbar = fig.colorbar(
        ScalarMappable(norm=dike_norm, cmap=cmc.lajolla_r),
        cax=colorbar_axis,
        extend=extend,
    )
    colorbar.set_ticks([0, 1, 2, DIKE_LIMIT_M])
    colorbar.set_label("U$_\\perp$ (m)", fontsize=7)
    colorbar.ax.tick_params(labelsize=6.5, length=2)

    fig.suptitle(
        "Manda Hararo 2026 — dike opening on the Wright et al. (2006) "
        "2005 plane",
        x=0.5,
        y=0.985,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.86,
        "reprocessed 30-m (GSLC) inputs · 18 July reference · "
        "same fit as the map-view figure · s = 0 at the plane's NW end",
        ha="center",
        va="bottom",
        fontsize=7.2,
        color="#3a3a3a",
    )

    png = RESULTS / "summary_uperp_wright2006.png"
    pdf = RESULTS / "summary_uperp_wright2006.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


# ---------------------------------------------------------------------------
# Figure 2: free East/North/Up decomposition with observational support.


def validate_recomputed_fields(fields: dict[str, np.ndarray], summary: dict) -> None:
    solved = np.isfinite(fields["e"])
    expected = summary["free_3d"]
    checks = {
        "solved cells": (int(solved.sum()), expected["n_cells"]),
        "three-geometry cells": (
            int((fields["n_geometries"] == 3).sum()),
            expected["n_cells_3geo"],
        ),
        "four-geometry cells": (
            int((fields["n_geometries"] == 4).sum()),
            expected["n_cells_4geo"],
        ),
    }
    for label, (actual, stored) in checks.items():
        if actual != stored:
            raise ValueError(f"{label} changed: recomputed {actual}, stored {stored}")
    for component in ("east", "north", "up"):
        key = component[0] if component != "up" else "u"
        actual = float(100 * np.nanmedian(fields[f"sigma_{key}"]))
        stored = float(expected["median_sigma_cm"][component])
        if not np.isclose(actual, stored, atol=1e-9):
            raise ValueError(
                f"median sigma {component} changed: recomputed {actual}, "
                f"stored {stored}"
            )


def enu_display_mask(centres_e: np.ndarray, centres_n: np.ndarray) -> np.ndarray:
    east, north = np.meshgrid(centres_e, centres_n)
    return (east >= X_MIN) & (east <= X_MAX) & (north >= Y_MIN) & (north <= Y_MAX)


def cell_edges(centres: np.ndarray) -> np.ndarray:
    step = float(np.median(np.diff(centres)))
    return np.r_[centres - step / 2, centres[-1] + step / 2]


def draw_enu_background(
    ax: object, terrain: np.ndarray, terrain_x: np.ndarray, terrain_y: np.ndarray
) -> None:
    shade = hillshade(terrain)
    ax.imshow(
        shade,
        cmap="gray",
        vmin=0.70,
        vmax=1.0,
        extent=[X_MIN, X_MAX, Y_MIN, Y_MAX],
        origin="upper",
        interpolation="nearest",
        zorder=0,
    )
    levels = np.arange(
        np.ceil(np.nanmin(terrain) / 100) * 100,
        np.floor(np.nanmax(terrain) / 100) * 100 + 1,
        200,
    )
    ax.contour(
        terrain_x, terrain_y, terrain, levels=levels,
        colors="#303030", linewidths=0.28, alpha=0.24, zorder=2,
    )


def decomposition_inputs() -> tuple[dict, dict, np.ndarray, dict, list[dict]]:
    """Recompute the free E–N–U solve from the GSLC inputs and check it."""
    summary = json.loads((RESULTS / "decomposition_joint_gslc.json").read_text())
    _, geometry, _ = load_mapped_gslc_fit()
    model = np.asarray(
        json.loads((RESULTS / "okada_fit.json").read_text())["values"][:7],
        dtype=float,
    )
    decay = json.loads((RESULTS / "s1_tracks.json").read_text())["post_event_decay"]
    cells = load_cells(model, decay, nisar_source="gslc")
    fields = solve_cells_3d(cells)
    validate_recomputed_fields(fields, summary)
    return summary, geometry, model, fields, cells


def draw_frame_axes(ax: object, strike_deg: float) -> None:
    """Arrows for +U⊥ and +U∥ of the dike frame, lower right of a map."""
    along, perp = _strike_vectors(strike_deg)
    x0, y0 = X_MAX - 7_500.0, Y_MIN + 7_000.0
    length = 4_500.0
    for (v_n, v_e), label in ((perp, "U⊥"), (along, "U∥")):
        ax.annotate(
            "",
            xy=(x0 + length * v_e, y0 + length * v_n),
            xytext=(x0, y0),
            arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "black"},
            zorder=8,
        )
        ax.text(
            x0 + 1.4 * length * v_e,
            y0 + 1.4 * length * v_n,
            label,
            ha="center",
            va="center",
            fontsize=6.4,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                  "pad": 0.8},
            zorder=9,
        )


def three_component_figure(
    fields: dict,
    geometry: dict,
    components: tuple[tuple[str, str, str, str | None], ...],
    *,
    stem: str,
    title: str,
    subtitle: str,
    footnote: str,
    frame_strike_deg: float | None = None,
    look_row: tuple[tuple[str, np.ndarray], ...] | None = None,
    weak_mask_sigma_cm: float | None = None,
) -> list[Path]:
    """Map figure: three components with formal σ, coverage, closure.

    ``components`` gives (field key, panel title, σ panel title, weak-note)
    per column. The column with a weak-note is stippled where its formal σ
    exceeds WEAK_SIGMA_STIPPLE_CM, or, when ``weak_mask_sigma_cm`` is set,
    masked above that σ instead. ``look_row`` prepends a row with the four
    LOS inputs (title, display-grid field in metres); the figure is then
    three rows instead of two.
    """
    centres_e, centres_n = lattice_axes()
    terrain_x, terrain_y, terrain_transform = display_grid()
    terrain = load_terrain(terrain_transform, (len(terrain_y), len(terrain_x)))

    plt.rcParams.update(RC_PARAMS)
    displacement_norm = TwoSlopeNorm(
        vmin=-DISPLACEMENT_LIMIT_CM, vcenter=0, vmax=DISPLACEMENT_LIMIT_CM
    )
    sigma_norm = Normalize(0, SIGMA_LIMIT_CM)
    closure_norm = Normalize(0, CLOSURE_LIMIT_SIGMA)
    coverage_cmap = ListedColormap(["#d79a32", "#287f8b"])
    coverage_norm = BoundaryNorm([2.5, 3.5, 4.5], coverage_cmap.N)
    lines = geometry_lines(geometry)
    visible_mask = enu_display_mask(centres_e, centres_n)
    n_rows = 2 if look_row is None else 3
    component_row = n_rows - 2
    sigma_row = n_rows - 1

    def draw_field(ax, values, cmap, norm):
        draw_enu_background(ax, terrain, terrain_x, terrain_y)
        return ax.pcolormesh(
            cell_edges(centres_e),
            cell_edges(centres_n),
            np.ma.masked_invalid(values),
            cmap=cmap,
            norm=norm,
            shading="flat",
            alpha=0.86,
            rasterized=True,
            zorder=3,
        )

    def format_map_axis(ax, row, column):
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.set_aspect("equal")
        ax.set_xticks([665_000, 675_000, 685_000, 695_000])
        ax.set_yticks([1_335_000, 1_345_000, 1_355_000, 1_365_000])
        if row == n_rows - 1:
            ax.set_xticklabels(["665", "675", "685", "695"])
            if column == 1:
                ax.set_xlabel("UTM 37N easting (km)", fontsize=7)
        else:
            ax.set_xticklabels([])
        if column == 0:
            ax.set_yticklabels(["1335", "1345", "1355", "1365"])
            ax.set_ylabel("UTM 37N northing (km)", fontsize=7)
        else:
            ax.set_yticklabels([])
        ax.tick_params(length=2.5, width=0.6, labelsize=6.4, pad=1.3)
        for spine in ax.spines.values():
            spine.set_linewidth(0.65)
            spine.set_color("#303030")

    def add_annotation(ax, text):
        ax.text(
            0.975, 0.025, text,
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.0,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84,
                  "pad": 1.3},
            zorder=9,
        )

    fig = plt.figure(figsize=(7.5, 5.5 if n_rows == 2 else 7.7))
    grid = fig.add_gridspec(
        n_rows, 4,
        left=0.075,
        # room for the latitude labels on the rightmost panel column
        right=0.925,
        bottom=0.21 if n_rows == 2 else 0.135,
        top=0.88 if n_rows == 2 else 0.905,
        hspace=0.32 if n_rows == 2 else 0.30,
        wspace=0.10,
    )
    axes = [[fig.add_subplot(grid[row, column]) for column in range(4)]
            for row in range(n_rows)]

    if look_row is not None:
        shade = hillshade(terrain)
        for column, (look_title, look_field) in enumerate(look_row):
            ax = axes[0][column]
            ax.imshow(
                blended_map(100 * look_field, shade, cmc.vik,
                            displacement_norm, alpha=0.86),
                extent=[X_MIN, X_MAX, Y_MIN, Y_MAX],
                origin="upper",
                interpolation="nearest",
                zorder=1,
            )
            ax.set_title(look_title, fontweight="bold", pad=3)

    for column, (key, panel_title, sigma_title, weak_note) in enumerate(components):
        sigma = 100 * fields[f"sigma_{key}"]
        values = 100 * fields[key]
        if weak_note is not None and weak_mask_sigma_cm is not None:
            values = np.where(sigma <= weak_mask_sigma_cm, values, np.nan)
        draw_field(axes[component_row][column], values, cmc.vik,
                   displacement_norm)
        axes[component_row][column].set_title(panel_title, fontweight="bold",
                                              pad=3)
        draw_field(axes[sigma_row][column], sigma, cmc.batlowW, sigma_norm)
        axes[sigma_row][column].set_title(sigma_title, fontweight="bold", pad=3)
        median = np.nanmedian(np.where(visible_mask, sigma, np.nan))
        add_annotation(axes[sigma_row][column], f"display median {median:.2f} cm")
        if weak_note is None:
            continue
        if weak_mask_sigma_cm is not None:
            note = f"{weak_note}\nmasked: σ > {weak_mask_sigma_cm:.0f} cm"
        else:
            note = f"{weak_note}\ndots: σ > {WEAK_SIGMA_STIPPLE_CM:.0f} cm"
            stipple = np.isfinite(fields[key]) & (sigma > WEAK_SIGMA_STIPPLE_CM)
            east_grid, north_grid = np.meshgrid(centres_e, centres_n)
            axes[component_row][column].plot(
                east_grid[stipple],
                north_grid[stipple],
                ls="none",
                marker=".",
                ms=0.9,
                color="#454545",
                alpha=0.65,
                zorder=4,
            )
        axes[component_row][column].text(
            0.975,
            0.025,
            note,
            transform=axes[component_row][column].transAxes,
            ha="right",
            va="bottom",
            fontsize=6.0,
            linespacing=1.35,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84,
                  "pad": 1.3},
            zorder=9,
        )

    draw_field(
        axes[component_row][3],
        np.where(fields["n_geometries"] >= 3, fields["n_geometries"], np.nan),
        coverage_cmap,
        coverage_norm,
    )
    axes[component_row][3].set_title("Geometry coverage", fontweight="bold",
                                     pad=3)
    visible_three = int(((fields["n_geometries"] == 3) & visible_mask).sum())
    visible_four = int(((fields["n_geometries"] == 4) & visible_mask).sum())
    add_annotation(
        axes[component_row][3],
        f"display: {visible_three:,} three · {visible_four:,} four",
    )
    draw_field(axes[sigma_row][3], fields["closure"], cmc.lajolla, closure_norm)
    axes[sigma_row][3].set_title("Four-geometry closure", fontweight="bold",
                                 pad=3)
    closure_visible = visible_mask & np.isfinite(fields["closure"])
    add_annotation(
        axes[sigma_row][3],
        f"display median {np.nanmedian(fields['closure'][closure_visible]):.2f} σ",
    )

    panel_index = 0
    for row in range(n_rows):
        for column in range(4):
            ax = axes[row][column]
            format_map_axis(ax, row, column)
            if row == 0 or column == 3:
                add_lonlat_axes(ax, top=row == 0, right=column == 3)
            draw_geometry(ax, lines)
            add_panel_label(ax, panel_index)
            if row == 0 and column == 0:
                add_scale_and_north(ax)
            if row == component_row and column == 0 \
                    and frame_strike_deg is not None:
                draw_frame_axes(ax, frame_strike_deg)
            panel_index += 1

    fig.canvas.draw()

    def add_horizontal_colorbar(axes_row, y, mappable, label, ticks,
                                ticklabels=None, *, label_top=False):
        left = axes_row[0].get_position().x0
        right = axes_row[-1].get_position().x1
        cax = fig.add_axes([left, y, right - left, 0.015])
        bar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        bar.set_ticks(ticks)
        if ticklabels is not None:
            bar.set_ticklabels(ticklabels)
        bar.set_label(label, fontsize=6.8, labelpad=1.5)
        if label_top:
            bar.ax.xaxis.set_label_position("top")
        bar.ax.tick_params(labelsize=6.2, length=2, pad=1.2)

    # The displacement bar sits in the gap under the first row: in the
    # two-row layout that is the component row, in the three-row layout the
    # look row, which shares the same norm.
    top_colorbar_y = (min(ax.get_position().y0 for ax in axes[0])
                      - (0.038 if n_rows == 2 else 0.034))
    add_horizontal_colorbar(
        axes[0][:3],
        top_colorbar_y,
        ScalarMappable(norm=displacement_norm, cmap=cmc.vik),
        "displacement (cm; clipped)",
        [-60, -30, 0, 30, 60],
        ["−60", "−30", "0", "+30", "+60"],
        label_top=True,
    )
    coverage_position = axes[component_row][3].get_position()
    coverage_legend_y = (top_colorbar_y + 0.0075 if n_rows == 2
                         else coverage_position.y0 - 0.020)
    fig.legend(
        handles=[
            Patch(facecolor=coverage_cmap.colors[0], edgecolor="#303030", lw=0.4,
                  label="3 geometries"),
            Patch(facecolor=coverage_cmap.colors[1], edgecolor="#303030", lw=0.4,
                  label="4 geometries"),
        ],
        loc="center",
        bbox_to_anchor=(
            0.5 * (coverage_position.x0 + coverage_position.x1),
            coverage_legend_y,
        ),
        bbox_transform=fig.transFigure,
        ncols=2,
        frameon=False,
        fontsize=6.2,
        handlelength=1.15,
        handleheight=1.0,
        columnspacing=0.9,
        handletextpad=0.45,
    )
    bottom_colorbar_y = 0.085 if n_rows == 2 else 0.060
    add_horizontal_colorbar(
        axes[sigma_row][:3],
        bottom_colorbar_y,
        ScalarMappable(norm=sigma_norm, cmap=cmc.batlowW),
        "formal σ (cm; clipped)",
        [0, 2, 4, 6, 8],
        label_top=n_rows == 3,
    )
    add_horizontal_colorbar(
        [axes[sigma_row][3]],
        bottom_colorbar_y,
        ScalarMappable(norm=closure_norm, cmap=cmc.lajolla),
        "weighted RMS (σ units; clipped)",
        [0, 1, 2],
        ["0", "1", "2"],
        label_top=n_rows == 3,
    )

    fig.suptitle(title, x=0.5, y=0.992, fontweight="bold")
    fig.text(
        0.5,
        0.947,
        subtitle,
        ha="center",
        va="bottom",
        fontsize=7.0,
        color="#3a3a3a",
    )
    fig.text(
        0.5,
        0.006 if n_rows == 2 else 0.002,
        footnote,
        ha="center",
        va="bottom",
        fontsize=6.0,
        color="#4a4a4a",
    )

    png = RESULTS / f"{stem}.png"
    pdf = RESULTS / f"{stem}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def build_decomposition_figure() -> list[Path]:
    _, geometry, _, fields, _ = decomposition_inputs()
    return three_component_figure(
        fields,
        geometry,
        (
            ("e", "East", "Formal σ east", None),
            ("n", "North", "Formal σ north", "weak LOS sensitivity"),
            ("u", "Up", "Formal σ up", None),
        ),
        stem="summary_decomposition_enu",
        title="Manda Hararo 2026 — free E–N–U decomposition",
        subtitle=(
            "reprocessed 30-m (GSLC) inputs · 1.28-km cells · weighted solve on "
            "≥3 of 4 LOS geometries · positive east, north, and up"
        ),
        footnote=(
            "No smoothing; gray = no solution; in the closure panel gray also "
            "covers three-geometry cells, which have no closure value.\n"
            "Formal σ uses far-field scatter only. Red/dashed lines show the "
            "mapped dike and east-ridge context, not a solve constraint."
        ),
    )


def build_dike_frame_figure() -> list[Path]:
    summary, geometry, model, fields, cells = decomposition_inputs()
    strike = float(model[2])
    perp_azimuth = (strike - 90.0) % 360.0
    par_azimuth = (strike + 180.0) % 360.0
    if not np.isclose(perp_azimuth, float(summary["perp_azimuth_deg"]), atol=1e-6):
        raise ValueError(
            f"strike frame changed: figure azimuth {perp_azimuth}, "
            f"saved {summary['perp_azimuth_deg']}"
        )
    frame = dike_frame(fields, strike)

    # The four decay-normalized LOS inputs, reprojected to the display grid.
    factors = {c["name"]: c["decay_factor"] for c in cells}
    xs, ys, transform = display_grid()
    shape = (len(ys), len(xs))
    X, Y = np.meshgrid(xs, ys)
    points = np.column_stack([Y.ravel(), X.ravel()])
    look_row = []
    for look in LOOKS:
        if look["kind"] == "nisar":
            observation, _ = load_nisar_look(
                look, factors[look["label"]], transform, shape, points
            )
        else:
            observation, _ = load_s1_look(
                look, model, factors[look["label"]], transform, shape
            )
        look_row.append((look["label"], observation))
        print(f"loaded {look['label']}")

    return three_component_figure(
        frame,
        geometry,
        (
            ("perp", "Dike-perpendicular U⊥", "Formal σ U⊥", None),
            ("par", "Dike-parallel U∥", "Formal σ U∥", "weak LOS sensitivity"),
            ("u", "Up", "Formal σ up", None),
        ),
        stem="summary_decomposition_dike_frame",
        title="Manda Hararo 2026 — dike-frame decomposition",
        subtitle=(
            "four LOS inputs (top row; positive toward the satellite, "
            "decay-normalized to 18 July)\n"
            f"free solve rotated to strike {strike:.1f}° · U⊥ positive "
            f"toward {perp_azimuth:.1f}° (ENE) · U∥ positive toward "
            f"{par_azimuth:.1f}° (NNW, near-north) · up positive"
        ),
        footnote=(
            "Top row: the LOS inputs of the solve. The components are an "
            "exact rotation of the free E–N–U solve; σ includes the "
            "east–north cross-covariance.\n"
            "U∥ is nearly north–south, which near-polar SAR observes weakly; "
            "U∥ cells with formal σ > 5 cm are masked. No smoothing; gray = "
            "no solution\n"
            "(three-geometry cells in the closure panel). Red/dashed lines "
            "show the mapped dike and east-ridge context, not a solve "
            "constraint."
        ),
        frame_strike_deg=strike,
        look_row=tuple(look_row),
        weak_mask_sigma_cm=5.0,
    )


# ---------------------------------------------------------------------------
# Figure 3: the saved source in physical three-dimensional geometry.


def validate_source(
    geometry: dict, magnitudes: np.ndarray, payload: dict, polygons: np.ndarray
) -> None:
    if polygons.shape != (len(magnitudes), 4, 3):
        raise ValueError("source patch corners are incomplete")
    dike = geometry["kinds"] == "dike"
    east = geometry["kinds"] == "fault_east"
    if int(dike.sum()) != int(np.prod(geometry["dike_shape"])):
        raise ValueError("saved dike patch count changed")
    if int(east.sum()) != int(geometry["fault_masks"]["east"].sum()):
        raise ValueError("saved east-fault patch count changed")
    if np.any(geometry["kinds"] == "fault_west"):
        raise ValueError("the trace-supported fit unexpectedly contains west-fault patches")
    if not np.isclose(magnitudes[dike].max(), payload["dike"]["max_opening_m"],
                      atol=1e-12):
        raise ValueError("saved dike maximum does not match its summary")
    if not np.isclose(magnitudes[east].max(),
                      payload["dike"]["fault_east_max_slip_m"], atol=1e-12):
        raise ValueError("saved fault maximum does not match its summary")


def local_surface(
    geometry: dict,
    terrain_x: np.ndarray,
    terrain_y: np.ndarray,
    terrain: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    s = np.linspace(-15_000, 15_000, 61)
    p = np.linspace(-3_500, 5_000, 35)
    S, P = np.meshgrid(s, p)
    strike = np.radians(geometry["strike"])
    centre_e, centre_n = geometry["trace_centre"]
    east = centre_e + S * np.sin(strike) - P * np.cos(strike)
    north = centre_n + S * np.cos(strike) + P * np.sin(strike)
    interpolator = RegularGridInterpolator(
        (terrain_y[::-1], terrain_x),
        terrain[::-1],
        bounds_error=False,
        fill_value=np.nan,
    )
    values = interpolator(np.column_stack([north.ravel(), east.ravel()]))
    return S / 1_000, P / 1_000, values.reshape(S.shape)


def build_source_figure() -> list[Path]:
    payload, geometry, magnitudes = load_mapped_gslc_fit()
    polygons = patch_corners_local(geometry)
    validate_source(geometry, magnitudes, payload, polygons)
    section_s = peak_fault_section_s_m(geometry, magnitudes) / 1_000

    terrain_x, terrain_y, terrain_transform = display_grid()
    terrain = load_terrain(terrain_transform, (len(terrain_y), len(terrain_x)))

    plt.rcParams.update(RC_PARAMS)
    dike_norm = Normalize(0, DIKE_LIMIT_M)
    fault_norm = Normalize(0, FAULT_LIMIT_M)
    dike = geometry["kinds"] == "dike"
    fault = geometry["kinds"] == "fault_east"

    fig = plt.figure(figsize=(7.5, 6.6))
    grid = fig.add_gridspec(
        3, 3,
        width_ratios=[1, 1, 0.76],
        height_ratios=[1.0, 1.0, 0.62],
        left=0.035,
        right=0.97,
        bottom=0.155,
        top=0.90,
        hspace=0.30,
        wspace=0.13,
    )
    main_ax = fig.add_subplot(grid[:2, :2], projection="3d")
    plan_ax = fig.add_subplot(grid[0, 2])
    section_ax = fig.add_subplot(grid[1:, 2])
    strip_ax = fig.add_subplot(grid[2, :2])

    # 3-D view.
    S, P, local_terrain = local_surface(geometry, terrain_x, terrain_y, terrain)
    shade = hillshade(local_terrain)
    surface = main_ax.plot_surface(
        S, P, np.zeros_like(S),
        facecolors=plt.cm.gray(0.72 + 0.20 * shade),
        rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False,
        # low alpha keeps the terrain context without darkening the dike wall
        alpha=0.16, zorder=0,
    )
    surface.set_rasterized(True)
    dike_collection = Poly3DCollection(
        polygons[dike],
        facecolors=cmc.lajolla_r(dike_norm(magnitudes[dike])),
        edgecolors=(0.25, 0.25, 0.25, 0.42),
        linewidths=0.28,
        alpha=0.98,
        zsort="average",
    )
    fault_collection = Poly3DCollection(
        polygons[fault],
        facecolors=FAULT_CMAP(fault_norm(magnitudes[fault])),
        # drawn after the dike so the near-side fault tint stays visible
        edgecolors=(0.04, 0.18, 0.33, 0.55),
        linewidths=0.40,
        alpha=0.78,
        zsort="average",
    )
    main_ax.add_collection3d(dike_collection)
    main_ax.add_collection3d(fault_collection)

    dike_s = np.asarray(geometry["dike_s_edges_m"]) / 1_000
    fault_s = np.asarray(geometry["fault_s_centres_m"]) / 1_000
    fault_p = np.asarray(geometry["fault_column_offsets_m"]["east"]) / 1_000
    main_ax.plot(dike_s, np.zeros_like(dike_s), np.zeros_like(dike_s),
                 color="#b21f2d", lw=2.0, zorder=9)
    main_ax.plot(fault_s, fault_p, np.zeros_like(fault_s),
                 color="#17466f", lw=1.6, ls=(0, (4, 2)), zorder=9)
    fault_polygons = polygons[fault]
    fault_centres_s = fault_polygons[:, :, 0].mean(axis=1)
    bottom_s, bottom_p, bottom_z = [], [], []
    for s_centre in fault_s:
        selected = np.isclose(fault_centres_s, s_centre, atol=1e-6)
        if not selected.any():
            continue
        column = fault_polygons[selected].reshape(-1, 3)
        z_max = column[:, 2].max()
        at_bottom = np.isclose(column[:, 2], z_max, atol=1e-6)
        bottom_s.append(s_centre)
        bottom_p.append(float(column[at_bottom, 1].mean()))
        bottom_z.append(float(z_max))
    main_ax.plot(bottom_s, bottom_p, bottom_z, color="#17466f", lw=0.72,
                 alpha=0.74, zorder=9)
    for end in (0, -1):
        main_ax.plot(
            [bottom_s[end], bottom_s[end]],
            [fault_p[end], bottom_p[end]],
            [0, bottom_z[end]],
            color="#17466f", lw=0.68, alpha=0.74, zorder=9,
        )

    main_ax.set_xlim(-15, 15)
    main_ax.set_ylim(-3.5, 5.0)
    main_ax.set_zlim(10, 0)
    main_ax.set_xticks([-10, 0, 10])
    main_ax.set_yticks([-2, 0, 2, 4])
    main_ax.set_zticks([0, 2, 4, 6, 8, 10])
    main_ax.set_xlabel("along strike, NNW → SSE (km)", labelpad=6)
    main_ax.set_ylabel("east of dike (km)", labelpad=5)
    # mplot3d pushes the z-label far outside this projection; the negative
    # labelpad pulls it back beside its own ticks.
    main_ax.zaxis.set_rotate_label(False)
    main_ax.set_zlabel("depth (km)", labelpad=-4, rotation=90, va="center")
    main_ax.tick_params(labelsize=6.3, pad=0)
    main_ax.set_box_aspect((30, 8.5, 10), zoom=1.05)
    main_ax.set_proj_type("ortho")
    main_ax.view_init(elev=33, azim=123)
    for axis in (main_ax.xaxis, main_ax.yaxis, main_ax.zaxis):
        axis.pane.set_facecolor((1, 1, 1, 0))
        axis.pane.set_edgecolor((0.65, 0.65, 0.65, 0.55))
        axis._axinfo["grid"]["color"] = (0.55, 0.55, 0.55, 0.22)
        axis._axinfo["grid"]["linewidth"] = 0.4
    main_ax.set_title("Physical source geometry", fontweight="bold", pad=1)
    main_ax.text2D(0.015, 0.965, "(a)", transform=main_ax.transAxes,
                   fontsize=7.4, fontweight="bold")
    main_ax.text2D(
        0.02, 0.055,
        f"dike: {int(dike.sum())} opening patches\n"
        f"east fault: {int(fault.sum())} retained normal-slip patches",
        transform=main_ax.transAxes, ha="left", va="bottom", fontsize=6.0,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
    )
    main_ax.text2D(
        0.985, 0.02, "no vertical exaggeration (box 1:1:1)",
        transform=main_ax.transAxes, ha="right", va="bottom",
        fontsize=6.0, color="#4a4a4a",
    )

    # Plan view.
    plan_shade = hillshade(terrain)
    plan_ax.imshow(
        plan_shade,
        cmap="gray",
        vmin=0.70,
        vmax=1.0,
        extent=[terrain_x[0] - 100, terrain_x[-1] + 100,
                terrain_y[-1] - 100, terrain_y[0] + 100],
        origin="upper",
        interpolation="nearest",
    )
    levels = np.arange(
        np.ceil(np.nanmin(terrain) / 100) * 100,
        np.floor(np.nanmax(terrain) / 100) * 100 + 1,
        200,
    )
    plan_ax.contour(terrain_x, terrain_y, terrain, levels=levels,
                    colors="#303030", linewidths=0.28, alpha=0.30)
    draw_geometry(plan_ax, geometry_lines(geometry), legend=True)
    plan_ax.set_xlim(PLAN_X_MIN, PLAN_X_MAX)
    plan_ax.set_ylim(PLAN_Y_MIN, PLAN_Y_MAX)
    plan_ax.set_aspect("equal")
    plan_ax.set_xticks([675_000, 685_000, 695_000])
    plan_ax.set_xticklabels(["675", "685", "695"])
    plan_ax.set_yticks([1_340_000, 1_350_000, 1_360_000])
    plan_ax.set_yticklabels(["1340", "1350", "1360"])
    plan_ax.tick_params(labelsize=6.0, length=2, width=0.6, pad=1)
    plan_ax.set_xlabel("UTM 37N easting (km)", fontsize=6.2, labelpad=1)
    plan_ax.set_ylabel("northing (km)", fontsize=6.2, labelpad=1)
    bar_x, bar_y = 673_500, 1_337_000
    plan_ax.plot([bar_x, bar_x + 5_000], [bar_y, bar_y], color="black", lw=2.0)
    plan_ax.text(bar_x + 2_500, bar_y + 650, "5 km", ha="center", fontsize=6.0)
    plan_ax.annotate(
        "N",
        xy=(673_800, 1_363_800),
        xytext=(673_800, 1_359_700),
        ha="center",
        fontsize=6.2,
        fontweight="bold",
        arrowprops={"arrowstyle": "-|>", "lw": 0.9},
    )
    plan_ax.set_title("Plan view", fontweight="bold", pad=2)
    plan_ax.text(0.02, 0.97, "(b)", transform=plan_ax.transAxes, ha="left",
                 va="top", fontsize=7.2, fontweight="bold",
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                       "pad": 1.2})

    # Cross-section at the peak-slip column.
    kinds = geometry["kinds"]
    fault_indices = np.flatnonzero(kinds == "fault_east")
    peak_fault_index = fault_indices[np.argmax(magnitudes[fault_indices])]
    centres_s = polygons[:, :, 0].mean(axis=1)
    spacing = float(np.median(np.diff(geometry["dike_s_centres"]))) / 2_000
    selected_dike = np.flatnonzero(
        (kinds == "dike") & (np.abs(centres_s - section_s) < spacing)
    )
    selected_fault = np.flatnonzero(
        (kinds == "fault_east") & (np.abs(centres_s - section_s) < spacing)
    )
    section_ax.axhline(0, color="#3d3d3d", lw=1.0)
    for index in selected_dike:
        z = polygons[index, :, 2]
        color = cmc.lajolla_r(dike_norm(magnitudes[index]))
        section_ax.plot([0, 0], [z.min(), z.max()], color="#484848", lw=6.5,
                        solid_capstyle="butt")
        section_ax.plot([0, 0], [z.min(), z.max()], color=color, lw=4.8,
                        solid_capstyle="butt")
    for index in selected_fault:
        p = polygons[index, :, 1]
        z = polygons[index, :, 2]
        order = np.argsort(z)
        low, high = order[0], order[-1]
        color = FAULT_CMAP(fault_norm(magnitudes[index]))
        section_ax.plot([p[low], p[high]], [z[low], z[high]], color="#484848",
                        lw=6.5, solid_capstyle="butt")
        section_ax.plot([p[low], p[high]], [z[low], z[high]], color=color,
                        lw=4.8, solid_capstyle="butt")
    section_ax.set_xlim(-2.5, 4.5)
    section_ax.set_ylim(10, -0.4)
    section_ax.set_aspect("equal")
    section_ax.set_xticks([-2, 0, 2, 4])
    section_ax.set_yticks([0, 2, 4, 6, 8, 10])
    section_ax.tick_params(labelsize=6.0, length=2, width=0.6, pad=1)
    section_ax.set_xlabel("east of dike (km)", fontsize=6.2, labelpad=1)
    section_ax.set_ylabel("depth (km)", fontsize=6.2, labelpad=1)
    section_ax.set_title(f"Section at s = {section_s:+.1f} km",
                         fontweight="bold", pad=2)
    section_ax.text(0.02, 0.97, "(c)", transform=section_ax.transAxes,
                    ha="left", va="top", fontsize=7.2, fontweight="bold")
    peak_p = float(polygons[peak_fault_index, :, 1].mean())
    peak_z = float(polygons[peak_fault_index, :, 2].mean())
    section_ax.annotate(
        "peak fault slip",
        xy=(peak_p, peak_z),
        xytext=(0.97, 0.46),
        textcoords="axes fraction",
        ha="right",
        va="top",
        fontsize=6.0,
        color="#3f3f3f",
        arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "#3f3f3f",
                    "shrinkA": 2, "shrinkB": 3},
    )
    section_ax.text(0.97, 0.03, "no vertical\nexaggeration",
                    transform=section_ax.transAxes, ha="right", va="bottom",
                    fontsize=6.0, color="#4a4a4a", linespacing=1.35)
    for spine in section_ax.spines.values():
        spine.set_linewidth(0.65)

    # Unfolded east-fault slip. The dipping patches project to axis-aligned
    # rectangles here because the along-strike patch edges are horizontal.
    for index in fault_indices:
        s = polygons[index, :, 0]
        z = polygons[index, :, 2]
        strip_ax.add_patch(
            plt.Rectangle(
                (s.min(), z.min()),
                s.max() - s.min(),
                z.max() - z.min(),
                facecolor=FAULT_CMAP(fault_norm(magnitudes[index])),
                edgecolor="#484848",
                lw=0.3,
            )
        )
    strip_ax.axhline(0, color="#3d3d3d", lw=0.8)
    strip_ax.axvline(section_s, color="#666666", lw=0.6, ls=(0, (3, 2)))
    max_depth = float(np.ceil(polygons[fault_indices, :, 2].max()))
    strip_ax.set_xlim(-15, 15)
    strip_ax.set_ylim(max_depth, -0.2)
    strip_ax.set_aspect("equal")
    strip_ax.set_xticks([-10, 0, 10])
    strip_ax.set_yticks(np.arange(0, max_depth + 1, 2))
    strip_ax.tick_params(labelsize=6.0, length=2, width=0.6, pad=1)
    strip_ax.set_xlabel("along strike, NNW → SSE (km)", fontsize=6.2, labelpad=1)
    strip_ax.set_ylabel("depth (km)", fontsize=6.2, labelpad=1)
    strip_ax.set_title("East-fault slip (unfolded)", fontweight="bold", pad=2)
    strip_ax.text(0.012, 0.94, "(d)", transform=strip_ax.transAxes, ha="left",
                  va="top", fontsize=7.2, fontweight="bold")
    strip_ax.text(
        0.012, 0.07,
        "no vertical exaggeration\ndashed line = section (c)",
        transform=strip_ax.transAxes, ha="left", va="bottom",
        fontsize=6.0, color="#4a4a4a", linespacing=1.35,
    )
    for spine in strip_ax.spines.values():
        spine.set_linewidth(0.65)

    def add_colorbar(position, norm, cmap, label, ticks):
        cax = fig.add_axes(position)
        bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                           orientation="horizontal")
        bar.set_ticks(ticks)
        bar.set_label(label, fontsize=6.8, labelpad=1.5)
        bar.ax.tick_params(labelsize=6.2, length=2, pad=1.2)

    add_colorbar([0.075, 0.088, 0.40, 0.014], dike_norm, cmc.lajolla_r,
                 "dike opening (m)", [0, 1, 2, 2.6])
    add_colorbar([0.555, 0.088, 0.36, 0.014], fault_norm, FAULT_CMAP,
                 "east-fault normal slip (m)", [0, 0.2, 0.4, 0.65])

    dike_summary = payload["dike"]
    potency = payload["potency_partition"]
    fig.suptitle(
        "Manda Hararo 2026 — physical distributed source",
        x=0.5,
        y=0.992,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.947,
        "trace-supported fit to the 30-m (GSLC) inputs · "
        f"strike {geometry['strike']:.1f}° · east fault dips "
        f"{geometry['fault_dip']:.0f}° toward the axis",
        ha="center",
        va="bottom",
        fontsize=7.1,
        color="#3a3a3a",
    )
    fig.text(
        0.5,
        0.010,
        f"Opening max {dike_summary['max_opening_m']:.2f} m · volume "
        f"{dike_summary['dike_volume_km3']:.3f} km³ · fault max "
        f"{dike_summary['fault_east_max_slip_m']:.2f} m · fault/dike potency "
        f"{potency['fault_to_dike_potency_ratio']:.2f} "
        "(area × slip over area × opening; not a moment-tensor ratio).\n"
        "The 2005-plane refit is in the Wright-2006 summary figures; the "
        "2006–2010 geometry and opening files are still awaited.",
        ha="center",
        va="bottom",
        fontsize=6.0,
        color="#4a4a4a",
        linespacing=1.35,
    )

    png = RESULTS / "summary_source_geometry.png"
    pdf = RESULTS / "summary_source_geometry.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


# ---------------------------------------------------------------------------
# Figure 4: unwrapped LOS and the signed cross-rift phase change.


def strike_coordinates(
    x: np.ndarray, y: np.ndarray, geometry: dict
) -> tuple[np.ndarray, np.ndarray]:
    X, Y = np.meshgrid(x, y)
    centre_e, centre_n = geometry["trace_centre"]
    strike = np.radians(geometry["strike"])
    d_e = X - centre_e
    d_n = Y - centre_n
    along = d_n * np.cos(strike) + d_e * np.sin(strike)
    perp = d_n * np.sin(strike) - d_e * np.cos(strike)
    return along, perp


def binned_profile(
    field: np.ndarray,
    along: np.ndarray,
    perp: np.ndarray,
    section_s_m: float,
    *,
    half_width_m: float = 500.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    edges = np.arange(-SECTION_P_LIMIT_M, SECTION_P_LIMIT_M + 100, 100.0)
    centres = 0.5 * (edges[:-1] + edges[1:])
    near = np.abs(along - section_s_m) < half_width_m
    median = np.full(len(centres), np.nan)
    p16 = np.full(len(centres), np.nan)
    p84 = np.full(len(centres), np.nan)
    indices = np.digitize(perp[near], edges)
    values = field[near]
    for index in range(1, len(edges)):
        selected = (indices == index) & np.isfinite(values)
        if selected.sum() >= 6:
            block = values[selected]
            median[index - 1] = np.median(block)
            p16[index - 1] = np.percentile(block, 16)
            p84[index - 1] = np.percentile(block, 84)
    return centres, median, p16, p84


def binned_alias_fraction(
    flag: np.ndarray,
    along: np.ndarray,
    perp: np.ndarray,
    section_s_m: float,
    *,
    half_width_m: float = 500.0,
) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(-SECTION_P_LIMIT_M, SECTION_P_LIMIT_M + 100, 100.0)
    centres = 0.5 * (edges[:-1] + edges[1:])
    near = np.abs(along - section_s_m) < half_width_m
    indices = np.digitize(perp[near], edges)
    values = flag[near].astype(float)
    fraction = np.zeros(len(centres))
    for index in range(1, len(edges)):
        selected = indices == index
        if selected.sum() >= 6:
            fraction[index - 1] = float(values[selected].mean())
    return centres, fraction


def ridge_offsets(traces: dict, section_s_m: float) -> dict[str, float]:
    s = np.asarray(traces["s_centres_m"], dtype=float)
    output = {}
    for side in ("east", "west"):
        values = np.array(
            [np.nan if value is None else value for value in traces[side]["offset_m"]],
            dtype=float,
        )
        finite = np.isfinite(values)
        if not finite.any():
            output[side] = np.nan
            continue
        index = np.argmin(np.abs(s[finite] - section_s_m))
        output[side] = float(values[finite][index])
    return output


def trace_lines(
    traces: dict, geometry: dict, section_s_m: float
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    centre_e, centre_n = geometry["trace_centre"]
    strike = np.radians(geometry["strike"])
    s = np.asarray(traces["s_centres_m"], dtype=float)
    lines = {}
    for side in ("east", "west"):
        p = np.array(
            [np.nan if value is None else value for value in traces[side]["offset_m"]],
            dtype=float,
        )
        lines[side] = (
            centre_e + s * np.sin(strike) - p * np.cos(strike),
            centre_n + s * np.cos(strike) + p * np.sin(strike),
        )
    dike_s = np.asarray(geometry["dike_s_edges_m"], dtype=float)
    lines["dike"] = (
        centre_e + dike_s * np.sin(strike),
        centre_n + dike_s * np.cos(strike),
    )
    p = np.array([-SECTION_P_LIMIT_M, SECTION_P_LIMIT_M])
    lines["section"] = (
        centre_e + section_s_m * np.sin(strike) - p * np.cos(strike),
        centre_n + section_s_m * np.cos(strike) + p * np.sin(strike),
    )
    return lines


def image_extent(x: np.ndarray, y: np.ndarray) -> list[float]:
    dx = float(abs(np.median(np.diff(x))))
    dy = float(abs(np.median(np.diff(y))))
    return [x.min() - dx / 2, x.max() + dx / 2, y.min() - dy / 2, y.max() + dy / 2]


def build_gradient_figure() -> list[Path]:
    traces = json.loads((RESULTS / "fault_traces.json").read_text())
    payload, geometry, magnitudes = load_mapped_gslc_fit()
    section_s_m = peak_fault_section_s_m(geometry, magnitudes)
    # |Δφ| = π per pixel corresponds to λ/4 of LOS change per pixel.
    alias_limit_m = WAVELENGTH_L / 4.0

    terrain_x, terrain_y, terrain_transform = display_grid()
    terrain = load_terrain(terrain_transform, (len(terrain_y), len(terrain_x)))
    terrain_shade = hillshade(terrain)

    looks = []
    for tag in INVERSION_PAIRS:
        interferogram = np.load(DATA_GSLC / f"{tag}_ifg.npz")
        unwrapped = np.load(DATA_GSLC / f"{tag}_unw.npz")
        if not (
            np.array_equal(interferogram["x"], unwrapped["x"])
            and np.array_equal(interferogram["y"], unwrapped["y"])
        ):
            raise ValueError(f"{tag}: wrapped and unwrapped grids differ")
        los = np.where(
            unwrapped["coherence"] > COHERENCE_MIN_LOS, unwrapped["los"], np.nan
        )
        looks.append(
            {
                "tag": tag,
                "label": NISAR_LABEL[tag],
                "x": interferogram["x"],
                "y": interferogram["y"],
                "phase": signed_cross_rift_increment(
                    interferogram["phase"],
                    interferogram["coherence"],
                    interferogram["x"],
                    interferogram["y"],
                    geometry["strike"],
                ),
                "los": los,
                "alias": aliasing_flag(
                    los, unwrapped["y"], geometry["strike"], alias_limit_m
                ),
            }
        )
        print(f"loaded {NISAR_LABEL[tag]}")

    plt.rcParams.update(RC_PARAMS)
    los_norm = TwoSlopeNorm(vmin=-LOS_LIMIT_CM, vcenter=0, vmax=LOS_LIMIT_CM)
    phase_norm = TwoSlopeNorm(
        vmin=-PHASE_INCREMENT_LIMIT_RAD, vcenter=0, vmax=PHASE_INCREMENT_LIMIT_RAD
    )
    lines = trace_lines(traces, geometry, section_s_m)
    section_ridges = ridge_offsets(traces, section_s_m)

    def draw_terrain_layer(ax):
        ax.imshow(
            terrain_shade,
            cmap="gray",
            vmin=0.70,
            vmax=1.0,
            extent=[terrain_x[0] - 100, terrain_x[-1] + 100,
                    terrain_y[-1] - 100, terrain_y[0] + 100],
            origin="upper",
            interpolation="nearest",
            zorder=0,
        )
        levels = np.arange(
            np.ceil(np.nanmin(terrain) / 100) * 100,
            np.floor(np.nanmax(terrain) / 100) * 100 + 1,
            200,
        )
        ax.contour(terrain_x, terrain_y, terrain, levels=levels,
                   colors="#303030", linewidths=0.28, alpha=0.24, zorder=2)

    def draw_lines(ax, *, labels=False):
        halo = [path_effects.Stroke(linewidth=2.3, foreground="white"),
                path_effects.Normal()]
        styles = {
            "dike": ("#b21f2d", "-", 1.05),
            "east": ("#202020", (0, (4, 2)), 0.95),
            "west": ("#2f4f6f", (0, (2.5, 1.5)), 1.05),
        }
        for name in ("dike", "east", "west"):
            color, style, width = styles[name]
            line = ax.plot(*lines[name], color=color, ls=style, lw=width, zorder=7)[0]
            line.set_path_effects(halo)
        section = ax.plot(*lines["section"], color="#2b2b2b", ls=(0, (3, 2)),
                          lw=0.95, zorder=6)[0]
        section.set_path_effects(halo)
        if labels:
            ax.text(lines["section"][0][0], lines["section"][1][0], "B",
                    ha="right", va="top", fontsize=6.3, fontweight="bold")
            ax.text(lines["section"][0][-1], lines["section"][1][-1], "B′",
                    ha="left", va="bottom", fontsize=6.3, fontweight="bold")

    def format_map(ax, row, column):
        ax.set_xlim(PLAN_X_MIN, PLAN_X_MAX)
        ax.set_ylim(PLAN_Y_MIN, PLAN_Y_MAX)
        ax.set_aspect("equal")
        ax.set_xticks([675_000, 685_000, 695_000])
        ax.set_yticks([1_340_000, 1_350_000, 1_360_000])
        if row == 1:
            ax.set_xticklabels(["675", "685", "695"])
            if column == 1:
                ax.set_xlabel("UTM 37N easting (km)", fontsize=7)
        else:
            ax.set_xticklabels([])
        if column == 0:
            ax.set_yticklabels(["1340", "1350", "1360"])
            ax.set_ylabel("UTM 37N northing (km)", fontsize=7)
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=6.4, length=2.5, width=0.6, pad=1.2)
        for spine in ax.spines.values():
            spine.set_linewidth(0.65)

    def add_map_annotation(ax, text):
        ax.text(
            0.025, 0.025, text,
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.0,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84,
                  "pad": 1.3},
            zorder=10,
        )

    fig = plt.figure(figsize=(7.5, 7.2))
    # The narrow third column holds the °N labels of the signed-gradient maps
    # and the profile's left axis text.
    grid = fig.add_gridspec(
        2, 4,
        width_ratios=[1, 1, 0.18, 1.24],
        left=0.07,
        right=0.93,
        bottom=0.16,
        top=0.885,
        hspace=0.12,
        wspace=0.15,
    )
    axes = [
        [fig.add_subplot(grid[row, column]) for column in (0, 1, 3)]
        for row in range(2)
    ]
    panel_index = 0
    for row, look in enumerate(looks):
        x, y = look["x"], look["y"]
        along, perp = strike_coordinates(x, y, geometry)

        los_ax = axes[row][0]
        draw_terrain_layer(los_ax)
        los_ax.imshow(
            100 * look["los"],
            cmap=cmc.vik,
            norm=los_norm,
            extent=image_extent(x, y),
            origin="upper",
            interpolation="nearest",
            alpha=0.84,
            zorder=3,
        )
        draw_lines(los_ax, labels=(row == 0))
        add_map_annotation(los_ax, look["label"])
        format_map(los_ax, row, 0)
        if row == 0:
            add_lonlat_axes(los_ax, top=True)
            los_ax.set_title("Unwrapped LOS", fontweight="bold", pad=12)
            bar_x, bar_y = 689_500.0, 1_336_200.0
            halo = [path_effects.Stroke(linewidth=3.2, foreground="white"),
                    path_effects.Normal()]
            bar = los_ax.plot([bar_x, bar_x + 5_000], [bar_y, bar_y],
                              color="black", lw=2.0, zorder=8)[0]
            bar.set_path_effects(halo)
            los_ax.text(
                bar_x + 2_500, bar_y + 500, "5 km",
                ha="center", va="bottom", fontsize=6.0,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75,
                      "pad": 0.8},
                zorder=8,
            )
        add_panel_label(los_ax, panel_index)
        panel_index += 1

        phase_ax = axes[row][1]
        draw_terrain_layer(phase_ax)
        phase_ax.imshow(
            look["phase"],
            cmap=cmc.broc,
            norm=phase_norm,
            extent=image_extent(x, y),
            origin="upper",
            interpolation="nearest",
            alpha=0.88,
            zorder=3,
        )
        veil = np.zeros((*look["alias"].shape, 4))
        veil[look["alias"]] = (1.0, 1.0, 1.0, VEIL_ALPHA)
        phase_ax.imshow(
            veil,
            extent=image_extent(x, y),
            origin="upper",
            interpolation="nearest",
            zorder=4,
        )
        draw_lines(phase_ax)
        add_map_annotation(
            phase_ax,
            "east: step-supported · west: ridge only"
            if row == 0
            else "same merged ridges in both geometries",
        )
        format_map(phase_ax, row, 1)
        add_lonlat_axes(phase_ax, top=(row == 0), right=True)
        if row == 0:
            phase_ax.set_title("Signed cross-rift phase change",
                               fontweight="bold", pad=12)
            phase_ax.text(
                0.975, 0.975, "white veil = aliased",
                transform=phase_ax.transAxes, ha="right", va="top", fontsize=6.0,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84,
                      "pad": 1.2},
                zorder=10,
            )
        add_panel_label(phase_ax, panel_index)
        panel_index += 1

        profile_ax = axes[row][2]
        distance, los_profile, _, _ = binned_profile(
            look["los"], along, perp, section_s_m
        )
        _, phase_median, phase_p16, phase_p84 = binned_profile(
            look["phase"], along, perp, section_s_m
        )
        _, alias_fraction = binned_alias_fraction(
            look["alias"], along, perp, section_s_m
        )
        aliased_bins = alias_fraction > ALIAS_FRACTION_CUT
        # Inside the aliased zone the binned medians are not interpretable; a
        # labelled gray band replaces the trace there.
        phase_median = np.where(aliased_bins, np.nan, phase_median)
        phase_p16 = np.where(aliased_bins, np.nan, phase_p16)
        phase_p84 = np.where(aliased_bins, np.nan, phase_p84)
        phase_profile_ax = profile_ax.twinx()
        bin_half_km = 0.05
        band_edges = []
        start = None
        for index, flagged in enumerate([*aliased_bins, False]):
            if flagged and start is None:
                start = index
            elif not flagged and start is not None:
                band_edges.append((distance[start], distance[index - 1]))
                start = None
        for low, high in band_edges:
            profile_ax.axvspan(
                low / 1_000 - bin_half_km,
                high / 1_000 + bin_half_km,
                color="#d4d4d4",
                alpha=0.55,
                zorder=0.5,
                linewidth=0,
            )
        if band_edges:
            widest = max(band_edges, key=lambda pair: pair[1] - pair[0])
            # A sliver band gets no in-panel label; the footnote explains the
            # gray band for both rows.
            if widest[1] - widest[0] >= 250.0:
                phase_profile_ax.text(
                    0.5 * (widest[0] + widest[1]) / 1_000,
                    0.262,
                    "aliased",
                    ha="center",
                    va="top",
                    fontsize=6.0,
                    color="#787878",
                    zorder=6,
                )
        profile_ax.plot(distance / 1_000, 100 * los_profile, color="#202020",
                        lw=1.15, label="LOS")
        clip_level = 0.985 * PHASE_INCREMENT_LIMIT_RAD
        phase_profile_ax.plot(
            distance / 1_000,
            np.clip(phase_median, -clip_level, clip_level),
            color="#b33a33",
            lw=0.95,
            label="phase Δ⊥",
        )
        phase_profile_ax.fill_between(
            distance / 1_000,
            np.clip(phase_p16, -clip_level, clip_level),
            np.clip(phase_p84, -clip_level, clip_level),
            color="#b33a33",
            alpha=0.15,
            linewidth=0,
        )
        clipped_high = phase_median > clip_level
        clipped_low = phase_median < -clip_level
        for mask, level, marker in (
            (clipped_high, clip_level, "^"),
            (clipped_low, -clip_level, "v"),
        ):
            if mask.any():
                phase_profile_ax.plot(
                    distance[mask] / 1_000,
                    np.full(int(mask.sum()), level),
                    marker=marker,
                    ms=2.3,
                    ls="none",
                    color="#b33a33",
                    zorder=5,
                )
        profile_ax.axvline(0, color="#b21f2d", lw=0.9)
        profile_ax.axvline(section_ridges["east"] / 1_000, color="#202020",
                           lw=0.85, ls=(0, (4, 2)))
        profile_ax.axvline(section_ridges["west"] / 1_000, color="#2f4f6f",
                           lw=0.95, ls=(0, (2.5, 1.5)))
        profile_ax.axhline(0, color="#777777", lw=0.45)
        phase_profile_ax.axhline(0, color="#b33a33", lw=0.35, alpha=0.35)
        profile_ax.set_xlim(-5, 5)
        profile_ax.set_ylim(-LOS_LIMIT_CM, LOS_LIMIT_CM)
        phase_profile_ax.set_ylim(-PHASE_INCREMENT_LIMIT_RAD,
                                  PHASE_INCREMENT_LIMIT_RAD)
        profile_ax.set_yticks([-60, -30, 0, 30, 60])
        phase_profile_ax.set_yticks([-0.3, -0.15, 0, 0.15, 0.3])
        profile_ax.tick_params(labelsize=6.4, length=2.5, width=0.6, pad=1.2)
        phase_profile_ax.tick_params(labelsize=6.1, length=2.3, width=0.55,
                                     pad=1.0, colors="#8b322f")
        profile_ax.set_ylabel("LOS (cm)", fontsize=6.5, labelpad=1)
        phase_profile_ax.set_ylabel("phase Δ⊥ (rad/30 m)", fontsize=6.5,
                                    labelpad=2.5, color="#8b322f")
        if row == 1:
            profile_ax.set_xlabel("distance from dike, ENE positive (km)", fontsize=7)
        else:
            profile_ax.set_xticklabels([])
            profile_ax.set_title(
                f"Profile BB′ at s = {section_s_m / 1_000:+.1f} km",
                fontweight="bold",
                pad=3,
            )
            handles = [profile_ax.lines[0], phase_profile_ax.lines[0]]
            profile_ax.legend(
                handles, ["LOS", "phase Δ⊥"],
                loc="upper left", frameon=False, fontsize=6.1, ncols=2,
                handlelength=1.8, columnspacing=0.9,
            )
        step = traces["los_step_across_ridges"]["per_geometry"][look["tag"]]
        add_map_annotation(
            profile_ax,
            f"ridge-median step: east {step['east']['median_step_cm']:+.1f} cm "
            f"· west {step['west']['median_step_cm']:+.1f} cm",
        )
        profile_ax.grid(color="#a8a8a8", alpha=0.23, lw=0.45)
        for spine in profile_ax.spines.values():
            spine.set_linewidth(0.65)
        for spine in phase_profile_ax.spines.values():
            spine.set_linewidth(0.65)
        profile_ax.set_box_aspect(
            (PLAN_Y_MAX - PLAN_Y_MIN) / (PLAN_X_MAX - PLAN_X_MIN) / 1.24
        )
        add_panel_label(profile_ax, panel_index, right=True)
        panel_index += 1

    fig.canvas.draw()

    def add_colorbar(ax, norm, cmap, y, label, ticks, ticklabels):
        position = ax.get_position()
        cax = fig.add_axes([position.x0, y, position.width, 0.015])
        bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                           orientation="horizontal")
        bar.set_ticks(ticks)
        bar.set_ticklabels(ticklabels)
        bar.set_label(label, fontsize=6.7, labelpad=1.5)
        bar.ax.tick_params(labelsize=6.1, length=2, pad=1.1)

    add_colorbar(
        axes[1][0], los_norm, cmc.vik, 0.104,
        "LOS displacement (cm; clipped)", [-70, 0, 70], ["−70", "0", "+70"],
    )
    add_colorbar(
        axes[1][1], phase_norm, cmc.broc, 0.104,
        "signed phase change (rad per 30 m; clipped)",
        [-0.3, 0, 0.3], ["−0.30", "0", "+0.30"],
    )
    profile_position = axes[1][2].get_position()
    fig.text(
        profile_position.x0,
        0.082,
        "red = dike · black dashed = step-supported east ridge\n"
        "blue dotted = west ridge; its LOS-step sign is inconsistent",
        ha="left",
        va="center",
        fontsize=6.0,
        color="#4a4a4a",
    )

    fig.suptitle(
        "Manda Hararo 2026 — signed phase-gradient evidence for shallow faulting",
        x=0.5,
        y=0.992,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.947,
        "30-m Goldstein-filtered interferograms · signed wrapped-phase "
        "difference, positive ENE across the rift",
        ha="center",
        va="bottom",
        fontsize=7.0,
        color="#3a3a3a",
    )
    fig.text(
        0.5,
        0.002,
        "Signed rift-perpendicular wrapped-phase differences after Sandwell & "
        "Price (1998), doi:10.1029/1998JB900008; cf. Xu et al. (2020), "
        "doi:10.1126/science.abd1690.\n"
        "The ridge picker (results/fault_traces.json) uses the unsigned "
        "gradient magnitude |∇φ|.\n"
        "The signed difference aliases where the cross-rift phase rate exceeds "
        "π per 30-m pixel (λ/4 of LOS per pixel); the white map veil and the "
        "gray profile band\n"
        "mark that zone, estimated from the unwrapped LOS. The "
        "forward-difference stencil shifts the signed map ≈15 m WNW of true "
        "gradient positions.",
        ha="center",
        va="bottom",
        fontsize=6.0,
        color="#4a4a4a",
        linespacing=1.3,
    )

    png = RESULTS / "summary_phase_gradient.png"
    pdf = RESULTS / "summary_phase_gradient.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


# ---------------------------------------------------------------------------
# Figure 5: dike-induced stress on the east-fault receivers.


def stress_fields(payload: dict) -> dict:
    geometry = geometry_from_json(payload["geometry"])
    magnitudes = np.asarray(payload["dislocation_magnitudes_m"], dtype=float)
    if magnitudes.shape != geometry["kinds"].shape:
        raise ValueError("saved stress source vector does not match its geometry")

    dike = geometry["kinds"] == "dike"
    source_patches = geometry["patches"][dike]
    source_dislocations = geometry["disl"][dike] * magnitudes[dike, None]
    _, perpendicular = _strike_vectors(geometry["strike"])

    receiver = fault_profile_receivers(geometry, MAIN_DIP)["east"]
    keep = np.abs(receiver["offsets_m"]) > DIKE_STANDOFF_M
    points = receiver["points"][keep]
    tensor = stress_at(source_patches, source_dislocations, points)
    sigma_n, cfs = resolve_normal_slip(tensor, perpendicular, +1.0, MAIN_DIP)
    sigma_n = sigma_n / 1e6
    cfs = cfs / 1e6

    depth_edges = np.arange(
        0.0,
        float(np.asarray(geometry["fault_z_edges_m"])[-1]) + PROFILE_DZ_M,
        PROFILE_DZ_M,
    )
    depth_centres = 0.5 * (depth_edges[:-1] + depth_edges[1:])
    along_centres = np.asarray(geometry["fault_s_centres_m"], dtype=float)
    shape = (len(depth_centres), len(along_centres))
    maps = {key: np.full(shape, np.nan) for key in ("sigma_n", "cfs")}
    depths = receiver["depths_m"][keep]
    columns = receiver["columns"][keep]
    rows = np.rint((depths - PROFILE_DZ_M / 2) / PROFILE_DZ_M).astype(int)
    if np.any(rows < 0) or np.any(rows >= len(depth_centres)):
        raise ValueError("profile receiver depth is outside the map grid")
    maps["sigma_n"][rows, columns] = sigma_n
    maps["cfs"][rows, columns] = cfs

    profiles = {"depth_km": [], "sigma_n": [], "sigma_p16": [], "sigma_p84": [],
                "cfs": []}
    for row, depth in enumerate(depth_centres):
        valid = np.isfinite(maps["sigma_n"][row])
        if valid.sum() < 4:
            continue
        profiles["depth_km"].append(float(depth / 1_000))
        profiles["sigma_n"].append(float(np.median(maps["sigma_n"][row, valid])))
        profiles["sigma_p16"].append(
            float(np.percentile(maps["sigma_n"][row, valid], 16))
        )
        profiles["sigma_p84"].append(
            float(np.percentile(maps["sigma_n"][row, valid], 84))
        )
        profiles["cfs"].append(float(np.median(maps["cfs"][row, valid])))
    profiles = {key: np.asarray(value) for key, value in profiles.items()}

    dip_summary = {}
    for dip in RECEIVER_DIPS:
        patch_receiver = fault_patch_receivers(geometry, dip)["east"]
        selected = np.abs(patch_receiver["offsets_m"]) > DIKE_STANDOFF_M
        patch_points = patch_receiver["points"][selected]
        patch_tensor = stress_at(source_patches, source_dislocations, patch_points)
        patch_sigma, patch_cfs = resolve_normal_slip(
            patch_tensor, perpendicular, +1.0, dip
        )
        dip_summary[f"{dip:.0f}"] = side_stats(patch_sigma / 1e6, patch_cfs / 1e6)

    return {
        "geometry": geometry,
        "receiver_points": points,
        "along_edges_km": np.asarray(geometry["fault_s_edges_m"]) / 1_000,
        "depth_edges_km": depth_edges / 1_000,
        "maps": maps,
        "profiles": profiles,
        "dip_summary": dip_summary,
    }


def map_space_stress(payload: dict) -> dict:
    """Map-space stress at a fixed depth on the trace-supported east-fault
    receiver orientation (dike strike, 65° dip toward the axis, normal rake,
    friction from ``mh26.stress``)."""
    geometry = geometry_from_json(payload["geometry"])
    magnitudes = np.asarray(payload["dislocation_magnitudes_m"], dtype=float)
    dike = geometry["kinds"] == "dike"
    source_patches = geometry["patches"][dike]
    source_dislocations = geometry["disl"][dike] * magnitudes[dike, None]
    along, perpendicular = _strike_vectors(geometry["strike"])

    xs, ys, _ = display_grid()
    easting, northing = np.meshgrid(xs, ys)
    points = np.column_stack(
        [northing.ravel(), easting.ravel(), np.full(easting.size, MAP_DEPTH_M)]
    )
    tensor = stress_at(source_patches, source_dislocations, points)
    sigma_n, cfs = resolve_normal_slip(tensor, perpendicular, +1.0, MAIN_DIP)
    shape = easting.shape
    fields = {
        "sigma_n": (sigma_n / 1e6).reshape(shape),
        "cfs": (cfs / 1e6).reshape(shape),
    }

    e_c, n_c = geometry["trace_centre"]
    delta_n = northing - n_c
    delta_e = easting - e_c
    along_m = delta_n * along[0] + delta_e * along[1]
    perp_m = delta_n * perpendicular[0] + delta_e * perpendicular[1]
    dike_s = np.asarray(geometry["dike_s_edges_m"], dtype=float)
    strip = (
        (np.abs(perp_m) < DIKE_STANDOFF_M)
        & (along_m >= dike_s[0])
        & (along_m <= dike_s[-1])
    )
    for field in fields.values():
        field[strip] = np.nan
    return {"geometry": geometry, "fields": fields}


def compare_stress_fields(gunw: dict, gslc: dict) -> dict:
    if not np.allclose(gunw["receiver_points"], gslc["receiver_points"]):
        raise ValueError("GUNW and GSLC stress cases do not use identical receivers")
    metrics = {}
    differences = {}
    for key in ("sigma_n", "cfs"):
        first = gunw["maps"][key]
        second = gslc["maps"][key]
        if not np.array_equal(np.isfinite(first), np.isfinite(second)):
            raise ValueError(f"{key}: GUNW and GSLC receiver masks differ")
        valid = np.isfinite(first)
        difference = second - first
        differences[key] = difference
        metrics[key] = {
            "relative_l2": float(
                np.linalg.norm(difference[valid]) / np.linalg.norm(first[valid])
            ),
            "sign_agreement_fraction": float(
                np.mean(np.sign(second[valid]) == np.sign(first[valid]))
            ),
        }
    return {"differences": differences, "metrics": metrics}


def build_stress_figure() -> list[Path]:
    if not selftest():
        raise SystemExit("stress derivative self-check failed")
    gunw_payload = json.loads((RESULTS / "distributed_fit_mapped.json").read_text())
    gslc_payload = json.loads(
        (RESULTS / "distributed_fit_mapped_gslc.json").read_text()
    )
    if gunw_payload["geometry_sha256"] != gslc_payload["geometry_sha256"]:
        raise ValueError("GUNW and GSLC fits use different receiver geometry")
    gunw = stress_fields(gunw_payload)
    gslc = stress_fields(gslc_payload)
    comparison = compare_stress_fields(gunw, gslc)
    map_space = map_space_stress(gslc_payload)

    xs, ys, terrain_transform = display_grid()
    terrain = load_terrain(terrain_transform, (len(ys), len(xs)))
    shade = hillshade(terrain)

    plt.rcParams.update(RC_PARAMS)
    stress_norm = TwoSlopeNorm(vmin=-STRESS_LIMIT_MPA, vcenter=0,
                               vmax=STRESS_LIMIT_MPA)
    difference_norm = TwoSlopeNorm(vmin=-DIFFERENCE_LIMIT_MPA, vcenter=0,
                                   vmax=DIFFERENCE_LIMIT_MPA)
    map_norm = TwoSlopeNorm(vmin=-MAP_STRESS_LIMIT_MPA, vcenter=0,
                            vmax=MAP_STRESS_LIMIT_MPA)
    lines = geometry_lines(map_space["geometry"])
    map_images = {
        key: blended_map(map_space["fields"][key], shade, cmc.vik, map_norm,
                         alpha=0.79)
        for key in ("sigma_n", "cfs")
    }
    along_span_km = float(gslc["along_edges_km"][-1] - gslc["along_edges_km"][0])
    depth_span_km = float(gslc["depth_edges_km"][-1] - gslc["depth_edges_km"][0])
    exaggeration = MAP_BOX_ASPECT * along_span_km / depth_span_km

    fig = plt.figure(figsize=(7.5, 7.3))
    grid = fig.add_gridspec(
        3, 3,
        width_ratios=[1, 1, 1.08],
        height_ratios=[1.5, 1, 1],
        left=0.072,
        right=0.96,
        bottom=0.15,
        top=0.895,
        hspace=0.38,
        wspace=0.20,
    )
    axes = [[fig.add_subplot(grid[row, column]) for column in range(3)]
            for row in range(3)]

    def draw_stress_map(ax, image, title, annotation, letter, *, first,
                        right_ticks=False):
        ax.imshow(
            image,
            extent=[X_MIN, X_MAX, Y_MIN, Y_MAX],
            origin="upper",
            interpolation="nearest",
            rasterized=True,
            zorder=1,
        )
        draw_geometry(ax, lines, legend=first)
        if first:
            add_scale_and_north(ax)
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.xaxis.set_major_locator(MultipleLocator(10_000))
        ax.yaxis.set_major_locator(MultipleLocator(10_000))
        formatter = FuncFormatter(lambda value, _: f"{value / 1_000:.0f}")
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        ax.set_xlabel("UTM 37N easting (km)", fontsize=6.6, labelpad=1)
        if first:
            ax.set_ylabel("UTM 37N northing (km)", fontsize=6.5, labelpad=1)
        else:
            ax.set_yticklabels([])
        add_lonlat_axes(ax, top=True, right=right_ticks)
        ax.tick_params(labelsize=6.2, length=2.3, width=0.6, pad=1.1)
        ax.set_title(title, fontweight="bold", pad=11)
        ax.text(0.018, 0.972, f"({letter})", transform=ax.transAxes, ha="left",
                va="top", fontsize=7.2, fontweight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                      "pad": 1.2})
        if annotation:
            ax.text(0.977, 0.028, annotation, transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=6.0,
                    bbox={"facecolor": "white", "edgecolor": "none",
                          "alpha": 0.84, "pad": 1.2})
        for spine in ax.spines.values():
            spine.set_linewidth(0.65)

    def draw_receiver_map(ax, values, norm, title, annotation, panel_index,
                          *, bottom_row):
        ax.pcolormesh(
            gslc["along_edges_km"],
            gslc["depth_edges_km"],
            np.ma.masked_invalid(values),
            cmap=cmc.vik,
            norm=norm,
            shading="flat",
            rasterized=True,
        )
        ax.set_xlim(gslc["along_edges_km"][0], gslc["along_edges_km"][-1])
        ax.set_ylim(gslc["depth_edges_km"][-1], gslc["depth_edges_km"][0])
        ax.set_xticks([-10, 0, 10])
        ax.set_yticks([0, 1, 2, 3, 4])
        if bottom_row:
            ax.set_xlabel("along strike, NNW → SSE (km)", fontsize=6.6, labelpad=1)
        else:
            ax.set_xticklabels([])
        ax.set_ylabel("depth (km)", fontsize=6.5, labelpad=1)
        ax.tick_params(labelsize=6.2, length=2.3, width=0.6, pad=1.1)
        ax.set_title(title, fontweight="bold", pad=3)
        ax.set_box_aspect(MAP_BOX_ASPECT)
        ax.text(0.018, 0.965, f"({string.ascii_lowercase[panel_index]})",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.2,
                fontweight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                      "pad": 1.2})
        ax.text(0.975, 0.04, annotation, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=6.0,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84,
                      "pad": 1.2})
        for spine in ax.spines.values():
            spine.set_linewidth(0.65)

    draw_stress_map(axes[0][0], map_images["sigma_n"],
                    "Normal stress at 1.0 km", "", "a", first=True)
    draw_stress_map(axes[0][1], map_images["cfs"], "Coulomb stress at 1.0 km",
                    "GSLC inputs · depth 1.0 km", "b", first=False,
                    right_ticks=True)

    depth_ax = axes[0][2]
    colors = {"sigma_n": "#b33a33", "cfs": "#285d88"}
    for key, label in (("sigma_n", "normal"), ("cfs", "Coulomb")):
        depth_ax.plot(gslc["profiles"][key], gslc["profiles"]["depth_km"],
                      color=colors[key], lw=1.25, label=f"GSLC {label}")
        depth_ax.plot(gunw["profiles"][key], gunw["profiles"]["depth_km"],
                      color=colors[key], lw=0.95, ls=(0, (4, 2)),
                      label=f"GUNW {label}")
    depth_ax.fill_betweenx(
        gslc["profiles"]["depth_km"],
        gslc["profiles"]["sigma_p16"],
        gslc["profiles"]["sigma_p84"],
        color=colors["sigma_n"],
        alpha=0.11,
        linewidth=0,
    )
    depth_ax.axvline(0, color="#666666", lw=0.6, ls=":")
    depth_ax.set_ylim(4, 0)
    depth_ax.set_xlim(-20, 15)
    depth_ax.set_xlabel("along-strike median (MPa)", fontsize=6.6, labelpad=1)
    # The depth axis sits on the right so panel (b)'s outer °N labels have
    # room on this panel's left side.
    depth_ax.yaxis.tick_right()
    depth_ax.yaxis.set_label_position("right")
    depth_ax.set_ylabel("receiver depth (km)", fontsize=6.5, labelpad=2)
    depth_ax.set_title("Depth summary", fontweight="bold", pad=3)
    depth_ax.tick_params(labelsize=6.2, length=2.3, width=0.6, pad=1.1)
    depth_ax.grid(color="#a8a8a8", alpha=0.23, lw=0.45)
    depth_ax.legend(fontsize=6.0, frameon=False, loc="lower left", ncols=2,
                    handlelength=1.8, columnspacing=0.8)
    depth_ax.set_box_aspect(1.0)
    depth_ax.text(0.975, 0.965, "(c)", transform=depth_ax.transAxes,
                  ha="right", va="top", fontsize=7.2, fontweight="bold")

    sigma_metrics = comparison["metrics"]["sigma_n"]
    cfs_metrics = comparison["metrics"]["cfs"]
    draw_receiver_map(
        axes[1][0], gslc["maps"]["sigma_n"], stress_norm, "GSLC normal stress",
        f"median {np.nanmedian(gslc['maps']['sigma_n']):+.2f} MPa", 3,
        bottom_row=False,
    )
    draw_receiver_map(
        axes[1][1], gslc["maps"]["cfs"], stress_norm, "GSLC Coulomb stress",
        f"median {np.nanmedian(gslc['maps']['cfs']):+.2f} MPa", 4,
        bottom_row=False,
    )

    dip_ax = axes[1][2]
    dips = np.array(RECEIVER_DIPS, dtype=float)
    for key, label in (("median_sigma_n_mpa", "normal"),
                       ("median_cfs_mpa", "Coulomb")):
        color = colors["sigma_n" if "sigma" in key else "cfs"]
        for fields, chain, style in ((gslc, "GSLC", "-"),
                                     (gunw, "GUNW", (0, (4, 2)))):
            values = [fields["dip_summary"][f"{dip:.0f}"][key] for dip in dips]
            dip_ax.plot(
                dips, values,
                color=color,
                lw=1.15 if chain == "GSLC" else 0.9,
                ls=style,
                marker="o" if chain == "GSLC" else None,
                ms=2.7,
                label=f"{chain} {label}",
            )
    dip_ax.axhline(0, color="#666666", lw=0.6, ls=":")
    dip_ax.set_xlim(52, 78)
    dip_ax.set_xticks(dips)
    dip_ax.set_xlabel("receiver dip (degrees)", fontsize=6.6, labelpad=1)
    dip_ax.set_ylabel("overall median (MPa)", fontsize=6.5, labelpad=1)
    dip_ax.set_title("Receiver-dip sensitivity", fontweight="bold", pad=3)
    dip_ax.tick_params(labelsize=6.2, length=2.3, width=0.6, pad=1.1)
    dip_ax.grid(color="#a8a8a8", alpha=0.23, lw=0.45)
    dip_ax.legend(fontsize=6.0, frameon=False, loc="best", ncols=2,
                  handlelength=1.7, columnspacing=0.7)
    dip_ax.set_box_aspect(MAP_BOX_ASPECT)
    dip_ax.text(0.975, 0.965, "(f)", transform=dip_ax.transAxes, ha="right",
                va="top", fontsize=7.2, fontweight="bold")

    draw_receiver_map(
        axes[2][0], comparison["differences"]["sigma_n"], difference_norm,
        "GSLC − GUNW normal",
        f"relative L2 {100 * sigma_metrics['relative_l2']:.1f}% · signs "
        f"{100 * sigma_metrics['sign_agreement_fraction']:.1f}%",
        6, bottom_row=True,
    )
    draw_receiver_map(
        axes[2][1], comparison["differences"]["cfs"], difference_norm,
        "GSLC − GUNW Coulomb",
        f"relative L2 {100 * cfs_metrics['relative_l2']:.1f}% · signs "
        f"{100 * cfs_metrics['sign_agreement_fraction']:.1f}%",
        7, bottom_row=True,
    )
    notes_ax = axes[2][2]
    notes_ax.axis("off")
    notes_ax.text(
        0.0,
        0.98,
        "\n".join(
            [
                "+ = unclamping / normal-slip CFS",
                "map receivers: strike "
                f"{map_space['geometry']['strike']:.1f}°, dip {MAIN_DIP:.0f}° "
                "toward the dike,",
                f"normal rake, friction {FRICTION:.1f}, depth "
                f"{MAP_DEPTH_M / 1_000:.1f} km",
                "dike-plane strip (±"
                f"{DIKE_STANDOFF_M / 1_000:.1f} km) masked",
                f"receiver-plane depth exaggerated ×{exaggeration:.1f}",
            ]
        ),
        transform=notes_ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.0,
        color="#3a3a3a",
        linespacing=1.55,
    )

    def add_colorbar(position, norm, label, ticks, ticklabels):
        cax = fig.add_axes(position)
        bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmc.vik), cax=cax,
                           orientation="horizontal")
        bar.set_ticks(ticks)
        bar.set_ticklabels(ticklabels)
        bar.set_label(label, fontsize=6.7, labelpad=1.5)
        bar.ax.tick_params(labelsize=6.1, length=2, pad=1.1)

    add_colorbar([0.075, 0.072, 0.235, 0.013], map_norm,
                 "map stress change (MPa; clipped)", [-5, 0, 5],
                 ["−5", "0", "+5"])
    add_colorbar([0.385, 0.072, 0.235, 0.013], stress_norm,
                 "receiver-plane stress (MPa; clipped)", [-20, 0, 20],
                 ["−20", "0", "+20"])
    add_colorbar([0.695, 0.072, 0.20, 0.013], difference_norm,
                 "GSLC − GUNW (MPa; clipped)", [-2, 0, 2], ["−2", "0", "+2"])

    fig.suptitle(
        "Manda Hararo 2026 — dike-induced stress on the east-fault receivers",
        x=0.5,
        y=0.992,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.955,
        "map stress at 1.0 km depth · trace-supported east-fault receivers, "
        "65° main case · source: dike opening only · identical GSLC/GUNW "
        "receivers",
        ha="center",
        va="bottom",
        fontsize=6.7,
        color="#3a3a3a",
    )
    fig.text(
        0.5,
        0.008,
        "All panels resolve stress on the trace-supported east-fault "
        "orientation with the parameters listed above.\n"
        "The stress source is the saved distributed dike opening only; the "
        "fitted fault slip is excluded.",
        ha="center",
        va="bottom",
        fontsize=6.0,
        color="#4a4a4a",
        linespacing=1.45,
    )

    png = RESULTS / "summary_fault_stress.png"
    pdf = RESULTS / "summary_fault_stress.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


# ---------------------------------------------------------------------------


BUILDERS = {
    "los": build_los_figure,
    "decomposition": build_decomposition_figure,
    "source": build_source_figure,
    "gradient": build_gradient_figure,
    "stress": build_stress_figure,
    "dike_frame": build_dike_frame_figure,
    "wright_los": build_wright_los_figure,
    "wright_uperp": build_wright_uperp_figure,
}

RESULT_INPUTS = (
    "distributed_fit_mapped.json",
    "distributed_fit_mapped_gslc.json",
    "distributed_fit_geometry_4a7b2386c8b5_gslc.json",
    "decomposition_joint_gslc.json",
    "fault_traces.json",
    "okada_fit.json",
    "s1_tracks.json",
)


def write_manifest(outputs: list[Path]) -> Path:
    manifest = {
        "scope": (
            "Identity record for the summary figure set. Raw data files are "
            "covered by the archive checksum manifests in data/."
        ),
        "generated_by": "mh26.summary_figures",
        "results_inputs": {
            name: sha256_file(RESULTS / name) for name in RESULT_INPUTS
        },
        "terrain": {
            "path": str(DEM_PATH.relative_to(REPO)),
            "sha256": DEM_SHA256,
            "url": DEM_URL,
        },
        "outputs": {
            path.name: sha256_file(path) for path in sorted(outputs)
        },
    }
    path = RESULTS / "summary_figures.json"
    path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=sorted(BUILDERS))
    args = parser.parse_args()

    require_dem()
    names = [args.only] if args.only else list(BUILDERS)
    outputs: list[Path] = []
    for name in names:
        print(f"building {name} ...")
        outputs.extend(BUILDERS[name]())
        for path in outputs[-2:]:
            print(f"figure: {path}")
    if not args.only:
        manifest = write_manifest(outputs)
        print(f"manifest: {manifest}")


if __name__ == "__main__":
    main()
