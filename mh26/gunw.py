"""Loading and correcting the NISAR GUNW interferogram crops (data/nisar/*.npz).

Each .npz holds a 90 x 90 km crop of one L2 GUNW product: unwrappedPhase,
coherence, ionospherePhaseScreen, connectedComponents, mask, the x/y UTM
coordinate vectors, the EPSG code, a metadata summary string, and the
radar-grid geometry cubes (incidence angle, LOS unit vectors).

Conventions: LOS displacement d = -phi * lambda / (4 pi), positive = motion
TOWARD the satellite (range decrease). x ascends east; y descends, so row 0
is the northern edge.
"""

from __future__ import annotations

import re

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .params import (
    CENTER_E,
    CENTER_N,
    COHERENCE_MIN,
    BLOCK_PIXELS,
    DATA_NISAR,
    FARFIELD_RADIUS,
    LOS_HEIGHT_AE,
    WAVELENGTH_L,
)


CONNECTED_COMPONENT_FILL = np.iinfo(np.uint16).max
VALID_SAMPLE_FILL = 255


def load_pair(tag: str) -> dict | None:
    """Load data/nisar/<tag>.npz as a dict of arrays, or None if absent."""
    path = DATA_NISAR / f"{tag}.npz"
    if not path.exists():
        return None
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def pair_window(d: dict) -> str | None:
    """'YYYY-MM-DD→YYYY-MM-DD' reference/secondary window from the metadata."""
    dates = re.findall(r"(\d{4}-\d{2}-\d{2})T", str(d.get("meta", "")))
    return f"{dates[0]}→{dates[1]}" if len(dates) >= 2 else None


def valid_sample_mask(mask: np.ndarray) -> np.ndarray:
    """Return pixels valid in both the reference and secondary RSLCs.

    The GUNW product specification encodes water/reference-subswath/
    secondary-subswath as a three-digit value.  A zero in either of the last
    two digits means that the corresponding RSLC sample is invalid; 255 is the
    fill value.  Current provisional products can carry additional flags above
    the low byte, so decode the documented value from that byte.

    Applying this mask is important for the provisional products: otherwise
    exposed subswath-edge samples can pass both coherence and connected-
    component checks and appear as narrow displacement stripes.
    """
    code = np.asarray(mask, dtype=np.uint32) & 0xFF
    reference_subswath = (code // 10) % 10
    secondary_subswath = code % 10
    return (
        (code != VALID_SAMPLE_FILL)
        & (reference_subswath > 0)
        & (secondary_subswath > 0)
    )


def corrected_los(d: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return analysis-screened LOS displacement (m) and raw coherence.

    The returned LOS field is ionosphere-corrected and keeps only finite pixels
    with coherence above the project threshold, a real unwrapping connected
    component, and valid samples in both contributing RSLCs.
    """
    required = (
        "unwrappedPhase",
        "coherence",
        "ionospherePhaseScreen",
        "connectedComponents",
        "mask",
    )
    missing = [name for name in required if name not in d]
    if missing:
        raise KeyError(f"required GUNW layers missing: {', '.join(missing)}")

    phi = d["unwrappedPhase"].astype(np.float64)
    ionosphere = np.asarray(d["ionospherePhaseScreen"], dtype=np.float64)
    if not np.isfinite(ionosphere).any():
        raise ValueError("ionospherePhaseScreen contains no finite samples")
    phi = phi - ionosphere

    coh = np.asarray(d["coherence"])
    cc = np.asarray(d["connectedComponents"])
    shapes = {np.asarray(d[name]).shape for name in required}
    if len(shapes) != 1:
        raise ValueError(f"GUNW analysis layers have inconsistent shapes: {sorted(shapes)}")
    good = (
        np.isfinite(phi)
        & np.isfinite(coh)
        & (coh > COHERENCE_MIN)
        & (cc > 0)
        & (cc != CONNECTED_COMPONENT_FILL)
        & valid_sample_mask(d["mask"])
    )
    phi = np.where(good, phi, np.nan)
    los = -phi * WAVELENGTH_L / (4 * np.pi)
    return los, coh


def reference_farfield(los: np.ndarray, x: np.ndarray, y: np.ndarray,
                       cx: float = CENTER_E, cy: float = CENTER_N,
                       r_in: float = FARFIELD_RADIUS) -> np.ndarray:
    """Remove a plane (per-pair offset + ramp) fitted outside radius r_in of (cx, cy)."""
    if los.shape != (len(y), len(x)):
        raise ValueError("LOS grid shape does not match x/y coordinate vectors")
    X, Y = np.meshgrid(x, y)
    far = (np.hypot(X - cx, Y - cy) > r_in) & np.isfinite(los)
    if far.sum() < 500:
        return los - np.nanmedian(los)
    A = np.c_[np.ones(far.sum()), (X[far] - cx) / 1e4, (Y[far] - cy) / 1e4]
    coef, *_ = np.linalg.lstsq(A, los[far], rcond=None)
    plane = coef[0] + coef[1] * (X - cx) / 1e4 + coef[2] * (Y - cy) / 1e4
    return los - plane


def remove_farfield_plane(pts: np.ndarray, values: np.ndarray,
                          r_in: float = FARFIELD_RADIUS,
                          cx: float = CENTER_E, cy: float = CENTER_N):
    """Apply the observation reference operator to model values.

    ``reference_farfield`` removes a plane fitted on the far field from every
    observation, and part of the true deformation signal goes with it. A
    model compared against those observations must lose its own far-field
    plane through the same operator; otherwise the missing plane is absorbed
    into the source parameters and biases the source low. Accepts a vector
    or column-stacked values (e.g. Green-function columns). Returns the
    referenced values and the plane coefficients [constant, east, north].
    """
    radius = np.hypot(pts[:, 1] - cx, pts[:, 0] - cy)
    far = radius > r_in
    if int(far.sum()) < 3:
        raise ValueError("too few far-field samples for the reference plane")
    basis = np.column_stack([np.ones(len(pts)), pts[:, 1] - cx,
                             pts[:, 0] - cy])
    coef, *_ = np.linalg.lstsq(basis[far], values[far], rcond=None)
    return values - basis @ coef, coef


def block_average(arr: np.ndarray, x: np.ndarray, y: np.ndarray,
                  block: int = BLOCK_PIXELS):
    """Block-average a grid; blocks under 30% valid become NaN.

    Returns (pts, vals, (xs, ys, grid)): pts as (north, east) rows for the
    valid blocks, their values, and the full block grid with coordinates.
    """
    if block <= 0:
        raise ValueError("block must be a positive integer")
    ny, nx = arr.shape
    by, bx = ny // block, nx // block
    if by == 0 or bx == 0:
        raise ValueError("block is larger than the input grid")
    a = arr[: by * block, : bx * block].reshape(by, block, bx, block)
    n_ok = np.isfinite(a).sum(axis=(1, 3))
    sums = np.nansum(a, axis=(1, 3))
    m = np.full(sums.shape, np.nan, dtype=np.float64)
    np.divide(sums, n_ok, out=m, where=n_ok > 0)
    m[n_ok < block * block * 0.3] = np.nan
    xs = x[: bx * block].reshape(bx, block).mean(axis=1)
    ys = y[: by * block].reshape(by, block).mean(axis=1)
    X, Y = np.meshgrid(xs, ys)
    ok = np.isfinite(m)
    return np.column_stack([Y[ok], X[ok]]), m[ok], (xs, ys, m)


def los_unit_neu(d: dict, pts: np.ndarray,
                 height: float = LOS_HEIGHT_AE) -> np.ndarray:
    """Local LOS unit vectors (north, east, up) at ``pts=(north, east)``.

    East and north are trilinearly interpolated from the product geometry
    cubes at ``height`` metres above the ellipsoid. Raises if the required
    cubes are missing or invalid; a nominal geometry would silently change
    the science.
    """
    required = ("cube_heightAboveEllipsoid", "cube_yCoordinates",
                "cube_xCoordinates", "cube_losUnitVectorX",
                "cube_losUnitVectorY")
    if any(name not in d for name in required):
        raise KeyError("LOS unit-vector cubes missing from this crop; "
                       "re-extract with the geometry cubes included")
    pts = np.asarray(pts, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("pts must have shape (n, 2) as (north, east)")

    h = np.asarray(d["cube_heightAboveEllipsoid"], dtype=float)
    y = np.asarray(d["cube_yCoordinates"], dtype=float)
    x = np.asarray(d["cube_xCoordinates"], dtype=float)
    ho, yo, xo = np.argsort(h), np.argsort(y), np.argsort(x)
    query = np.column_stack([np.full(len(pts), height), pts[:, 0], pts[:, 1]])
    components = []
    for name in ("cube_losUnitVectorX", "cube_losUnitVectorY"):
        cube = np.asarray(d[name], dtype=float)[np.ix_(ho, yo, xo)]
        interpolator = RegularGridInterpolator(
            (h[ho], y[yo], x[xo]), cube,
            bounds_error=True,
        )
        components.append(interpolator(query))
    ue, un = components
    if not (np.isfinite(ue).all() and np.isfinite(un).all()):
        raise ValueError("LOS unit-vector interpolation produced non-finite values")
    horizontal_norm_sq = ue**2 + un**2
    if np.any(horizontal_norm_sq > 1.0 + 1e-5):
        raise ValueError("LOS horizontal components have norm greater than one")
    uu = np.sqrt(np.clip(1.0 - horizontal_norm_sq, 0.0, None))
    return np.column_stack([un, ue, uu])
