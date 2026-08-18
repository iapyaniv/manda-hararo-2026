"""Okada rectangular tensile dislocation (dike) in an elastic half-space.

Wraps pyrocko's okada_ext with this project's parameterisation:
    params = [east_c, north_c, strike_deg, length, top, bottom, opening]
(metres and degrees; top/bottom are depths of the dike's upper and lower
edges; dip is fixed at 90).

Run `python -m mh26.okada` to verify the wrapped sign conventions:
an inflating horizontal sill must uplift the surface, and a vertical
opening dike must push its two sides apart.
"""

from __future__ import annotations

import numpy as np
from pyrocko.modelling import okada_ext

from .params import LAME_LAMBDA, SHEAR_MODULUS


def forward(params: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Surface displacement rows (north, east, up), metres, at pts=(north, east)."""
    e_c, n_c, strike, length, top, bottom, opening = params
    if bottom <= top:
        raise ValueError("bottom depth must be greater than top depth")
    depth_c = 0.5 * (top + bottom)
    hw = 0.5 * (bottom - top)
    patch = np.array([[n_c, e_c, depth_c, strike, 90.0,
                       -length / 2, length / 2, -hw, hw]])
    disl = np.array([[0.0, 0.0, opening]])
    rec = np.zeros((len(pts), 3))
    rec[:, 0] = pts[:, 0]  # north
    rec[:, 1] = pts[:, 1]  # east
    res = okada_ext.okada(patch, disl, rec, LAME_LAMBDA, SHEAR_MODULUS,
                          nthreads=0, rotate_sdn=0, stack_sources=1)
    un, ue, ud = res[:, 0], res[:, 1], res[:, 2]
    return np.column_stack([un, ue, -ud])  # convert down -> up


def selftest() -> bool:
    """Check the sign conventions of the wrapped Okada implementation."""
    disl = np.array([[0.0, 0.0, 1.0]])
    # Inflating horizontal sheet (dip 0) at 2 km depth must UPLIFT the surface above.
    patch = np.array([[0.0, 0.0, 2000.0, 0.0, 0.0, -3000, 3000, -3000, 3000]])
    rec = np.array([[0.0, 0.0, 0.0]])
    res = okada_ext.okada(patch, disl, rec, LAME_LAMBDA, SHEAR_MODULUS,
                          nthreads=0, rotate_sdn=0, stack_sources=1)
    up = -res[0, 2]
    ok_sill = up > 0
    print(f"selftest: surface motion above inflating sill = {up:+.3f} m "
          f"({'PASS (uplift)' if ok_sill else 'FAIL - sign convention wrong'})")
    # Vertical opening dike: stations on both sides must move APART.
    patch = np.array([[0.0, 0.0, 3000.0, 0.0, 90.0, -5000, 5000, -2500, 2500]])
    rec = np.array([[0.0, 3000.0, 0.0], [0.0, -3000.0, 0.0]])
    res = okada_ext.okada(patch, disl, rec, LAME_LAMBDA, SHEAR_MODULUS,
                          nthreads=0, rotate_sdn=0, stack_sources=1)
    div = res[0, 1] - res[1, 1]
    ok_dike = div > 0
    print(f"selftest: cross-dike divergence = {div:+.3f} m "
          f"({'PASS' if ok_dike else 'FAIL'})")
    return ok_sill and ok_dike


if __name__ == "__main__":
    raise SystemExit(0 if selftest() else 1)
