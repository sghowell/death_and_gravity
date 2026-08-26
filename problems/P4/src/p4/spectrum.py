"""Zero counting and refinement for the matching function E(kappa) of ``p4.perturb``.

* ``winding_number``: argument-principle count Z - P on a rectangle by adaptive
  sampling of E on the boundary (phase increments kept below pi/2).
* ``locate_zeros``: recursive bisection of a rectangle down to small boxes.
* ``refine_zero``: complex secant iteration started from a box centre.
E is analytic in {Re kappa > -1.1} (the sonic-point resonances sit at
kappa = -0.099 - 1.099 n, n = 1, 2, ...), so Z - P = Z there.
"""
from __future__ import annotations

import numpy as np


def _rect_boundary(rect, n_per_side):
    x0, x1, y0, y1 = rect
    s = np.linspace(0, 1, n_per_side, endpoint=False)
    bottom = x0 + (x1 - x0) * s + 1j * y0
    right = x1 + 1j * (y0 + (y1 - y0) * s)
    top = x1 - (x1 - x0) * s + 1j * y1
    left = x0 + 1j * (y1 - (y1 - y0) * s)
    return np.concatenate([bottom, right, top, left])


def winding_number(prob, rect, n_per_side=40, max_rounds=6, max_dphase=np.pi / 2, verbose=False):
    """Winding number of E along the (closed) boundary of rect = (re0, re1, im0, im1).
    Returns (winding, kappas, E_values)."""
    ks = _rect_boundary(rect, n_per_side)
    Es = prob.E(ks)
    for _ in range(max_rounds):
        ph = np.angle(np.roll(Es, -1) / Es)
        bad = np.where(np.abs(ph) > max_dphase)[0]
        if len(bad) == 0:
            break
        new = 0.5 * (ks[bad] + np.roll(ks, -1)[bad])
        Enew = prob.E(new)
        order = np.argsort(np.concatenate([np.arange(len(ks)), bad + 0.5]))
        ks = np.concatenate([ks, new])[order]
        Es = np.concatenate([Es, Enew])[order]
        if verbose:
            print(f"  refined {len(bad)} segments -> {len(ks)} points")
    ph = np.angle(np.roll(Es, -1) / Es)
    if np.any(np.abs(ph) > max_dphase):
        raise RuntimeError("phase increments not resolved on the contour")
    return int(np.rint(ph.sum() / (2 * np.pi))), ks, Es


def locate_zeros(prob, rect, n_per_side=24, min_size=0.05, verbose=False):
    """Recursively bisect rect until each box with nonzero winding has size < min_size.
    Returns list of (box, winding)."""
    w, _, _ = winding_number(prob, rect, n_per_side)
    if verbose:
        print(f"box {rect}: winding {w}")
    if w == 0:
        return []
    x0, x1, y0, y1 = rect
    if max(x1 - x0, y1 - y0) < min_size:
        return [(rect, w)]
    # split (slightly off-centre if the split line runs through a zero, e.g. Im = 0)
    for frac in (0.5, 0.46, 0.54, 0.42, 0.58):
        if x1 - x0 >= y1 - y0:
            xm = x0 + frac * (x1 - x0)
            subs = [(x0, xm, y0, y1), (xm, x1, y0, y1)]
        else:
            ym = y0 + frac * (y1 - y0)
            subs = [(x0, x1, y0, ym), (x0, x1, ym, y1)]
        try:
            out = []
            for sb in subs:
                out += locate_zeros(prob, sb, n_per_side, min_size, verbose)
            return out
        except RuntimeError:
            continue
    raise RuntimeError(f"could not split {rect} away from zeros of E")


def refine_zero(prob, k0, k1=None, tol=1e-12, maxiter=40, real=False):
    """Complex secant iteration for a zero of E near k0; returns (kappa, |E|, iterations)."""
    k0 = complex(k0)
    k1 = k0 + 1e-3 if k1 is None else complex(k1)
    if real:
        k0, k1 = k0.real, k1.real
    E0, E1 = [complex(e) for e in prob.E(np.array([k0, k1]))]
    if real:
        E0, E1 = E0.real, E1.real
    for it in range(maxiter):
        if E1 == E0:
            break
        k2 = k1 - E1 * (k1 - k0) / (E1 - E0)
        E2 = complex(prob.E(np.array([k2]))[0])
        if real:
            E2 = E2.real
        k0, E0, k1, E1 = k1, E1, k2, E2
        if abs(k1 - k0) < tol * max(1.0, abs(k1)):
            break
    return k1, abs(E1), it + 1


def real_axis_scan(prob, k_lo, k_hi, n=41):
    """E on a real grid; returns (kappas, E.real, brackets with sign change)."""
    ks = np.linspace(k_lo, k_hi, n)
    E = prob.E(ks).real
    br = [(ks[i], ks[i + 1]) for i in range(n - 1) if np.sign(E[i]) != np.sign(E[i + 1])]
    return ks, E, br
