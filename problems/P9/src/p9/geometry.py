"""Exact linear geometry for piecewise-linear u(x), x = ln(1+z), on a general grid.

u(x) = c/(r_d H(z)) is linear on each grid segment [x_k, x_{k+1}]. Then
  D_H(z)/r_d = u(x(z))                                (linear interpolation)
  D_M(z)/r_d = int_0^{x(z)} u(x') e^{x'} dx'          (closed form, linear in nodes)
All coefficients are nonnegative.
"""

from __future__ import annotations

import numpy as np


def grid(N: int, z_max: float = 2.5) -> np.ndarray:
    """Uniform grid in x with N segments (kept for tests and the v1 class)."""
    return np.linspace(0.0, np.log1p(z_max), N + 1)


def geometric_grid(z_max: float = 2.5, z_first: float = 0.005, ratio: float = 1.1,
                   h_max: float = 0.02) -> np.ndarray:
    """0, then geometric nodes from ln(1+z_first) with the given ratio until the spacing
    reaches h_max, then uniform spacing h_max up to X = ln(1+z_max) (last segment may be shorter)."""
    X = np.log1p(z_max)
    nodes = [0.0, np.log1p(z_first)]
    while True:
        x = nodes[-1]
        h = x * (ratio - 1.0)
        if h >= h_max:
            break
        nodes.append(x * ratio)
    x = nodes[-1]
    n_uni = int(np.ceil((X - x) / h_max))
    nodes += list(np.linspace(x, X, n_uni + 1)[1:])
    return np.array(nodes)


def segment_index(x_nodes: np.ndarray, x: np.ndarray) -> np.ndarray:
    """k with x_nodes[k] <= x <= x_nodes[k+1] (right endpoint assigned to the last segment)."""
    N = len(x_nodes) - 1
    k = np.searchsorted(x_nodes, x, side="right") - 1
    return np.clip(k, 0, N - 1)


def dm_matrix(x_nodes: np.ndarray, z) -> np.ndarray:
    """A with (A @ u)[r] = D_M(z_r)/r_d for node values u (piecewise linear in x)."""
    x = np.log1p(np.asarray(z, dtype=float))
    if np.any(x < 0) or np.any(x > x_nodes[-1] + 1e-12):
        raise ValueError("redshift outside grid")
    N = len(x_nodes) - 1
    hs = np.diff(x_nodes)
    ex = np.exp(x_nodes)
    # full-segment contributions to (u_k, u_{k+1})
    full_lo = (ex[1:] - (hs + 1.0) * ex[:-1]) / hs
    full_hi = ((hs - 1.0) * ex[1:] + ex[:-1]) / hs
    A = np.zeros((len(x), N + 1))
    k = segment_index(x_nodes, x)
    for r in range(len(x)):
        kk = k[r]
        A[r, :kk] += full_lo[:kk]
        A[r, 1 : kk + 1] += full_hi[:kk]
        xa, xb, h = x_nodes[kk], x_nodes[kk + 1], hs[kk]
        ea, exv = ex[kk], np.exp(x[r])
        A[r, kk] += ((xb - x[r] + 1.0) * exv - (h + 1.0) * ea) / h
        A[r, kk + 1] += ((x[r] - xa - 1.0) * exv + ea) / h
    return A


def dh_matrix(x_nodes: np.ndarray, z) -> np.ndarray:
    """A with (A @ u)[r] = D_H(z_r)/r_d = u(x(z_r)) by linear interpolation."""
    x = np.log1p(np.asarray(z, dtype=float))
    N = len(x_nodes) - 1
    hs = np.diff(x_nodes)
    k = segment_index(x_nodes, x)
    A = np.zeros((len(x), N + 1))
    t = (x - x_nodes[k]) / hs[k]
    for r in range(len(x)):
        A[r, k[r]] += 1.0 - t[r]
        A[r, k[r] + 1] += t[r]
    return A


def segment_increment_factor(x_nodes: np.ndarray, z) -> tuple[np.ndarray, np.ndarray]:
    """For each z: segment index k and (e^{x} - e^{x_k}); D_M(z) - D_M(z_k) lies between
    min(u_k,u_{k+1}) and max(u_k,u_{k+1}) times this factor."""
    x = np.log1p(np.asarray(z, dtype=float))
    k = segment_index(x_nodes, x)
    return k, np.exp(x) - np.exp(x_nodes[k])


def lcdm_u_nodes(x_nodes: np.ndarray, omega_m: float, h_rd: float) -> np.ndarray:
    """Node values u_i = c/(r_d H(z_i)) for flat LCDM with h*r_d = h_rd [Mpc]."""
    from . import C_KM_S

    z = np.expm1(x_nodes)
    E = np.sqrt(omega_m * (1 + z) ** 3 + (1 - omega_m))
    return C_KM_S / (100.0 * h_rd * E)


# ---------------------------------------------------------------------------
# Interpolation-error bound for kappa(x) = log10[ D_M(x) / (e^x - 1) ]  (FORMULATION v2, §3.2)
# ---------------------------------------------------------------------------

def _I_pm(x: float, L: float, sign: int) -> float:
    """int_0^x e^{sign*L*(x-x')} e^{x'} dx'  (closed form)."""
    a = 1.0 - sign * L
    if abs(a) < 1e-12:
        return np.exp(sign * L * x) * x
    return np.exp(sign * L * x) * (np.exp(a * x) - 1.0) / a


def theta_bounds(x: float, L: float) -> tuple[float, float]:
    """Class-only bounds on theta(x) = u(x)(e^x-1)/D_M(x) from |d ln u/dx| <= L."""
    em1 = np.expm1(x)
    return em1 / _I_pm(x, L, +1), em1 / _I_pm(x, L, -1)


def kappa_second_derivative_bound(xa: float, xb: float, L: float) -> float:
    """Rigorous bound on sup_{x in [xa,xb]} |F''(x)|, F = ln[D_M/(e^x-1)], over the class.

    F'' = q [ theta a - (theta-1)(q(1+theta) - 1) ],  q = e^x/(e^x-1), a = u'/u in [-L, L].
    q is decreasing in x; the theta-range widens with x. Each monotone factor is evaluated at
    its worst endpoint. For xa = 0 the first segment is handled by the exact formula
    F = ln(u_0 + b*phi(x)), |b|/u <= L, phi'' <= 1/6, phi' <= 1/2 on [0, xb] (xb <= 0.1 assumed).
    """
    if xa <= 0.0:
        if xb > 0.1:
            raise ValueError("first segment too long for the closed-form bound")
        return L * (1.0 / 6.0) * 1.02 + L * L * 0.25 * 1.02
    q = np.exp(xa) / np.expm1(xa)
    th_lo, th_hi = theta_bounds(xb, L)
    dth = max(abs(th_lo - 1.0), abs(th_hi - 1.0))
    return q * (L * th_hi + dth * (q * (1.0 + th_hi) - 1.0))


def kappa_interp_slack(x_nodes: np.ndarray, L: float) -> np.ndarray:
    """Per-segment bound e_k on |kappa(x) - linear interpolation of kappa between nodes| (log10 units)."""
    hs = np.diff(x_nodes)
    out = np.empty(len(hs))
    for k in range(len(hs)):
        B = kappa_second_derivative_bound(x_nodes[k], x_nodes[k + 1], L)
        out[k] = hs[k] ** 2 / 8.0 * B / np.log(10.0)
    return out
