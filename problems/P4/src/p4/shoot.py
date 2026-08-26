"""Stage B: shooting on the sonic-point velocity V0 for regular-centre CSS solutions.

For each V0 and analytic branch we start at x = -delta from the sonic-point
Taylor series, integrate the constraint-reduced, centre-scaled system
(N, W, V) toward x -> -inf, and read off
    F(V0) = m~(x_end) = (1 - 1/A) e^x |_{x_end}       (central-mass mismatch),
which tends to the coefficient of the irregular mode  A - 1 ~ m~ e^{-x}
(a point mass at the centre; equivalently NV + 1/2 ~ e^{-3x}); regular
centre  <=>  F = 0.  Since m~ is monotone in x (dm~/dx > 0) its sign at an
early stop is still meaningful for bracketing.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import css, taylor

SQ3 = css.SQ3


@dataclass
class Shot:
    V0: float
    branch: int
    x_stop: float
    reason: str          # 'end' | 'A_blowup' | 'V_lightspeed' | 'sonic' | 'fail'
    mass: float          # m~ at x_stop
    y_end: np.ndarray    # scaled variables (Ah, Nh, Wh, Vh) at x_stop
    n_zeros_V: int       # sign changes of V on (x_stop, -delta]
    delta: float
    sol: object = None


def _events():
    def ev_A(x, yh3):           # 2m/r -> 1 (A -> infinity)
        return 0.999 - css.scaled3_to_full(x, yh3)[0] * np.exp(2 * x)
    ev_A.terminal = True

    def ev_V(x, yh3):           # |V| -> 1
        return 0.9995 - abs(yh3[2] * np.exp(x))
    ev_V.terminal = True

    def ev_sonic(x, yh3):       # another sonic point: |Delta| -> 0 (relative), from either side
        y = css.from_scaled(x, css.scaled3_to_full(x, yh3))
        q = css.coeffs(*y)
        D, _, _ = css.det_and_numerators(q)
        return abs(D) / (abs(q["a"] * q["d"]) + abs(q["b"] * q["c"])) - 1e-7
    ev_sonic.terminal = True
    return [ev_A, ev_V, ev_sonic]


def start_data(V0, branch, K=40, delta=None):
    """Series start point: returns (delta, y_plain(-delta), coef)."""
    coef = taylor.background_series(V0, branch, K=K)
    if delta is None:
        R = taylor.radius_estimate(coef)
        if not np.isfinite(R) or R < 2e-2:       # ill-conditioned series: linear data, tiny step
            coef = coef[:, :2]
            delta = 1e-5
        else:
            delta = float(np.clip(R / 6, 1e-4, 0.1))
    return delta, taylor.eval_series(coef, -delta), coef


def shoot(V0, branch, x_end=-15.0, K=40, delta=None, rtol=1e-13, atol=1e-16,
          keep_sol=False, method="DOP853"):
    """Integrate from the sonic point toward the centre; see module docstring."""
    delta, y_start, _ = start_data(V0, branch, K, delta)
    yh0 = css.to_scaled(-delta, y_start)[1:]
    evs = _events()
    try:
        sol = solve_ivp(css.rhs_scaled3, (-delta, x_end), yh0, method=method, rtol=rtol,
                        atol=atol, events=evs, dense_output=keep_sol)
    except Exception:                        # noqa: BLE001
        return Shot(V0, branch, -delta, "fail", np.nan, css.scaled3_to_full(-delta, yh0), 0, delta)
    reasons = ["A_blowup", "V_lightspeed", "sonic"]
    reason = "end" if sol.status == 0 else ("fail" if sol.status < 0 else
                                            reasons[[len(t) > 0 for t in sol.t_events].index(True)])
    x_stop = sol.t[-1]
    yh = css.scaled3_to_full(x_stop, sol.y[:, -1])
    V = sol.y[2] * np.exp(sol.t)
    nz = int(np.sum(np.sign(V[1:]) != np.sign(V[:-1])))
    return Shot(V0, branch, x_stop, reason, css.mass_function(x_stop, yh), yh, nz, delta,
                sol if keep_sol else None)


def mismatch(V0, branch, **kw):
    return shoot(V0, branch, **kw).mass


def scan(V0_grid, branch, x_end=-8.0, rtol=1e-10, atol=1e-13):
    """Coarse scan; returns list of Shot (cheap tolerances, early x_end)."""
    return [shoot(v, branch, x_end=x_end, rtol=rtol, atol=atol) for v in V0_grid]


def brackets_from_scan(shots):
    """Consecutive grid points where the mismatch changes sign."""
    out = []
    for s0, s1 in zip(shots[:-1], shots[1:]):
        if np.isfinite(s0.mass) and np.isfinite(s1.mass) and np.sign(s0.mass) != np.sign(s1.mass):
            out.append((s0.V0, s1.V0))
    return out


def refine_root(a, b, branch, x_end=-15.0, xtol=1e-14, **kw):
    """brentq on the mismatch; returns (V0, Shot at the root)."""
    f = lambda v: mismatch(v, branch, x_end=x_end, **kw)      # noqa: E731
    V0 = brentq(f, a, b, xtol=xtol, rtol=4 * np.finfo(float).eps, maxiter=200)
    return V0, shoot(V0, branch, x_end=x_end, keep_sol=True, **kw)


def gauge_kappa(V0):
    """Pure-gauge eigenvalue in the sonic-point gauge: kappa_bar = -dN_bar/dx(0) = -F_N(y0)."""
    y0 = css.sonic_data(V0)
    return -css.coeffs(*y0)["FN"]


def continue_outward(V0, branch, x_max=25.0, K=40, delta=None, rtol=1e-12, atol=1e-15):
    """Integrate the reduced plain system from +delta outward (t -> 0^- is x -> +inf).
    Returns the solve_ivp result with y = (N, W, V); A from the constraint."""
    delta, _, coef = start_data(V0, branch, K, delta)
    y0 = taylor.eval_series(coef, delta)[1:]

    def ev_V(x, y):
        return 0.9999 - abs(y[2])
    ev_V.terminal = True

    def ev_sonic(x, y):
        Am1 = css.Am1_constraint(*y)
        q = css.coeffs(1 + Am1, *y, Am1=Am1)
        D, _, _ = css.det_and_numerators(q)
        return abs(D) / (abs(q["a"] * q["d"]) + abs(q["b"] * q["c"])) - 1e-7
    ev_sonic.terminal = True
    return solve_ivp(css.rhs_plain3, (delta, x_max), y0, method="DOP853", rtol=rtol, atol=atol,
                     events=[ev_V, ev_sonic], dense_output=True)


def count_zeros_V_full(V0, branch, x_end=-15.0, x_max=25.0):
    """Zeros of V on (x_end, x_max): inside (scaled integration) + outside."""
    inner = shoot(V0, branch, x_end=x_end, keep_sol=True)
    outer = continue_outward(V0, branch, x_max=x_max)
    Vo = outer.y[2]
    nz_out = int(np.sum(np.sign(Vo[1:]) != np.sign(Vo[:-1])))
    return inner.n_zeros_V, nz_out, inner, outer
