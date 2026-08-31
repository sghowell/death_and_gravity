"""Theorem C, S4-1: FLOAT validation of the exact (chi, eta) reduction (tc_reduce) and
the route-(b) sign survey.  Nothing here is validated arithmetic.

Modes (run from the repo root, PYTHONPATH=problems/P4/src):
    uv run python -m p4.validated.tc_reduced_eigs check    # reduced vs full system
    uv run python -m p4.validated.tc_reduced_eigs eigs     # kappa1, kbar reproduction
    uv run python -m p4.validated.tc_reduced_eigs survey   # sign survey on the tube grid
The reduced 2x2 system is integrated with the coefficients evaluated from the exact
fmpq_mpoly data of tc_reduce; for real kappa in (0.35, 1.05) the N-clock chart has an
apparent singularity at F_N(x) = -kappa, crossed by a complex-x contour detour (the
background and the true solution are analytic in x; only the chart coefficients have
the pole).  E_red(kappa) := chi(x_end) e^{3 x_end}; its zeros near kappa1 and kbar are
the reduced eigenvalues (at kbar the admissible solution is the gauge mode, chi == 0).
Survey grid: the certified tube's step abscissae (P4_TUBE_CACHE, decoded as floats);
fallback: a dense grid on [-8, -0.05].
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

from p4 import css, perturb, taylor
from p4.validated.tc_reduce import Reduction

V0 = 0.112439401388092
DELTA = 0.05
BUMP = (-0.25, -0.05, 0.04)           # complex-x detour (a, b, height); ends exactly at x0
_red = _F = None


def compile_poly(p):
    ex, co = [], []
    for exps, coef in p.terms():
        ex.append(list(map(int, exps)))
        co.append(float(coef.p) / float(coef.q))
    ex, co = np.array(ex, dtype=np.int64), np.array(co, dtype=complex)
    used = [i for i in range(10) if ex[:, i].any()]

    def f(vals):
        acc = co.copy()
        for i in used:
            pw = np.empty(ex[:, i].max() + 1, dtype=complex)
            pw[0] = 1.0
            for j in range(1, len(pw)):
                pw[j] = pw[j - 1] * vals[i]
            acc = acc * pw[ex[:, i]]
        return acc.sum()
    return f


def setup():
    global _red, _F
    if _F is None:
        _red = Reduction()
        rp = _red.reduce_pair()
        _F = {n: compile_poly(p) for n, p in
              (("Achi", rp["chi"][0]), ("Bchi", rp["chi"][1]), ("Aeta", rp["eta"][0]),
               ("Beta", rp["eta"][1]), ("D2", rp["D2"]), ("PV", _red.PV), ("PW", _red.PW),
               ("e", _red.e))}
    return _F


def vals_of(u, kap):
    return [0.0, u[0], u[1], u[2], u[3], 0, 0, 0, 0, kap]


def xpath(s, h):
    a, b, _ = BUMP
    if h == 0.0 or not (a < s < b):
        return s, 1.0
    u = (s - a) / (b - a)
    return s + 1j * h * np.sin(np.pi * u), 1.0 + 1j * h * np.pi / (b - a) * np.cos(np.pi * u)


def rhs(s_, y, kap, h=0.0):
    F = setup()
    x, dxds = xpath(s_, h)
    Nh, Wh, Vh = y[0], y[1], y[2]
    e_x = np.exp(x)
    N, W, V = Nh / e_x, Wh * e_x**2, Vh * e_x
    Am1 = css.Am1_constraint(N, W, V)
    q = css.coeffs(1 + Am1, N, W, V, Am1=Am1)
    Delta, P, Q = css.det_and_numerators(q)
    vals = vals_of([1 + Am1, N, W, V], kap)
    D2 = F["D2"](vals)
    dchi = (F["Achi"](vals) * y[3] + F["Bchi"](vals) * y[4]) / D2
    deta = (F["Aeta"](vals) * y[3] + F["Beta"](vals) * y[4]) / D2
    return np.array([Nh * (Am1 - 2 * W / 3), (P / Delta) / e_x**2 - 2 * Wh,
                     (Q / Delta) / e_x - Vh, dchi, deta], dtype=complex) * dxds


def chi_eta_of_p(u, kap, p):
    F = setup()
    vals = vals_of(u, kap)
    e_, PV, PW = F["e"](vals), F["PV"](vals), F["PW"](vals)
    return e_ * p[3] - PV * p[1], e_ * p[2] - PW * p[1]


def reduced_final(kap, x_end=-8.0, rtol=1e-12, detour=None, bg=None):
    kap = complex(kap)
    bg = taylor.background_series(V0, 1, K=36) if bg is None else bg
    if detour is None:
        detour = BUMP[2] if (abs(kap.imag) < 1e-3 and 0.3 < kap.real < 1.05) else 0.0
    pc = taylor.perturbation_series(bg, kap, K=36)
    u0, p0 = taylor.eval_series(bg, -DELTA), taylor.eval_series(pc, -DELTA)
    chi0, eta0 = chi_eta_of_p(u0, kap, p0)
    yh = css.to_scaled(-DELTA, u0)[1:]
    y0 = np.array([yh[0], yh[1], yh[2], chi0, eta0], dtype=complex)
    sol = solve_ivp(rhs, (-DELTA, x_end), y0, args=(kap, detour), method="DOP853",
                    rtol=rtol, atol=1e-300)
    if sol.status != 0:
        raise RuntimeError(sol.message)
    return sol.y[3, -1], sol.y[4, -1]


def direct_final(kap, x_end=-8.0):
    prob = perturb.Problem(V0, x_end=x_end, delta=DELTA, K=36)
    x, ybg, Ap, Np, Wp, Vp = prob.perturbation_profile(kap)
    e_x = np.exp(x[-1])
    N, W, V = ybg[0][-1] / e_x, ybg[1][-1] * e_x**2, ybg[2][-1] * e_x
    u = [1 + css.Am1_constraint(N, W, V), N, W, V]
    return chi_eta_of_p(u, kap, [Ap[-1], Np[-1], Wp[-1], Vp[-1]])


def find_zero(k0, k1, x_end=-11.0, bg=None):
    E = lambda k: reduced_final(k, x_end=x_end, bg=bg)[0] * np.exp(3 * x_end)
    a, b = complex(k0), complex(k1)
    fa, fb = E(a), E(b)
    for _ in range(60):
        c = b - fb * (b - a) / (fb - fa)
        a, fa, b, fb = b, fb, c, E(c)
        if abs(b - a) < 1e-13 * abs(b):
            break
    return b


def tube_grid():
    path = os.environ.get("P4_TUBE_CACHE", "")
    if path and os.path.exists(path):
        with open(path) as fh:
            data = json.load(fh)
        xs = sorted(float(int(s["x"][0])) * 2.0 ** int(s["x"][1]) for s in data["steps"])
        return np.array(xs), f"tube({len(xs)} steps)"
    return np.linspace(-8.0, -0.05, 345), "dense fallback"


def survey():
    F = setup()
    xs, tag = tube_grid()
    bg = taylor.background_series(V0, 1, K=36)
    from p4 import shoot
    sh = shoot.shoot(V0, 1, x_end=-9.0, keep_sol=True)

    def u_of(x):
        if x > -0.06:
            return taylor.eval_series(bg, x)
        y = sh.sol.sol(x)
        N, W, V = y[0] * np.exp(-x), y[1] * np.exp(2 * x), y[2] * np.exp(x)
        return np.array([1 + css.Am1_constraint(N, W, V), N, W, V])

    def TU(x, kap, h=1e-6):
        def abcd(xx):
            u = u_of(xx)
            vals = vals_of(u, kap)
            D2 = F["D2"](vals)
            return np.array([F[n](vals) for n in ("Achi", "Bchi", "Aeta", "Beta")]) / D2
        a, b, c2, d2 = abcd(x)
        ap, bp = (abcd(x + h)[:2] - abcd(x - h)[:2]) / (2 * h)
        T = a + d2 + bp / b
        return T, ap + b * c2 - a * (T - a), b

    rows = []
    for x in xs:
        A, N, W, V = u_of(x)
        du = css.rhs_plain(x, [A, N, W, V])
        D = 3 * N**2 * V**2 - N**2 + 4 * N * V - V**2 + 3
        Dp = (6 * N * V**2 - 2 * N + 4 * V) * du[1] + (6 * N**2 * V + 4 * N - 2 * V) * du[3]
        q = css.coeffs(A, N, W, V)
        C = np.linalg.solve(np.array([[q["a"], q["b"]], [q["c"], q["d"]]]),
                            np.array([[q["sa"], q["sb"]], [q["sc"], q["sd"]]]))
        trC, detC = float(np.trace(C)), float(np.linalg.det(C))
        vrel = (1 + N * V) / (N + V)
        rows.append({"x": float(x), "w": -D, "dw_over_w": Dp / D, "trC": trC,
                     "detC": detC, "threemV2": 3 - V * V,
                     "Hq": (3 - V * V) + D * trC**2 / 4,
                     "wGHJS": 1 / 3 - vrel**2, "factor": 3 * (N + V) ** 2})
    r = {k: np.array([row[k] for row in rows]) for k in rows[0]}
    ok_w = bool((r["w"] > 0).all())
    ok_dw = bool((r["dw_over_w"] < 0).all())
    ok_detC = bool((r["detC"] < 0).all())
    idw = np.abs(r["w"] - r["factor"] * r["wGHJS"]).max()
    print(f"survey grid: {tag};  w := -D  (weight; = 3(N+V)^2 (1/3 - v_rel^2), id.err {idw:.1e})")
    print(f"  w > 0 on grid: {ok_w};  range {r['w'].min():.4g} .. {r['w'].max():.4g}")
    print(f"  w'/w < 0 on grid: {ok_dw};  range {r['dw_over_w'].min():.4g} .. {r['dw_over_w'].max():.4g}")
    print(f"  detC < 0, so a0's kappa^2-density w detC = -(3-V^2) is sign-definite: {ok_detC}; "
          f"3-V^2 in {r['threemV2'].min():.4f}..{r['threemV2'].max():.4f}")
    neg = r["Hq"] < 0
    if neg.any():
        xH = r["x"][neg].min()
        single = bool((np.diff(neg.astype(int)) != 0).sum() <= 1 + (neg[0] != neg[-1]))
        inn = r["Hq"][r["x"] < xH]
        print(f"  Hq := (3-V^2) + D trC^2/4 (square-completion margin): Hq > 0 for x <= x_H = {xH:.3f},"
              f" Hq < 0 on ({xH:.3f}, 0) (single sign change: {single});")
        print(f"    the sonic shoulder, where trC ~ mu- ~ (1-sigma0)/x; interior min Hq = {inn.min():.4f}"
              f" (at x = {r[chr(34)+chr(34)] if False else r[chr(120)][r[chr(120)] < xH][np.argmin(inn)]:.3f})"
              f" -- the region GHJS treat with the w'/w commutator damping")
    else:
        print("  Hq > 0 on the whole grid")
    sub = xs[::len(xs) // 40]
    for tau in (15.0, 30.0, 60.0):
        dev_in, dev_all, sc_in, bmin = 0.0, 0.0, 0.0, np.inf
        for x in sub:
            _T, U, b = TU(float(x), 1j * tau)
            u2lim = U / (1j * tau) ** 2
            A, N, W, V = u_of(float(x))
            q = css.coeffs(A, N, W, V)
            C = np.linalg.solve(np.array([[q["a"], q["b"]], [q["c"], q["d"]]]),
                                np.array([[q["sa"], q["sb"]], [q["sc"], q["sd"]]]))
            d = abs(u2lim / (-np.linalg.det(C)) - 1)
            dev_all = max(dev_all, d)
            if x <= -0.5:
                dev_in = max(dev_in, d)
                sc_in = max(sc_in, d * tau * abs(x))
            bmin = min(bmin, abs(b) / tau)
        print(f"  kappa = {tau:.0f}i: |U/kappa^2/(-detC) - 1|: max {dev_all:.3f} (all x), "
              f"{dev_in:.4f} (x <= -0.5, scaled by tau|x|: {sc_in:.2f}); min |b|/tau = {bmin:.3g}")
    return r


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("check", "all"):
        for kap in (2.5, 2.8105525488, 0.34, 5.0 + 7.0j):
            cr, er = reduced_final(kap)
            cd, ed = direct_final(kap)
            print(f"kappa={kap}: |chi|={abs(cr):.3e}  rel(chi)={abs(cr - cd) / abs(cd):.2e}  "
                  f"rel(eta)={abs(er - ed) / abs(ed):.2e}")
    if which in ("eigs", "all"):
        bg = taylor.background_series(V0, 1, K=36)
        z1 = find_zero(2.81055, 2.81056, bg=bg)
        zg = find_zero(0.3556, 0.3558, bg=bg)
        print(f"kappa1_red = {z1.real:.13f}  (certified 2.8105525488271; diff {abs(z1 - 2.8105525488271472):.2e})")
        print(f"kbar_red   = {zg.real:.13f}  (certified 0.3556992037110; diff {abs(zg - 0.355699203710964):.2e})")
    if which in ("survey", "all"):
        survey()
