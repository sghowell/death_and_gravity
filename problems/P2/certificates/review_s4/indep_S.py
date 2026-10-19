"""Independent float check of the S certificates (S4 review).  No p2 code.

For every stored zeta-leaf (both sides) the contour x(t) = X_a + s (Y(t) - t_a)/k of the stored design is rebuilt
from the stored (X_a, k, w, y_a, phi_inf, p_c, p_w); r(t) is integrated in complex double precision (scipy DOP853,
rtol 1e-12) backward from the real exit r(X_a), which is found by an mpmath root solve of r*(r) = X_a.  q^ is then
evaluated on a t-grid of spacing 1/512 on [t_L, t_a] = [-8, t_a] at the 16 vertices of the certificate's
parameter box (Z = 1/zeta^2, eps|nu|, eps^2 nu^2, eps^2; q^ is affine in each), and on the tail t in [-40, -8].
The certified per-eps-box lower bound mu_cells must lie below this float minimum.

Closing boxes (zeta in (0, zeta_1]): no Taylor model.  For rho_T = lam rho_1 (lam on a grid in (0, 1]) the actual
contour x_T(zeta) + s Y(t)/kappa_h is integrated from its real exit, U~ = U(r)/U(r_T) exactly, and W~ is the true
W~(r); lam -> 0 is the exponential throat U~ = e^{2Y}, W~ = 0.
Also: the W1 crossing heights of the c contours, recomputed from the design.
"""
from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from indep_params import PB  # noqa: E402

CERT = Path(__file__).resolve().parents[1] / "part1"
p = PB()
lam3 = float(p.Lam / 3)
roots = [float(r) for r in p.roots]
OmP, km, d = float(p.OmegaP), float(p.km), float(p.d)
kp, kc = float(p.kp), float(p.kc)
X0 = float(p.X0)
WB = 99 / 12500


def Delta(r):
    return -lam3 * (r - roots[0]) * (r - roots[1]) * (r - roots[2]) * (r - roots[3])


def dDelta(r):
    return -4 * lam3 * r ** 3 + 2 * r - 2


def Ffun(r):
    return Delta(r) / r ** 2


def Ufun(r):
    return Delta(r) / r ** 4


def Wtfun(r):
    D, Dp = Delta(r), dDelta(r)
    F = D / r ** 2
    Fr = Dp / r ** 2 - 2 * D / r ** 3
    return F * Fr / r - F / r ** 2 / 4


def hidx(side):
    return 2 if side > 0 else 3


def Delta_s(s, side):
    """Delta at r = r_h + side*s, with the factor (r - r_h) = side*s taken exactly (no cancellation)."""
    i = hidx(side)
    rh = roots[i]
    r = rh + side * s
    prod = -lam3 * side * s
    for j, rj in enumerate(roots):
        if j != i:
            prod = prod * (r - rj)
    return prod, r


def UW_s(s, side):
    D, r = Delta_s(s, side)
    Dp = dDelta(r)
    F = D / r ** 2
    Fr = Dp / r ** 2 - 2 * D / r ** 3
    return D / r ** 4, F * Fr / r - F / r ** 2 / 4, r


# ---- design family -------------------------------------------------------------------------------
def m0(u):
    return 5 / 32 + u / 2 + 15 / 32 * u ** 2 - 5 / 32 * u ** 4 + u ** 6 / 32


def m0p(u):
    return 0.5 + 15 / 16 * u - 10 / 16 * u ** 3 + 3 / 16 * u ** 5


def m0pp(u):
    return 15 / 16 * (1 - u * u) ** 2


def eta_all(t, w):
    """eta, eta', eta'' (numpy arrays)."""
    u = np.clip(-t / w, -1, 1)
    inside = np.abs(t) <= w
    eta = np.where(t <= -w, t, np.where(t >= w, 0.0, -w * m0(u)))
    e1 = np.where(t <= -w, 1.0, np.where(t >= w, 0.0, m0p(u)))
    e2 = np.where(inside, -m0pp(u) / w, 0.0)
    return eta, e1, e2


def phase_all(t, phinf, pc, pw):
    p0 = pc - pw
    v = np.clip((t - p0) / (2 * pw), 0, 1)
    S = 10 * v ** 3 - 15 * v ** 4 + 6 * v ** 5
    dS = 30 * v ** 2 * (1 - v) ** 2
    inside = (t >= p0) & (t <= p0 + 2 * pw)
    phi = np.where(t <= p0, phinf, np.where(t >= p0 + 2 * pw, 0.0, phinf * (1 - S)))
    d1 = np.where(inside, -phinf * dS / (2 * pw), 0.0)
    return phi, d1


def le(cert_printed, float_min):
    """The certificates print bounds with 8 significant digits rounded to nearest (throat.lohi), so a printed
    lower bound can exceed the arb lower endpoint by half a unit in the 8th digit."""
    return cert_printed <= float_min + 1e-7 * abs(cert_printed)


def parse_design(key):
    w, ya = float(Fraction(key[0])), float(Fraction(key[1]))
    phinf = -np.pi / int(key[2].split("/")[1])
    return w, ya, phinf, float(Fraction(key[3])), float(Fraction(key[4]))


def r_real(X):
    return float(p.r_of_x(mp.mpf(X)))


def s_real(X, side):
    """s = side (r - r_h) at the real point r*(r) = X, computed in mpmath (no cancellation)."""
    with mp.workdps(60):
        r = p.r_of_x(mp.mpf(X))
        rh = p.rp if side > 0 else p.rc
        return float(side * (r - rh)), float(r)


def integrate_contour(side, Xa, k, w, ya, t_lo, n_per_unit=512):
    """s(t) = side (r(t) - r_h) along x(t), integrated backward from the real exit as u = log s (complex), so that s
    keeps full relative accuracy as it decays toward the horizon:  du/dt = x'(t) (-Lambda/3) prod_{i != h}(r - r_i)/r^2."""
    ta = w + ya
    sa, ra = s_real(Xa, side)
    i = hidx(side)
    rh = roots[i]
    others = [rj for j, rj in enumerate(roots) if j != i]

    def rhs(t, y):
        u = y[0] + 1j * y[1]
        r = rh + side * np.exp(u)
        _, e1, _ = eta_all(np.array([t]), w)
        xp = side * (1 + 1j * e1[0]) / k
        prod = -lam3
        for rj in others:
            prod = prod * (r - rj)
        v = xp * prod / r ** 2
        return [v.real, v.imag]

    n = int(round((ta - t_lo) * n_per_unit))
    tt = np.linspace(ta, t_lo, n + 1)
    bps = sorted({ta, t_lo} | {x for x in (w, -w) if t_lo < x < ta}, reverse=True)
    us, ts = [], []
    y0 = [float(np.log(sa)), 0.0]
    for a, b in itertools.pairwise(bps):
        sel = tt[(tt <= a) & (tt >= b)]
        sol = solve_ivp(rhs, (a, b), y0, method="DOP853", rtol=1e-12, atol=1e-12, t_eval=sel)
        assert sol.success, sol.message
        ts.append(sol.t)
        us.append(sol.y[0] + 1j * sol.y[1])
        y0 = [sol.y[0][-1], sol.y[1][-1]]
    t = np.concatenate(ts)
    s = np.exp(np.concatenate(us))
    order = np.argsort(t)
    t, s = t[order], s[order]
    keep = np.concatenate([[True], np.diff(t) > 0])
    return t[keep], s[keep], ra


def qhat_vertices(t, Ut_fun_vals, W, w, phinf, pc, pw, k, Zs, ENs, E2N2s, E2s):
    """min over the 16 vertices of q^ at each t.  Ut_fun_vals: U(r(t)) (Z multiplies it) or U~ (Zs = [1])."""
    eta, e1, e2 = eta_all(t, w)
    phi, d1 = phase_all(t, phinf, pc, pw)
    Yp = 1 + 1j * e1
    Ypp = 1j * e2
    G = np.exp(-1j * phi) * Yp ** 2
    K = np.abs(Ypp / Yp - 1j * d1) ** 2 / (4 * np.cos(phi))
    worst = np.full(len(t), np.inf)
    for Z, EN, E2N2, E2 in itertools.product(Zs, ENs, E2N2s, E2s):
        val = (G * (Z * Ut_fun_vals - 1 + 2j * EN + E2N2 + E2 * W)).real - k * k * E2 * K
        worst = np.minimum(worst, val)
    return worst


def check_leaf(side, x):
    w, ya, phinf, pc, pw = parse_design(x["design"])
    ta = w + ya
    Xa, k = float(Fraction(x["X_a"])), float(Fraction(x["k"]))
    zl, zh = (float(Fraction(v)) * OmP for v in x["zeta_over_OmegaP"])
    em = float(Fraction(x["eps_max"]))
    t, sv, ra = integrate_contour(side, Xa, k, w, ya, -40.0)
    U, W, r = UW_s(sv, side)
    res = {"zeta": x["zeta_over_OmegaP"], "design": x["design"], "eps_boxes": []}
    for e in x["eps_boxes"]:
        el, eh = float(Fraction(e["eps"][0])), float(Fraction(e["eps"][1]))
        q = qhat_vertices(t, U, W, w, phinf, pc, pw, k, [1 / zh ** 2, 1 / zl ** 2], [0, eh * d], [0, (eh * d) ** 2],
                          [el * el, eh * eh])
        main = t >= -8
        mu = float(e["mu_cells"].strip("[]").split(",")[0])
        mt = float(e["mu_tail"].strip("[]").split(",")[0])
        res["eps_boxes"].append({"eps": e["eps"], "mu_cells_cert": mu, "float_min_cells": float(q[main].min()),
                                 "t_at_min": float(t[main][np.argmin(q[main])]),
                                 "mu_tail_cert": mt, "float_min_tail": float(q[~main].min()),
                                 "ok": bool(le(mu, q[main].min()) and le(mt, q[~main].min()))})
    # exit and order (float)
    ex = Ufun(ra) / zh ** 2 - 1 - em ** 2 * WB
    res["exit_margin_float"] = float(ex)
    res["exit_order_ok"] = bool(ra < float(p.rP)) if side > 0 else bool(ra > float(p.rP))
    # W1: + side Re x <= X_a < X0; c side crossings
    if side > 0:
        res["w1"] = {"Xa_lt_X0": bool(Xa < X0)}
    else:
        cr = []
        for X, lim in ((X0, 0.1413), (0.0, 2.438)):
            if X <= Xa:
                continue
            ts = ta - k * (X - Xa)
            eta, _, _ = eta_all(np.array([ts]), w)
            cr.append({"line": X, "t_star": ts, "abs_Im_x": float(abs(eta[0] / k)), "ok": bool(abs(eta[0] / k) < lim)})
        res["w1"] = {"crossings": cr}
    # tail region: X_L = X_a + s (t_L - t_a)/k
    XL = Xa + side * (-8 - ta) / k
    res["X_L"] = XL
    res["tail_region_in_strip"] = bool(XL < X0) if side > 0 else bool(XL > 0)
    res["ok"] = all(e["ok"] for e in res["eps_boxes"]) and ex > 0 and res["exit_order_ok"] and res["tail_region_in_strip"]
    return res


def closing_leaf(side, x, lams):
    """Exact (Taylor-model-free) q^ for rho_T = lam rho_1."""
    w, ya, phinf, pc, pw = parse_design(x["design"])
    ta = w + ya
    kap = kp if side > 0 else kc
    rh = float(p.rp) if side > 0 else float(p.rc)
    ch = float(p.cp) if side > 0 else float(p.cc)
    rho1 = float(mp.mpf(x["near_horizon"]["rho_1"].split("+/-")[0].strip("[ ")))
    em = float(Fraction(x["eps_max"]))
    out = []
    worst_all = None
    for lam in lams:
        rhoT = lam * rho1
        # rho = e^{+2 kap (x - c)} (side +), e^{-2 kap (x - c)} (side c)  =>  x_T = c + side log(rho_T)/(2 kap)
        xT = ch + side * float(mp.log(rhoT)) / (2 * kap)
        sT, rT = s_real(xT, side)
        zeta2 = UW_s(np.array([sT]), side)[0][0]
        Xa = xT + side * ta / kap
        t, sv, ra = integrate_contour(side, Xa, kap, w, ya, -40.0)
        Uv, W, r = UW_s(sv, side)
        Ut = Uv / zeta2
        q = qhat_vertices(t, Ut, W, w, phinf, pc, pw, kap, [1.0], [0, em * d], [0, (em * d) ** 2], [0, em * em])
        main = t >= -8
        ex = Ut[-1].real - 1 - em ** 2 * WB
        rec = {"lam": lam, "zeta_over_OmegaP": float(np.sqrt(zeta2) / OmP), "float_min_cells": float(q[main].min()),
               "t_at_min": float(t[main][np.argmin(q[main])]), "float_min_tail": float(q[~main].min()), "exit_margin": float(ex)}
        out.append(rec)
        worst_all = rec["float_min_cells"] if worst_all is None else min(worst_all, rec["float_min_cells"])
    # lam -> 0: the exponential throat U~ = e^{2Y}, W~ = 0
    t = np.linspace(-40, ta, int((ta + 40) * 512) + 1)
    eta, _, _ = eta_all(t, w)
    Ut = np.exp(2 * (t + 1j * eta))
    q = qhat_vertices(t, Ut, np.zeros_like(t), w, phinf, pc, pw, kap, [1.0], [0, em * d], [0, (em * d) ** 2], [0, em * em])
    main = t >= -8
    lim = {"float_min_cells": float(q[main].min()), "t_at_min": float(t[main][np.argmin(q[main])]),
           "float_min_tail": float(q[~main].min()), "exit_margin": float(np.exp(2 * ta) - 1 - em ** 2 * WB)}
    mu = float(x["mu_cells"].strip("[]").split(",")[0])
    mt = float(x["mu_tail"].strip("[]").split(",")[0])
    ok = all(le(mu, o["float_min_cells"]) and le(mt, o["float_min_tail"]) and o["exit_margin"] > 0 for o in out) and \
        le(mu, lim["float_min_cells"]) and le(mt, lim["float_min_tail"])
    return {"mu_cells_cert": mu, "mu_tail_cert": mt, "rho_T_samples": out, "limit_lam_0": lim,
            "float_min_over_samples": min(worst_all, lim["float_min_cells"]), "ok": bool(ok)}


def main():
    out = {}
    for side, tag in ((1, "plus"), (-1, "c")):
        cert = json.loads((CERT / f"S_PB_{tag}.json").read_text())
        rows = []
        for x in cert["leaves"]:
            if x.get("closing"):
                lams = [1.0, 0.7, 0.5, 0.3, 0.1, 0.03, 1e-2, 1e-3, 1e-5]
                c = closing_leaf(side, x, lams)
                rows.append({"zeta": x["zeta_over_OmegaP"], "closing": True, **c})
                print(tag, "closing", c["mu_cells_cert"], c["float_min_over_samples"], c["ok"], flush=True)
            else:
                rr = check_leaf(side, x)
                rows.append(rr)
                e0 = rr["eps_boxes"][0]
                print(tag, rr["zeta"], e0["mu_cells_cert"], round(e0["float_min_cells"], 6), e0["t_at_min"], rr["ok"],
                      rr.get("w1"), flush=True)
        out[tag] = {"leaves": rows, "all_ok": all(r["ok"] for r in rows)}
    Path(sys.argv[1] if len(sys.argv) > 1 else "indep_S.json").write_text(json.dumps(out, indent=1) + "\n")
    print("plus all ok", out["plus"]["all_ok"], "c all ok", out["c"]["all_ok"])


if __name__ == "__main__":
    main()
