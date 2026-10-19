"""Independent float check of the T and A certificates (S4 review).

No p2 code.  r(s) on S2's contour x = x_P + e^{i pi/4} s is integrated in complex double precision (scipy DOP853,
rtol 1e-13) from r(0) = r_P.  On a fine s-grid (spacing 1/256 on [-64, 32], plus coarser tails to -160/+80):
  (a) the splitting bound F = E0(zeta_lo)/h0 + E1(zeta_hi, nu_lo) + h0 min(0, E2(nu_lo)) per stored zeta-box,
      against the certified mu (certified must be <= the float minimum), and E0 >= 0;
  (b) S2's E itself, computed from V_l = F(L/r^2 + F'/r) directly (no splitting), at actual l and sampled
      (zeta, nu): h E >= mu_cert(box) must hold pointwise;
  (c) A: the same at zeta_A, and the monotonicity dq/dRe(w) >= 0 on a grid of Re(w) >= zeta_A/h.
"""
from __future__ import annotations

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
Q2 = float(p.Q2)
roots = [float(r) for r in p.roots]
rP, xP = float(p.rP), float(p.xP)
OmP, g0, km, d = float(p.OmegaP), float(p.gamma0), float(p.km), float(p.d)
astar = float(p.astar)
h0 = 2 / 41
E4 = np.exp(1j * np.pi / 4)


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


def Vl(r, l):
    D, Dp = Delta(r), dDelta(r)
    F = D / r ** 2
    Fr = Dp / r ** 2 - 2 * D / r ** 3
    return F * (l * (l + 1) / r ** 2 + Fr / r)


def contour(s_end, n):
    """r(s) for s from 0 to s_end on n+1 equispaced points."""
    def rhs(s, y):
        r = y[0] + 1j * y[1]
        v = E4 * Ffun(r)
        return [v.real, v.imag]
    ss = np.linspace(0, s_end, n + 1)
    sol = solve_ivp(rhs, (0, s_end), [rP, 0.0], method="DOP853", rtol=1e-13, atol=1e-15, t_eval=ss)
    assert sol.success
    return ss, sol.y[0] + 1j * sol.y[1]


def grid():
    sL, rL = contour(-64.0, 64 * 256)
    sR, rR = contour(32.0, 32 * 256)
    s = np.concatenate([sL[::-1], sR[1:]])
    r = np.concatenate([rL[::-1], rR[1:]])
    # tails (coarse)
    sLt, rLt = contour(-160.0, 160 * 8)
    sRt, rRt = contour(80.0, 80 * 8)
    mL, mR = sLt < -64, sRt > 32
    st = np.concatenate([sLt[mL][::-1], sRt[mR]])
    rt = np.concatenate([rLt[mL][::-1], rRt[mR]])
    return s, r, st, rt


class Design:
    def __init__(self, phi0, phinf=-23 / 16, s1=7.0, eta=0.5, k=None):
        self.phi0, self.phinf, self.s1, self.eta, self.k = phi0, phinf, s1, eta, k

    def phase(self, s):
        s14 = self.s1 ** 4
        D = s14 + s ** 4
        a = self.phinf - self.phi0
        phi = self.phi0 + a * s ** 4 / D
        d1 = a * 4 * s ** 3 * s14 / D ** 2
        d2 = a * 4 * s14 * s ** 2 * (3 * s14 - 5 * s ** 4) / D ** 3
        return phi, d1, d2

    def gamma(self, s):
        if self.k is None:
            z = np.zeros_like(s)
            return z, z, z
        a = np.abs(s)
        sg = np.sign(s)
        th = np.where(a <= 2, 1.0, np.where(a >= 4, 0.0, (1 + np.cos(np.pi * (a - 2) / 2)) / 2))
        thp = np.where((a > 2) & (a < 4), -(np.pi / 4) * np.sin(np.pi * (a - 2) / 2) * sg, 0.0)
        thpp = np.where((a > 2) & (a < 4), -(np.pi ** 2 / 8) * np.cos(np.pi * (a - 2) / 2), 0.0)
        k = self.k
        return k * s * th, k * (th + s * thp), k * (2 * thp + s * thpp)

    def c_and_cp(self, s):
        phi, d1, d2 = self.phase(s)
        sgn = np.sign(d1)
        c = np.cos(phi) - self.eta * np.abs(d1)
        cp = -np.sin(phi) * d1 - self.eta * sgn * d2
        return c, cp


def E_parts(D, s, r, zeta_lo, zeta_hi, nu_lo):
    phi, d1, _ = D.phase(s)
    c, cp = D.c_and_cp(s)
    G, Gp, _ = D.gamma(s)
    rot = 1j * np.exp(-1j * phi)
    U, W = Ufun(r), Wtfun(r)
    E0 = (rot * (U - zeta_lo ** 2)).real - c * G * G
    E1 = 2 * zeta_hi * nu_lo * np.cos(phi) + cp * G + c * Gp
    E2 = (rot * (W + nu_lo ** 2)).real - np.abs(d1) / (4 * D.eta)
    return E0, E1, E2, c


def E_direct(D, s, r, l, zeta, nu):
    """S2's E = q + c'G' - c G'^2 + c G'' with G' = Gamma/h, from V_l directly."""
    h = 1 / (l + 0.5)
    w = zeta / h + 1j * nu
    phi, d1, _ = D.phase(s)
    c, cp = D.c_and_cp(s)
    G, Gp, _ = D.gamma(s)
    q = (np.exp(-1j * phi) * 1j * (Vl(r, l) - w * w)).real - np.abs(d1) / (4 * D.eta)
    Gx, Gxx = G / h, Gp / h
    return q + cp * Gx - c * Gx * Gx + c * Gxx


def main():
    s, r, st, rt = grid()
    out = {"grid": {"spacing": "1/256 on [-64, 32]; 1/8 on the tails to -160 and 80"}}
    TD = Design(0.0, k=0.95 * 32178552 / 2 ** 30)
    AD = Design(-0.2)
    cert = json.loads((CERT / "T_PB.json").read_text())
    rows = []
    worst_ratio = None
    for b in cert["boxes"]:
        zl = float(Fraction(b["zeta_over_OmegaP"][0])) * OmP
        zh = float(Fraction(b["zeta_over_OmegaP"][1])) * OmP
        E0, E1, E2, c = E_parts(TD, s, r, zl, zh, -d)
        F = E0 / h0 + E1 + h0 * np.minimum(0, E2)
        E0t, E1t, E2t, ct = E_parts(TD, st, rt, zl, zh, -d)
        Ft = E0t / h0 + E1t + h0 * np.minimum(0, E2t)
        mu_c = float(mp.mpf(b["mu_cert"].split("+/-")[0].strip("[ ")))
        i = int(np.argmin(F))
        row = {"zeta": b["zeta_over_OmegaP"], "mu_cert_over_astar": mu_c / astar,
               "float_min_F_over_astar": float(F.min() / astar), "at_s": float(s[i]),
               "float_min_F_tails_over_astar": float(Ft.min() / astar),
               "min_E0": float(E0.min()), "min_E0_tails": float(E0t.min()), "min_c": float(min(c.min(), ct.min())),
               "cert_le_float": bool(mu_c <= F.min() and mu_c <= Ft.min())}
        rows.append(row)
        ratio = mu_c / F.min()
        worst_ratio = ratio if worst_ratio is None else max(worst_ratio, ratio)
    out["T_boxes"] = rows
    out["T_all_cert_le_float"] = all(x["cert_le_float"] for x in rows)
    out["T_min_E0_over_grid"] = min(x["min_E0"] for x in rows)
    # (b) direct E at actual l: h E >= mu_cert of the box containing zeta
    rng = np.random.default_rng(20261019)
    direct = []
    for l in (20, 21, 25, 37, 50, 100, 400, 2000, 100000):
        for _ in range(6):
            zf = rng.uniform(0.8, 2.0)
            nu = -d * rng.choice([0.0, 0.5, 1.0, rng.uniform()])
            box = [bb for bb in cert["boxes"] if float(Fraction(bb["zeta_over_OmegaP"][0])) <= zf <= float(Fraction(bb["zeta_over_OmegaP"][1]))][0]
            mu_c = float(mp.mpf(box["mu_cert"].split("+/-")[0].strip("[ ")))
            h = 1 / (l + 0.5)
            hE = h * E_direct(TD, s, r, l, zf * OmP, nu)
            hEt = h * E_direct(TD, st, rt, l, zf * OmP, nu)
            m = float(min(hE.min(), hEt.min()))
            direct.append({"l": l, "zeta_over_OmegaP": zf, "nu_over_km": nu / km, "min_hE_over_astar": m / astar,
                           "mu_cert_over_astar": mu_c / astar, "ok": bool(m >= mu_c)})
    out["T_direct_E"] = direct
    out["T_direct_all_ok"] = all(x["ok"] for x in direct)
    # splitting identity check at one point: h^2 E == E0 + h E1 + h^2 E2 (zeta_lo = zeta_hi, nu_lo = nu)
    l, zf, nu = 37, 0.93, -0.017
    h = 1 / (l + 0.5)
    E0, E1, E2, _ = E_parts(TD, s, r, zf * OmP, zf * OmP, nu)
    Ed = E_direct(TD, s, r, l, zf * OmP, nu)
    out["T_splitting_max_rel_diff"] = float(np.max(np.abs(h * h * Ed - (E0 + h * E1 + h * h * E2)) / (np.abs(h * h * Ed) + 1e-12)))
    # A
    certA = json.loads((CERT / "A_PB.json").read_text())
    bA = certA["boxes"][0]
    zA = 1.5 * OmP
    E0, E1, E2, c = E_parts(AD, s, r, zA, zA, -d)
    F = E0 / h0 + E1 + h0 * np.minimum(0, E2)
    E0t, E1t, E2t, ct = E_parts(AD, st, rt, zA, zA, -d)
    Ft = E0t / h0 + E1t + h0 * np.minimum(0, E2t)
    muA = float(mp.mpf(bA["mu_cert"].split("+/-")[0].strip("[ ")))
    out["A"] = {"mu_cert_over_astar": muA / astar, "float_min_F_over_astar": float(min(F.min(), Ft.min()) / astar),
                "float_min_E0": float(min(E0.min(), E0t.min())), "min_c": float(min(c.min(), ct.min())),
                "cert_le_float": bool(muA <= min(F.min(), Ft.min()))}
    # A monotonicity: q(s; Re w) nondecreasing along Re w >= zeta_A/h (Gamma = 0, so E = q)
    mono = []
    for l in (20, 50, 1000):
        h = 1 / (l + 0.5)
        for nu in (0.0, -d / 2, -d):
            prev = None
            okm = True
            minE = None
            for fac in np.linspace(1.0, 6.0, 41):
                zeta = zA * fac
                Ed = E_direct(AD, np.concatenate([s, st]), np.concatenate([r, rt]), l, zeta, nu)
                if prev is not None and np.any(Ed < prev - 1e-12 * np.abs(prev)):
                    okm = False
                prev = Ed
                minE = Ed.min() if minE is None else min(minE, Ed.min())
            mono.append({"l": l, "nu_over_km": nu / km, "nondecreasing_on_grid": okm, "min_hE_over_astar": float(h * minE / astar)})
    out["A_monotone"] = mono
    out["A_monotone_all"] = all(x["nondecreasing_on_grid"] for x in mono)
    print(json.dumps({k: v for k, v in out.items() if k not in ("T_direct_E",)}, indent=1))
    print("direct E rows:", len(direct), "all ok:", out["T_direct_all_ok"],
          "min margin ratio (min hE / mu_cert):", min(x["min_hE_over_astar"] / x["mu_cert_over_astar"] for x in direct))
    Path(sys.argv[1] if len(sys.argv) > 1 else "indep_TA.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
