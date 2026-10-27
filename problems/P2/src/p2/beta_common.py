"""S6 (P2-6, the two-sided beta at P_B): constants and shared helpers.

Pre-registration: notes/s6-preregistration.md (c98ca16d, Amendments A1-A6 at f82114bf).
New modules only: the S2-S5 modules and certificates are imported, never edited.  Where an S4
routine reads a module-level design constant (the S regime's design grid), S6 overrides it at run
time in its own worker processes and records the values in its certificates.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx, fmpq

from .arbutil import A, git_commit
from .params import Point

ROOT = Path(__file__).resolve().parents[2]          # problems/P2
CERT = ROOT / "certificates" / "beta"
RES = ROOT / "results" / "beta"
PREC = 128
EDGE = Fraction(23, 25)                              # e = 0.92 (FORMULATION P2-6)
H0 = Fraction(2, 41)                                 # l0 = 20
DEN = 2 ** 30
MAIN_PIN = "300c22a7"                                # main with S5 (A2)

# part (i), T redesigned (pre-registration O4): S4's DESIGN_T with these three changes
T_DESIGN = {
    "name": "T6", "s1": "12", "p": 4, "phi_inf": "-1", "phi0": "0", "eta": "1/2",
    "eps_g": "1/20", "a_d": "32178552/1073741824", "taper": ["3", "6"],
    "zeta_over_OmegaP": ["4/5", "13/10"], "zeta_box_over_OmegaP": "1/20",
    "max_depth": 12, "target_over_a_star": "1/10",
}
T_FALLBACK = [("10", "-1"), ("14", "-1"), ("12", "-5/4"), ("12", "-3/4")]   # (s1, phi_inf), float maximin
A_DESIGN = {
    "name": "A6", "s1": "7", "p": 4, "phi_inf": "-23/16", "phi0": "-1/5", "eta": "1/2",
    "eps_g": None, "a_d": "0", "taper": ["2", "4"],
    "zeta_over_OmegaP": ["6/5", "6/5"], "zeta_box_over_OmegaP": "0",
    "max_depth": 12, "target_over_a_star": "1/10",
}
# part (i), S: S4's grid plus these values (pre-registration O4)
S_GRID_EXT = {"w": ["1/2", "3/2", "5/2"], "y_a": ["1/40", "3/10", "3/5"], "phi_n": [5, 6],
              "phase_centre": ["-3/2", "1/2"], "phase_width": ["3/4", "2"]}
# part (i), L: the band rows (kappa_- units)
L_RE = (Fraction(0), Fraction(8))
L_IM = (Fraction(-19, 20), Fraction(-11, 20))
L_BOX = Fraction(1, 20)
L_MIN_BOX = Fraction(1, 80)
# part (ii'), A at fixed h: the zeta_A rule
ZA_GRID = [Fraction(110 + 5 * k, 100) for k in range(219)]          # 1.10, 1.15, ..., 12.00
ZA_FLOAT_MARGIN = 0.02
ZA_RETRIES = 4
L_SMALL = list(range(20))


def deep_point(name: str = "PB") -> Point:
    """The point with the strip depth d = e kappa_- (every S2/S4 regime reads pt.beta)."""
    pt = Point(name, 256)
    pt.beta = A(EDGE) * pt.km
    return pt


def band_bounds(pt: Point):
    """(b', b0) exact dyadics: b' >= e kappa_- (rounded up), 0 < b0 <= kappa_-/2 (rounded down)."""
    e = A(EDGE) * pt.km
    bp = fmpq(math.ceil(float(e.upper()) * DEN), DEN)
    while not bool(arb(bp) >= e):
        bp += fmpq(1, DEN)
    h = pt.km / 2
    b0 = fmpq(math.floor(float(h.lower()) * DEN), DEN)
    while not bool(arb(b0) <= h):
        b0 -= fmpq(1, DEN)
    return bp, b0


def dy_up(x: arb) -> fmpq:
    q = fmpq(math.ceil(float(x.upper()) * DEN), DEN)
    while not bool(arb(q) >= x):
        q += fmpq(1, DEN)
    return q


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1) + "\n")


def commit_tag() -> str:
    return git_commit()


def use_extended_s_grid():
    """S4's S grid plus S_GRID_EXT (run-time override in this process; throat.py is not edited)."""
    from . import throat as T
    base = {"w": [str(x) for x in T.GRID_W], "y_a": [str(x) for x in T.GRID_YA], "phi_n": list(T.GRID_PHI),
            "phase_centre": [str(x) for x in T.GRID_PC], "phase_width": [str(x) for x in T.GRID_PW]}
    if not getattr(T, "_S6_EXTENDED", False):
        T.GRID_W = T.GRID_W + [Fraction(x) for x in S_GRID_EXT["w"]]
        T.GRID_YA = T.GRID_YA + [Fraction(x) for x in S_GRID_EXT["y_a"]]
        T.GRID_PHI = T.GRID_PHI + list(S_GRID_EXT["phi_n"])
        T.GRID_PC = T.GRID_PC + [Fraction(x) for x in S_GRID_EXT["phase_centre"]]
        T.GRID_PW = T.GRID_PW + [Fraction(x) for x in S_GRID_EXT["phase_width"]]
        T._S6_EXTENDED = True
    return {"base": base, "added": S_GRID_EXT,
            "used": {"w": [str(x) for x in T.GRID_W], "y_a": [str(x) for x in T.GRID_YA], "phi_n": list(T.GRID_PHI),
                     "phase_centre": [str(x) for x in T.GRID_PC], "phase_width": [str(x) for x in T.GRID_PW]}}


# ---------------------------------------------------------------------------------------------
# L: the pole-free series R~ = (1 + 2 rho) R (W6'), and a checker that uses it in every box
# ---------------------------------------------------------------------------------------------
def frob_tau_pf(A_, B_, C_, rho, tau, min_terms=30):
    """lowfreq.frob_tau for R~ = (2 rho + 1) R: a_0 = 2 rho + 1 and the n = 1 step with the factor
    (2 rho + 1) cancelled exactly, a_1 = -(A_3 rho (rho - 1) + B_2 rho + C_1)/A_2; n >= 2 unchanged.
    R~ is the solution regular at r_h for every omega (a nonzero multiple of R off omega0 = -i kappa_+,
    the exponent-(rho + 1) solution at omega0)."""
    from .affine import AF, up
    from .arbw import _geom_tails, _solve_q
    from .lowfreq import N_MAX
    A2 = A_[2]
    ta = arb(tau)
    a0 = rho * 2 + 1
    term0 = rho.const(0)
    if 3 < len(A_):
        term0 = term0 + A_[3] * (rho * (rho + (-1)))
    if 2 < len(B_):
        term0 = term0 + B_[2] * rho
    if 1 < len(C_):
        term0 = term0 + C_[1]
    a1 = -(term0 / A2)
    a = [a0, a1]
    S = a0 + a1 * ta
    Sd = a0 * rho * (1 / acb(ta)) + a1 * (rho + 1)
    tk = arb(ta)
    small = 0
    eps = arb(2) ** (-(ctx.prec - 12))
    for n in range(2, N_MAX):
        acc = rho.const(0)
        for k in range(max(0, n - 6), n):
            j2, j1, j0 = n - k + 2, n - k + 1, n - k
            term = rho.const(0)
            if j2 < len(A_):
                term = term + A_[j2] * ((rho + k) * (rho + (k - 1)))
            if j1 < len(B_):
                term = term + B_[j1] * (rho + k)
            if j0 < len(C_):
                term = term + C_[j0]
            acc = acc + a[k] * term
        an = -(acc / (A2 * ((rho * 2 + n) * n)))
        a.append(an)
        tk = tk * ta
        S = S + an * tk
        Sd = Sd + an * ((rho + n) * (tk / ta))
        if n >= min_terms and an.absu() * tk < eps * abs(S.c0).lower():
            small += 1
            if small >= 8:
                break
        else:
            small = 0
    M = len(a)
    Rr = rho.absu()
    Nn = arb(M)
    if not (Nn > 2 * Rr + 1):
        raise ArithmeticError("frob_tau_pf: N too small relative to |rho|")
    A2l = abs(A2.c0).lower() - A2.dev()
    if not (A2l > 0):
        raise ArithmeticError("frob_tau_pf: |A_2| not bounded below")
    g = []
    for j in range(1, 7):
        Aj = A_[j + 2].absu() if j + 2 < len(A_) else arb(0)
        Bj = B_[j + 1].absu() if j + 1 < len(B_) else arb(0)
        Cj = C_[j].absu() if j < len(C_) else arb(0)
        gj = (Aj * (1 + Rr / Nn) ** 2 / (1 - 2 * Rr / Nn) + Bj * (Nn + Rr) / (Nn * (Nn - 2 * Rr))
              + Cj / (Nn * (Nn - 2 * Rr))) / A2l
        g.append(up(gj))
    q = arb(_solve_q(g))
    if not (sum((gj * q ** (-(j + 1)) for j, gj in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("frob_tau_pf majorant root not verified")
    K = max((a[m].absu() / q ** m for m in range(M - 6, M)), key=lambda x: x.upper())
    s0, s1 = _geom_tails(up(K), q, abs(ta), M, Rr)
    S = S + AF(0, 0, 0, s0, S.Rw, S.Rh)
    Sd = Sd + AF(0, 0, 0, s1, S.Rw, S.Rh)
    return S, Sd, {"terms": M}


def y_exit_pf(pt: Point, side: int, wbox, h2_lo: Fraction, h2_hi: Fraction, tau: Fraction):
    """lowfreq.y_exit with frob_tau_pf (y = Phi_x/Phi does not depend on the normalisation)."""
    from . import lowfreq as LF
    w, h2 = LF.params_af(pt, wbox, h2_lo, h2_hi)
    kap = pt.kp if side > 0 else pt.kc
    rho = w * (-acb(0, 1) / (2 * kap))
    A_, B_, C_, Dc, (r0, r1) = LF.tau_polys(pt, side, w, h2)
    tq = fmpq(tau.numerator, tau.denominator)
    S, Sd, info = frob_tau_pf(A_, B_, C_, rho, tq)
    ta = arb(tq)
    Dt = Dc[1] * ta
    for k in range(2, 5):
        Dt = Dt + Dc[k] * ta ** k
    r = r0 + r1 * ta
    y = h2 * Dt / (r * r * r) + (Dt / (r * r)) * (Sd / S) * side
    return y, info


def make_lchecker(pt: Point, W_bound: Fraction):
    """S4's LChecker with the pole-free series in every box, and S6's pole condition: the only point of
    E_+ u E_c a band box can meet is omega0 = -i kappa_+, where R~ is used (W6')."""
    from . import lowfreq as LF

    class LCheckerPF(LF.LChecker):
        def side_check(self, side, wbox, tau):
            worst = None
            terms = 0
            for a, b in LF.h2_pieces():
                y, info = y_exit_pf(self.pt, side, wbox, a, b, tau)
                v = y.re_lower() if side > 0 else y.re_upper()
                terms = max(terms, info["terms"])
                if worst is None or (v < worst if side > 0 else v > worst):
                    worst = v
            return worst, terms

        def box(self, wbox, exit_scale=None):
            rec = super().box(wbox, exit_scale)
            km = self.pt.km
            im_lo = A(wbox[2]) * km
            # every box lies above -2 kappa_+ (the second point of E_+) and -kappa_c (E_c's first point)
            e_ok = bool(im_lo > -2 * self.pt.kp) and bool(im_lo > -self.pt.kc)
            rec.pop("poles_avoided", None)
            rec["e_points_only_omega0"] = e_ok
            rec["series"] = "pole-free R~ = (1 + 2 rho) R (W6')"
            ok = all(rec.get(s, {}).get("y_sign_ok") and rec.get(s, {}).get("exit_rule_ok") for s in ("plus", "c"))
            rec["pass"] = bool(ok and e_ok)
            return rec

    return LCheckerPF(pt, W_bound)


# ---------------------------------------------------------------------------------------------
# float evaluator of the T/A bound (design selection only; part1_T formulas)
# ---------------------------------------------------------------------------------------------
class FloatTA:
    """F(s) = E0/h0 + E1 + h0 min(0, E2) on S2's contour, in floats (part1_T.TChecker.parts)."""

    def __init__(self, pt: Point, h=1.0 / 256, s_lo=-64.0, s_hi=32.0):
        import numpy as np
        from .m3 import float_contour
        self.np = np
        s, r = float_contour(pt, s_lo, s_hi, h)
        self.s, self.r = s, r
        lam3 = float((pt.Lam / 3).mid())
        Q2 = float(pt.Q2.mid())
        D = -lam3 * r ** 4 + r ** 2 - 2 * r + Q2
        Dp = -4 * lam3 * r ** 3 + 2 * r - 2
        F = D / r ** 2
        Fr = (Dp * r - 2 * D) / r ** 3
        self.U = F / r ** 2
        self.W = F * Fr / r - self.U / 4
        self.OmP = float(pt.OmegaP.mid())
        self.km = float(pt.km.mid())
        self.astar = float((2 * pt.OmegaP * pt.gamma0).mid())

    def F(self, d: dict, zl, zh, nl, h0):
        np = self.np
        s = self.s
        s1, phinf, phi0, eta = (float(Fraction(d[k])) for k in ("s1", "phi_inf", "phi0", "eta"))
        phi = phi0 + (phinf - phi0) * s ** 4 / (s1 ** 4 + s ** 4)
        ds = s[1] - s[0]
        d1 = np.gradient(phi, ds)
        c = np.cos(phi) - eta * np.abs(d1)
        cp = np.gradient(c, ds)
        if d["eps_g"] is None:
            G = np.zeros_like(s)
        else:
            k = (1 - float(Fraction(d["eps_g"]))) * float(Fraction(d["a_d"]))
            tc, td = (float(Fraction(x)) for x in d["taper"])
            a = np.abs(s)
            th = np.where(a <= tc, 1.0, np.where(a >= td, 0.0, (1 + np.cos(np.pi * (a - tc) / (td - tc))) / 2))
            G = k * s * th
        Gp = np.gradient(G, ds)
        rot = np.sin(phi) + 1j * np.cos(phi)
        E0 = np.real(rot * (self.U - zl * zl)) - c * G * G
        E1 = 2 * zh * nl * np.cos(phi) + cp * G + c * Gp
        E2 = np.real(rot * (self.W + nl * nl)) - np.abs(d1) / (4 * eta)
        Fv = E0 / h0 + E1 + h0 * np.minimum(0, E2)
        i = int(np.argmin(Fv))
        return {"F_min_over_astar": float(Fv[i] / self.astar), "s_at_min": float(s[i]),
                "E0_min": float(E0.min()), "c_min": float(c.min())}

    def box(self, d: dict, zb, depth, h0):
        zl, zh = float(zb[0]) * self.OmP, float(zb[1]) * self.OmP
        nl = -float(depth) * self.km
        return self.F(d, zl, zh, nl, float(h0))
