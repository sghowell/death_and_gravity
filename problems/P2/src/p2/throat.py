"""S4 regime S (Omega_S <= Re omega <= zeta_SU/h): half-line M3 certificates at the two throats.
Pre-registration section 3.2, Amendments 1, 3 and 5.  Proof: notes/s4-proof.md (W1-W3, W5).

Normalised form.  With zeta = h Re(omega), eps = 1/Re(omega), nu = Im(omega) in [-d_delta, 0]:
    V - omega^2 = (Re omega)^2 [U/zeta^2 - 1 + 2 i eps |nu| + eps^2 (nu^2 + Wt)]      (exact).
Contour (side s = +1 at r_+, -1 at r_c): x(t) = X_a + s (Y(t) - t_a)/k, t <= t_a, with
    Y = t + i eta(t),  eta = t (t <= -w),  -w m0(-t/w) (|t| <= w),  0 (t >= w),  t_a = w + y_a,
m0 the C^3 smoothed max(0, u) on [-1, 1] (m0'' = (15/16)(1 - u^2)^2), so x(t_a) = X_a is real, x' is
real there, and for t <= -w the contour is a ray of angle pi/4 (toward the horizon).  Phase
phi = phi_inf for t <= p_c - p_w, phi_inf (1 - S5(v)) on the transition (quintic smoothstep), 0 after
(p_c + p_w <= t_a).  Half-line coercivity (W3) with the pointwise Young weight cos(phi)/kappa_tot
(kappa_tot = |Y''/Y' - i phi'|; kinetic coefficient 0):
    Re(e^{-i phi(t_a)} u' conj u)(t_a) >= (Re omega/k)^2 int q^ |u|^2,
    q^ = Re(e^{-i phi} Y'^2 [U/zeta^2 - 1 + 2 i eps |nu| + eps^2 (nu^2 + Wt)]) - (k eps)^2 |kappa|^2/(4 cos phi).
q^ > 0 on (-inf, t_a] gives Re y_+(X_a) > 0 (side +; x' = 1/k > 0 at t_a) resp. Re y_c(X_a) < 0
(side c; x' = -1/k).  The bulk (Riccati invariance) needs U(r_a)/zeta^2 - 1 > eps_max^2 W_bound and
X_a < x_P (side +) resp. X_a > x_P (side c).

Per zeta-box [zl, zh]: the design (grid of Amendment 5, float maximin), X_a and k (dyadics, from the
geometric-mean zeta) are fixed; q^ is affine in 1/zeta^2 and quadratic in eps, |nu| with separable
interval bounds.  Cells: mean-value form in t (t-derivatives through acb_series).  r(t) by validated
Taylor integration of dr/dt = x'(t) F(r) from the real exit (x' piecewise polynomial).  Tail
t <= t_L: S2's horizon disc bound.  Closing box zeta in (0, zeta_1]: the near-horizon map
rho = exp(+-2 kappa_h (x - c_h)) = rho_T e^{2Y} (k = kappa_h, x_T(zeta) exact), no integration.
"""
from __future__ import annotations

import itertools
import json
import math
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import acb, acb_series, arb, ctx

from .arbutil import Timer, git_commit, precision
from .m3 import dist_bounds, fa, real_point
from .m3_contour import F_of, apriori_box, rstar
from .params import Point

PREC = 128
H0 = Fraction(2, 41)
OMEGA_S_KM = Fraction(6)
ZETA_SU = Fraction(9, 10)
ZETA_1 = Fraction(1, 50)
ZBOX = Fraction(1, 50)
ZBOX_MIN = Fraction(1, 400)
EPS_LEVELS = 4
T_LEFT = Fraction(-8)
H_T = Fraction(1, 16)
N_TAYLOR = 20
MAX_DEPTH = 8
T_CLOSE_MAX = Fraction(21, 10)
D_R = {1: Fraction(3, 50), -1: Fraction(1, 1)}      # near-horizon r-discs (closing box)
GRID_W = [Fraction(3, 20), Fraction(3, 10), Fraction(3, 5), Fraction(1), Fraction(2), Fraction(3), Fraction(4)]
GRID_YA = [Fraction(1, 20), Fraction(1, 10), Fraction(1, 5), Fraction(2, 5), Fraction(4, 5)]
GRID_PHI = [4, 3]                                    # phi_inf = -pi/4, -pi/3
GRID_PC = [Fraction(-1), Fraction(-1, 2), Fraction(0)]
GRID_PW = [Fraction(1, 4), Fraction(1, 2), Fraction(1)]
ROOT = Path(__file__).resolve().parents[2]
CERT_DIR = ROOT / "certificates" / "part1"
RES_DIR = ROOT / "results" / "part1"
I = acb(0, 1)


def lohi(x, digits=8) -> str:
    """Outward-rounded endpoints of an arb as a readable string (Arb's short rendering can hide the sign)."""
    return f"[{float(arb(x.lower())):.{digits}g}, {float(arb(x.upper())):.{digits}g}]"


def pval(co, x):
    """Polynomial (ascending exact-rational coefficients co) at an arb x, in centred form: the
    coefficients are shifted to the midpoint m of x and summed against (x - m)^k (no dependency loss)."""
    if x.rad() == 0:
        acc = arb(0)
        for c in reversed(co):
            acc = acc * x + fa(c)
        return acc
    m = arb(x.mid())
    dlt = arb(0, x.rad())
    n = len(co)
    acc = arb(0)
    dp = arb(1)
    for k in range(n):
        dk = arb(0)
        mp = arb(1)
        for j in range(k, n):
            dk += math.comb(j, k) * fa(co[j]) * mp
            mp = mp * m
        acc += dk * dp
        dp = dp * dlt
    return acc


M0 = [Fraction(5, 32), Fraction(16, 32), Fraction(15, 32), 0, Fraction(-5, 32), 0, Fraction(1, 32)]
M1 = [Fraction(8, 16), Fraction(15, 16), 0, Fraction(-10, 16), 0, Fraction(3, 16)]
M2 = [Fraction(15, 16), 0, Fraction(-30, 16), 0, Fraction(15, 16)]
M3 = [0, Fraction(-15, 4), 0, Fraction(15, 4)]
S5_0 = [0, 0, 0, 10, -15, 6]
S5_1 = [0, 0, 30, -60, 30]
S5_2 = [0, 60, -180, 120]


# ---------------------------------------------------------------------------------------------
# the design family (exact rationals; phi_inf = -pi/n)
# ---------------------------------------------------------------------------------------------
class SDesign:
    def __init__(self, w, ya, phn, pc, pw):
        self.w, self.ya, self.phn, self.pc, self.pw = (Fraction(w), Fraction(ya), int(phn), Fraction(pc), Fraction(pw))
        self.ta = self.w + self.ya
        self.p0, self.p1 = self.pc - self.pw, self.pc + self.pw

    def key(self):
        return [str(self.w), str(self.ya), f"-pi/{self.phn}", str(self.pc), str(self.pw)]

    @staticmethod
    def from_key(k):
        return SDesign(Fraction(k[0]), Fraction(k[1]), int(k[2].split("/")[1]), Fraction(k[3]), Fraction(k[4]))

    def admissible(self):
        return self.p1 <= self.ta

    def breakpoints(self):
        return sorted({-self.w, self.w, self.p0, self.p1})

    def phinf(self):
        return -arb.pi() / self.phn

    # ---- contour: Y', Y'', Y''' on a t-interval inside one piece (arb t) ----
    def piece_bend(self, a: Fraction, b: Fraction) -> int:
        if b <= -self.w:
            return -1
        if a >= self.w:
            return 1
        if a >= -self.w and b <= self.w:
            return 0
        raise ValueError("cell straddles a bend breakpoint")

    def piece_phase(self, a: Fraction, b: Fraction) -> int:
        if b <= self.p0:
            return -1
        if a >= self.p1:
            return 1
        if a >= self.p0 and b <= self.p1:
            return 0
        raise ValueError("cell straddles a phase breakpoint")

    def eta_derivs(self, t, piece):
        """(eta', eta'', eta''') on the piece."""
        if piece < 0:
            return arb(1), arb(0), arb(0)
        if piece > 0:
            return arb(0), arb(0), arb(0)
        w = fa(self.w)
        u = -t / w
        return pval(M1, u), -pval(M2, u) / w, pval(M3, u) / (w * w)

    def eta_value(self, t, piece):
        if piece < 0:
            return t
        if piece > 0:
            return arb(0)
        w = fa(self.w)
        return -w * pval(M0, -t / w)

    def phase(self, t, piece):
        """(phi, phi', phi'') on the piece."""
        phinf = self.phinf()
        if piece < 0:
            return phinf, arb(0), arb(0)
        if piece > 0:
            return arb(0), arb(0), arb(0)
        pw2 = 2 * fa(self.pw)
        v = (t - fa(self.p0)) / pw2
        return phinf * (1 - pval(S5_0, v)), -phinf * pval(S5_1, v) / pw2, -phinf * pval(S5_2, v) / (pw2 * pw2)

    def xprime_poly(self, t0, side: int, k, piece):
        """Coefficients (acb) of x'(t0 + s) = s_side (1 + i eta'(t0 + s))/k as a polynomial in s
        (exact: eta' is a quintic on the bend); t0 may be an interval."""
        if piece < 0:
            return [acb(side, side) / k]
        if piece > 0:
            return [acb(side) / k]
        w = fa(self.w)
        u = acb_series([acb(-t0 / w), acb(-1 / w)], prec=6)
        m1 = (15 * u - 10 * u ** 3 + 3 * u ** 5 + 8) / 16
        co = m1.coeffs()
        co = [co[j] if j < len(co) else acb(0) for j in range(6)]
        return [side * ((1 if j == 0 else 0) + I * co[j]) / k for j in range(6)]


# ---------------------------------------------------------------------------------------------
# potential as series in the local parameter
# ---------------------------------------------------------------------------------------------
def pot_series(pt: Point, r):
    """(U, Wt) for an acb_series r (or acb)."""
    lam3 = pt.Lam / 3
    D = -lam3 * (r - pt.roots[0]) * (r - pt.roots[1]) * (r - pt.roots[2]) * (r - pt.roots[3])
    Dp = -4 * lam3 * r * r * r + 2 * r - 2
    r2 = r * r
    F = D / r2
    Fr = (Dp * r - 2 * D) / (r2 * r)
    U = F / r2
    W = F * Fr / r - U / 4
    return U, W


def _re(x):
    return [c.real for c in x.coeffs()] + [arb(0)] * (2 - len(x.coeffs()))


def _im(x):
    return [c.imag for c in x.coeffs()] + [arb(0)] * (2 - len(x.coeffs()))


def parts(pt: Point, D: SDesign, t, rr, rprime, side: int, k, pb: int, pp: int, model=None):
    """[value, d/dt] of P1 = Re(G U), P0 = -Re G, P2 = -2 Im G, P3 = Re G, P4 = Re(G Wt), K, |Y'|^2
    with G = e^{-i phi} Y'^2, K = |Y''/Y' - i phi'|^2/(4 cos phi); t a point or interval, rr (acb) the
    enclosure of r there, rprime of dr/dt (ignored when model is given: closing box, where model =
    (e2Y_series) replaces U by the exponential throat and P4 by a bound)."""
    e1, e2, e3 = D.eta_derivs(t, pb)
    ph, ph1, ph2 = D.phase(t, pp)
    Yp = acb_series([acb(1, e1), acb(0, e2)], prec=2)
    Ypp = acb_series([acb(0, e2), acb(0, e3)], prec=2)
    phs = acb_series([acb(ph), acb(ph1)], prec=2)
    ph1s = acb_series([acb(ph1), acb(ph2)], prec=2)
    G = (-I * phs).exp() * Yp * Yp
    out = {"P0": [-x for x in _re(G)], "P2": [-2 * x for x in _im(G)], "P3": _re(G)}
    if model is None:
        rs = acb_series([rr, rprime], prec=2)
        U, W = pot_series(pt, rs)
        out["P1"] = _re(G * U)
        out["P4"] = _re(G * W)
    else:
        out["P1"] = _re(G * model)          # Re(G e^{2Y})
    kap = Ypp / Yp - I * ph1s
    kc = kap.coeffs() + [acb(0)] * 2
    k0, k1 = kc[0], kc[1]
    m2 = [abs(k0) ** 2, 2 * (k0.conjugate() * k1).real]
    c0, c1 = ph.cos(), -ph.sin() * ph1
    out["K"] = [m2[0] / (4 * c0), (m2[1] * c0 - m2[0] * c1) / (4 * c0 * c0)]
    yp2 = [abs(acb(1, e1)) ** 2, 2 * e1 * e2]
    out["Yp2"] = yp2
    out["cosphi"] = c0
    return out


# ---------------------------------------------------------------------------------------------
# validated contour
# ---------------------------------------------------------------------------------------------
def taylor_coeffs_P(pt: Point, r0, P, N):
    """Taylor coefficients r_0..r_N of r' = P(s) F(r) at s = 0 (P a list of acb coefficients)."""
    lam3 = pt.Lam / 3
    Q2 = acb(pt.Q2)
    r = [acb(r0)]
    u = [1 / r[0]]
    Fk = []
    for k in range(N):
        r2k = sum((r[j] * r[k - j] for j in range(k + 1)), acb(0))
        u2k = sum((u[j] * u[k - j] for j in range(k + 1)), acb(0))
        Fk.append(-lam3 * r2k - 2 * u[k] + Q2 * u2k + (1 if k == 0 else 0))
        acc = acb(0)
        for j in range(min(k, len(P) - 1) + 1):
            acc += P[j] * Fk[k - j]
        r.append(acc / (k + 1))
        n = k + 1
        u.append(-sum((r[j] * u[n - j] for j in range(1, n + 1)), acb(0)) / r[0])
    return r


def _poly_at(P, s):
    acc = acb(0)
    for c in reversed(P):
        acc = acc * s + c
    return acc


class CStep:
    __slots__ = ("D0", "J", "N", "co", "h", "rem", "t0")

    def __init__(self, t0, h, co, rem, D0, J, N):
        self.t0, self.h, self.co, self.rem, self.D0, self.J, self.N = t0, h, co, rem, D0, J, N

    def eval_mid(self, s):
        acc = acb(0)
        for c in reversed(self.co):
            acc = acc * s + c
        sn = s ** (self.N + 1)
        return acc + acb(self.rem.real * sn, self.rem.imag * sn)

    def eval(self, s):
        return self.eval_mid(s) + self.J * self.D0


class CurveContour:
    """Validated r(t) along x(t) = X_a + s (Y(t) - t_a)/k for t in [t_L, t_a], integrated backward from
    the real exit r_a = r(X_a) (interval Newton)."""

    def __init__(self, pt: Point, side: int, D: SDesign, Xa: Fraction, k: Fraction, tL=T_LEFT, h=H_T, N=N_TAYLOR,
                 tol=None, hmin=Fraction(1, 2 ** 22)):
        tol = arb("1e-30") if tol is None else tol
        self.pt, self.side, self.D, self.N = pt, side, D, N
        self.Xa, self.k = Xa, k
        self.ka = fa(k)
        self.ra = real_point(pt, fa(Xa))
        bps = sorted({tL} | {x for x in (-D.w, D.w) if tL < x < D.ta})
        r0 = acb(self.ra)
        steps = []
        t = D.ta
        while t > tL:
            nb = max(x for x in bps if x < t)
            hh = min(h, t - nb)
            while True:
                t_lo = t - hh
                pb = D.piece_bend(t_lo, t)
                hs = fa(t_lo - t)
                t0 = fa(t)
                step_iv = (fa(t_lo) + fa(t)) / 2 + arb(0, 1) * (fa(hh) / 2)
                Pxi = D.xprime_poly(step_iv, side, self.ka, pb)
                Pstep = _poly_at(Pxi, acb(0))
                B = apriori_box(pt, r0, Pstep, hs)
                if B is not None:
                    rem = taylor_coeffs_P(pt, B, Pxi, N + 1)[N + 1]
                    if abs(rem).upper() * abs(hs) ** (N + 1) < tol:
                        break
                if hh / 2 < hmin:
                    raise ArithmeticError(f"contour step below {hmin} at t = {t}")
                hh = hh / 2
            P0 = D.xprime_poly(t0, side, self.ka, pb)
            mid = acb(r0.real.mid(), r0.imag.mid())
            D0 = r0 - mid
            co = taylor_coeffs_P(pt, mid, P0, N)
            F0 = F_of(pt, r0)
            J = F_of(pt, B) / F0
            st = CStep(t0, hs, co, rem, D0, J, N)
            Rh = st.eval_mid(hs) + J * D0
            Jh = F_of(pt, Rh) / F0
            r0 = st.eval_mid(hs) + Jh * D0
            steps.append((t_lo, t, st, pb))
            t = t_lo
        self.steps = list(reversed(steps))
        self.nodes = [s_[0] for s_ in self.steps] + [D.ta]
        self.r_L = r0
        self.tL = tL

    def find(self, a: Fraction, b: Fraction):
        lo_, hi_ = 0, len(self.steps) - 1
        while lo_ <= hi_:
            m = (lo_ + hi_) // 2
            t_lo, t_hi, st, pb = self.steps[m]
            if b <= t_lo:
                hi_ = m - 1
            elif a >= t_hi:
                lo_ = m + 1
            else:
                if a >= t_lo and b <= t_hi:
                    return st, pb
                raise ValueError("interval straddles a contour step")
        raise IndexError((a, b))

    def enclose(self, a: Fraction, b: Fraction):
        st, pb = self.find(a, b)
        s_lo, s_hi = fa(a) - st.t0, fa(b) - st.t0
        S = (s_lo + s_hi) / 2 + arb(0, 1) * ((s_hi - s_lo) / 2)
        return st.eval(S), pb

    def at(self, t: Fraction):
        st, pb = self.find(t, t)
        return st.eval(fa(t) - st.t0), pb

    def xprime(self, t, pb):
        """x'(t) (acb) at t (arb point or interval) inside a bend piece."""
        e1, _, _ = self.D.eta_derivs(t, pb)
        return self.side * acb(1, e1) / self.ka


# ---------------------------------------------------------------------------------------------
# the cell check
# ---------------------------------------------------------------------------------------------
class Params:
    """Interval factors of q^ for a (zeta, eps) box and nu in [-d_delta, 0]."""

    def __init__(self, pt: Point, zl, zh, el: Fraction, eh: Fraction, k=None):
        self.Z = (1 / (zh * zh)) .union(1 / (zl * zl)) if zl is not None else None
        d = pt.beta
        e_hi = fa(eh)
        e_lo = fa(el)
        self.EN = arb(0).union(e_hi * d)
        self.E2N2 = arb(0).union((e_hi * d) ** 2)
        self.E2 = (e_lo * e_lo).union(e_hi * e_hi)
        self.k2 = None if k is None else k * k

    def combine(self, p, idx: int, k2):
        q = p["P0"][idx] + self.EN * p["P2"][idx] + self.E2N2 * p["P3"][idx] - k2 * self.E2 * p["K"][idx]
        if self.Z is not None:
            q = q + self.Z * p["P1"][idx]
        else:
            q = q + p["P1"][idx]
        if "P4" in p:
            q = q + self.E2 * p["P4"][idx]
        return q


class SChecker:
    def __init__(self, pt: Point, side: int, D: SDesign, C: CurveContour, prm: Params):
        self.pt, self.side, self.D, self.C, self.prm = pt, side, D, C, prm
        self.k2 = C.ka * C.ka

    def cell(self, a: Fraction, b: Fraction):
        D, C = self.D, self.C
        pb = D.piece_bend(a, b)
        pp = D.piece_phase(a, b)
        a_, b_ = fa(a), fa(b)
        w = b_ - a_
        tm = (a_ + b_) / 2
        cell = tm + arb(0, 1) * (w / 2)
        rm, _ = C.at((a + b) / 2)
        R, _ = C.enclose(a, b)
        pm = parts(self.pt, D, tm, rm, C.xprime(tm, pb) * F_of(self.pt, rm), self.side, C.ka, pb, pp)
        pc = parts(self.pt, D, cell, R, C.xprime(cell, pb) * F_of(self.pt, R), self.side, C.ka, pb, pp)
        qm = self.prm.combine(pm, 0, self.k2)
        dq = self.prm.combine(pc, 1, self.k2)
        return arb((qm.lower() - abs(dq).upper() * (w / 2)).lower()), pc["cosphi"]

    def cells(self, target, max_depth=MAX_DEPTH):
        D = self.D
        nodes = set(self.C.nodes) | {x for x in D.breakpoints() if self.C.tL < x < D.ta}
        nodes = sorted(nodes)
        mu, worst, n_eval, fails, depths = None, None, 0, [], {}
        for a, b in itertools.pairwise(nodes):
            stack = [(a, b, 0)]
            while stack:
                a_, b_, d = stack.pop()
                low, cphi = self.cell(a_, b_)
                n_eval += 1
                if (low < target) and d < max_depth:
                    m = (a_ + b_) / 2
                    stack.append((m, b_, d + 1))
                    stack.append((a_, m, d + 1))
                    continue
                depths[d] = depths.get(d, 0) + 1
                if not (low > 0 and cphi > 0):
                    fails.append([str(a_), str(b_), low.str(5)])
                if mu is None or low < mu:
                    mu, worst = low, (float(a_), float(b_))
        return mu, worst, n_eval, fails, depths


# ---------------------------------------------------------------------------------------------
# tails, exit, W1
# ---------------------------------------------------------------------------------------------
def tail_disc(pt: Point, side: int, X_L, r_L, d0):
    """sup |r - r_h| over the contour tail (Re x <= X_L for side +, >= X_L for side c), or None."""
    rh = pt.rp if side > 0 else pt.rc
    Ah = pt.Acoef[2] if side > 0 else pt.Acoef[3]
    others = [(Ai, ri) for Ai, ri in zip(pt.Acoef, pt.roots) if not ri.overlaps(rh)]
    disc = acb(rh + arb(0, 1) * d0, arb(0, 1) * d0)
    Smin, Smax = arb(0), arb(0)
    for Ai, ri in others:
        lo_, hi_ = dist_bounds(disc, ri)
        if not (lo_ > 0):
            return None
        t1, t2 = Ai * lo_.log(), Ai * hi_.log()
        Smin = Smin + arb(t1.lower()).min(arb(t2.lower()))
        Smax = Smax + arb(t1.upper()).max(arb(t2.upper()))
    if not (abs(r_L - rh) < d0):
        return None
    if side > 0:
        d = ((X_L - arb(Smin.lower())) / Ah).exp()
    else:
        d = ((X_L - arb(Smax.upper())) / Ah).exp()
    d = arb(d.upper())
    return d if d < d0 else None


def tail_bound(pt: Point, side: int, D: SDesign, C: CurveContour, prm: Params, k2):
    """Lower bound of q^ on t <= t_L (straight ray, phi = phi_inf, K = 0)."""
    X_L = fa(C.Xa) + side * (fa(C.tL) - fa(D.ta)) / C.ka
    d = None
    for d0 in ("1e-8", "1e-6", "1e-4", "1e-3", "1e-2", "5e-2"):
        d = tail_disc(pt, side, X_L, C.r_L, arb(d0))
        if d is not None:
            break
    if d is None:
        return None, None
    rh = pt.rp if side > 0 else pt.rc
    disc = acb(rh + arb(0, 1) * d, arb(0, 1) * d)
    U, W = pot_series(pt, disc)
    Uu, Wu = abs(U).upper(), abs(W).upper()
    phinf = D.phinf()
    # G = 2 i e^{-i phi_inf}: P0 = -2 sin phi_inf, P2 = -4 cos phi_inf, P3 = 2 sin phi_inf, |G| = 2
    q = (-2 * phinf.sin() + prm.EN * (-4 * phinf.cos()) + prm.E2N2 * (2 * phinf.sin())
         - 2 * Uu * (prm.Z.upper() if prm.Z is not None else 1) - 2 * Wu * prm.E2.upper())
    return arb(q.lower()), d


def exit_checks(pt: Point, side: int, D: SDesign, C: CurveContour, zh, eps_max: Fraction, W_bound: Fraction):
    ra = C.ra
    U = pt.Delta(ra) / ra ** 4
    margin = U / (zh * zh) - 1 - fa(eps_max) ** 2 * fa(W_bound)
    order = bool(ra < pt.rP) if side > 0 else bool(ra > pt.rP)
    return {"r_a": ra.str(12), "Ut_minus_1_minus_epsmax2_Wbound": lohi(margin),
            "exit_forbidden_ok": bool(margin > 0), "exit_order_ok": order,
            "phi_zero_at_exit": bool(D.p1 <= D.ta), "eta_zero_at_exit": bool(D.w <= D.ta)}


def w1_check(pt: Point, side: int, D: SDesign, C: CurveContour):
    """Amendment 1: side + lies in the strip Re x < Re r*(0); side c crosses Re x = Re r*(0) and Re x = 0
    only within S2's checked segments (|Im x| <= 0.1413 resp. 2.438)."""
    x0re = sum((Ai * abs(ri).log() for Ai, ri in zip(pt.Acoef, pt.roots)), arb(0))
    Xa = fa(C.Xa)
    if side > 0:
        ok = bool(Xa < x0re)
        return {"side": "+", "X_a": str(C.Xa), "Re_rstar_0": x0re.str(10), "left_strip": ok, "pass": ok}
    out = {"side": "c", "X_a": str(C.Xa), "Re_rstar_0": x0re.str(10), "crossings": []}
    ok = True
    for X, lim in ((x0re, arb("0.1413")), (arb(0), arb("2.438"))):
        if X < Xa:
            continue                     # the contour (Re x >= X_a) never meets this line
        if not (X > Xa):
            ok = False
            out["crossings"].append({"line": X.str(8), "undecided": True})
            continue
        # Re x(t) = X_a - (t - t_a)/k = X  =>  t* = t_a - k (X - X_a)
        ts = fa(D.ta) - C.ka * (X - Xa)
        if ts < fa(-D.w):
            eta = ts
        elif ts > fa(D.w):
            eta = arb(0)
        else:
            eta = D.eta_value(ts, 0)
        imx = abs(eta / C.ka)
        good = bool(imx < lim)
        ok = ok and good
        out["crossings"].append({"line": X.str(8), "t_star": ts.str(8), "abs_Im_x": imx.str(8), "limit": lim.str(6), "ok": good})
    out["pass"] = ok
    return out


# ---------------------------------------------------------------------------------------------
# float design (numpy): the grid maximin of Amendment 5
# ---------------------------------------------------------------------------------------------
class FGeo:
    def __init__(self, pt: Point):
        self.roots = [float(x.mid()) for x in pt.roots]
        self.lam3 = float((pt.Lam / 3).mid())
        self.Q2 = float(pt.Q2.mid())
        self.A = [float(a.mid()) for a in pt.Acoef]
        self.rp, self.rc, self.rP = float(pt.rp.mid()), float(pt.rc.mid()), float(pt.rP.mid())
        self.km, self.kp, self.kc = float(pt.km.mid()), float(pt.kp.mid()), float(pt.kc.mid())
        self.OmP = float(pt.OmegaP.mid())
        self.beta = float(pt.beta.mid())

    def Delta(self, r):
        return -self.lam3 * r ** 4 + r ** 2 - 2 * r + self.Q2

    def F(self, r):
        return self.Delta(r) / r ** 2

    def U(self, r):
        return self.Delta(r) / r ** 4

    def Wt(self, r):
        D = self.Delta(r)
        Dp = -4 * self.lam3 * r ** 3 + 2 * r - 2
        F = D / r ** 2
        Fr = (Dp * r - 2 * D) / r ** 3
        return F * Fr / r - F / (4 * r * r)

    def Ux(self, r):
        D = self.Delta(r)
        Dp = -4 * self.lam3 * r ** 3 + 2 * r - 2
        return self.F(r) * (Dp * r - 4 * D) / r ** 5

    def rstar(self, r):
        return sum(a * math.log(abs(r - x)) for a, x in zip(self.A, self.roots))

    def turning(self, zeta, side):
        lo_, hi_ = (self.rp + 1e-13, self.rP) if side > 0 else (self.rP, self.rc - 1e-13)
        a, b = lo_, hi_
        for _ in range(200):
            m = 0.5 * (a + b)
            if (self.U(a) - zeta ** 2 < 0) == (self.U(m) - zeta ** 2 < 0):
                a = m
            else:
                b = m
        return 0.5 * (a + b)


def _fbend(t, w):
    u = np.clip(-t / w, -1, 1)
    inside = np.abs(t) <= w
    deta = np.where(t <= -w, 1.0, np.where(t >= w, 0.0, (15 * u - 10 * u ** 3 + 3 * u ** 5 + 8) / 16))
    d2eta = np.where(inside, -(15 / 16) * (1 - u ** 2) ** 2 / w, 0.0)
    return deta, d2eta


def _fphase(t, phinf, p0, pw):
    v = np.clip((t - p0) / (2 * pw), 0, 1)
    inside = (t >= p0) & (t <= p0 + 2 * pw)
    S = 10 * v ** 3 - 15 * v ** 4 + 6 * v ** 5
    dS = np.where(inside, 30 * v ** 2 * (1 - v) ** 2, 0.0)
    return phinf * (1 - S), -phinf * dS / (2 * pw)


_FCACHE: dict = {}


def _fcontour(g: FGeo, side: int, w: Fraction, ya: Fraction, zl: float, zh: float, closing: bool, tL: float, dt: float):
    """Float contour data (independent of the phase): t, Y, Y', Y'', U(r(t)), Wt(r(t)), k, r_a."""
    key = (side, w, ya, zl, zh, closing, tL, dt)
    if key in _FCACHE:
        return _FCACHE[key]
    ta = float(w + ya)
    wf = float(w)
    n = round((ta - tL) / dt)
    t = np.linspace(tL, ta, n + 1)
    deta, d2eta = _fbend(t, wf)
    u = np.clip(-t / wf, -1, 1)
    eta = np.where(t <= -wf, t, np.where(t >= wf, 0.0, -wf * (15 * u ** 2 - 5 * u ** 4 + u ** 6 + 16 * u + 5) / 32))
    Y = t + 1j * eta
    Yp = 1 + 1j * deta
    Ypp = 1j * d2eta
    if closing:
        k = g.kp if side > 0 else g.kc
        out = (t, Y, Yp, Ypp, None, 0 * t, k, None)
        _FCACHE[key] = out
        return out
    zref = math.sqrt(zl * zh)
    rT = g.turning(zref, side)
    k = abs(g.Ux(rT)) / (2 * zref ** 2)
    xT = g.rstar(rT)
    xa = xT + side * ta / k
    r = complex(rT)
    nn = 400
    for _ in range(nn):
        dx = (xa - xT) / nn
        k1 = g.F(r); k2 = g.F(r + dx / 2 * k1); k3 = g.F(r + dx / 2 * k2); k4 = g.F(r + dx * k3)
        r = r + dx / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    ra = r.real
    if (side > 0 and not (g.rp < ra < g.rP)) or (side < 0 and not (g.rP < ra < g.rc)):
        _FCACHE[key] = None
        return None
    rr = np.empty(len(t), complex)
    rr[-1] = ra
    rcur = complex(ra)
    xp = side * Yp / k
    for i in range(len(t) - 1, 0, -1):
        h = t[i - 1] - t[i]
        xm = 0.5 * (xp[i] + xp[i - 1])
        k1 = xp[i] * g.F(rcur); k2 = xm * g.F(rcur + h / 2 * k1); k3 = xm * g.F(rcur + h / 2 * k2)
        k4 = xp[i - 1] * g.F(rcur + h * k3)
        rcur = rcur + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        rr[i - 1] = rcur
    out = (t, Y, Yp, Ypp, g.U(rr), g.Wt(rr), k, ra)
    _FCACHE[key] = out
    return out


def float_margin(g: FGeo, side: int, D: SDesign, zl: float, zh: float, eps_max: float, W_bound: float,
                 closing=False, tL=-8.0, dt=0.01):
    """Float min_t q^ over the box (zeta endpoints, eps and |nu| samples) and the exit data, or None."""
    if not D.admissible():
        return None
    fc = _fcontour(g, side, D.w, D.ya, zl, zh, closing, tL, dt)
    if fc is None:
        return None
    t, Y, Yp, Ypp, Ur, Wr, k, ra = fc
    if not closing and not (g.U(ra) / zh ** 2 - 1 > eps_max ** 2 * W_bound):
        return None
    phi, dphi = _fphase(t, -math.pi / D.phn, float(D.p0), float(D.pw))
    G = np.exp(-1j * phi) * Yp ** 2
    kt = np.abs(Ypp / Yp - 1j * dphi) ** 2 / (4 * np.cos(phi))
    Ut = [np.exp(2 * Y)] if closing else [Ur / zl ** 2, Ur / zh ** 2]
    worst = np.full(len(t), np.inf)
    for Uz in Ut:
        base = np.real(G * (Uz - 1))
        for e in np.linspace(0, eps_max, 5):
            for nu in (0.0, g.beta / 2, g.beta):
                val = base + np.real(G * (2j * e * nu + e * e * (nu * nu + Wr))) - (k * e) ** 2 * kt
                worst = np.minimum(worst, val)
    m = float(worst.min())
    return {"margin": m, "k": k, "ra": ra, "t_min": float(t[int(np.argmin(worst))])}


def grid():
    for w in GRID_W:
        for ya in GRID_YA:
            for phn in GRID_PHI:
                for pc in GRID_PC:
                    for pw in GRID_PW:
                        yield SDesign(w, ya, phn, pc, pw)


def select_design(g: FGeo, side: int, zl: float, zh: float, eps_max: float, W_bound: float, closing=False):
    best = None
    for D in grid():
        if closing and D.ta > T_CLOSE_MAX:
            continue
        try:
            r = float_margin(g, side, D, zl, zh, eps_max, W_bound, closing=closing)
        except (FloatingPointError, OverflowError, ValueError):
            continue
        if r is None:
            continue
        if best is None or r["margin"] > best[1]["margin"]:
            best = (D, r)
    return best


def dyadic(x: float, bits=30) -> Fraction:
    return Fraction(round(x * 2 ** bits), 2 ** bits)


def contour_params(g: FGeo, side: int, D: SDesign, zl: float, zh: float):
    zref = math.sqrt(zl * zh)
    rT = g.turning(zref, side)
    k = abs(g.Ux(rT)) / (2 * zref ** 2)
    xa = g.rstar(rT) + side * float(D.ta) / k
    return dyadic(xa), dyadic(k)


# ---------------------------------------------------------------------------------------------
# certify a zeta-box (direct method)
# ---------------------------------------------------------------------------------------------
def eps_max_of(pt: Point, zl: Fraction) -> Fraction:
    """min(1/Omega_S, h0/zeta_lo) as an exact rational upper bound (Omega_S = 6 kappa_-, zeta = z Omega_P)."""
    a = 1 / (fa(OMEGA_S_KM) * pt.km)
    b = fa(H0) / (fa(zl) * pt.OmegaP) if zl > 0 else None
    v = a if b is None else a.min(b)
    return Fraction(int(float(v.upper()) * 2 ** 30) + 1, 2 ** 30)


def certify_eps(pt, side, D, C, zl_, zh_, el, eh, target, max_depth=MAX_DEPTH):
    prm = Params(pt, zl_, zh_, el, eh)
    ch = SChecker(pt, side, D, C, prm)
    mu, worst, n_eval, fails, depths = ch.cells(target, max_depth)
    tq, d = tail_bound(pt, side, D, C, prm, ch.k2)
    ok = bool(mu is not None and mu > 0 and not fails and tq is not None and tq > 0)
    return {"eps": [str(el), str(eh)], "mu_cells": lohi(mu) if mu is not None else None,
            "mu_tail": lohi(tq) if tq is not None else None, "tail_disc": d.str(4) if d is not None else None,
            "worst_cell": worst, "cells": n_eval, "depths": {str(k): v for k, v in sorted(depths.items())},
            "fails": fails[:10], "n_fails": len(fails), "pass": ok}


def certify_zbox(pt: Point, g: FGeo, side: int, zlo: Fraction, zhi: Fraction, W_bound: Fraction, design=None,
                 eps_levels=EPS_LEVELS, max_depth=MAX_DEPTH):
    zl_f, zh_f = float(zlo) * g.OmP, float(zhi) * g.OmP
    eps_max = eps_max_of(pt, zlo)
    rec = {"zeta_over_OmegaP": [str(zlo), str(zhi)], "eps_max": str(eps_max)}
    if design is None:
        best = select_design(g, side, zl_f, zh_f, float(eps_max), float(W_bound))
        if best is None:
            rec.update({"pass": False, "error": "no admissible design"})
            return rec
        D, fm = best
        Xa, k = contour_params(g, side, D, zl_f, zh_f)
        rec["float_margin"] = fm["margin"]
        if fm["margin"] <= 0:
            rec.update({"design": D.key(), "X_a": str(Xa), "k": str(k), "pass": False,
                        "error": "float design margin <= 0 (not attempted; the zeta-box is halved)"})
            return rec
    else:
        D, Xa, k = SDesign.from_key(design["design"]), Fraction(design["X_a"]), Fraction(design["k"])
    rec.update({"design": D.key(), "X_a": str(Xa), "k": str(k)})
    try:
        C = CurveContour(pt, side, D, Xa, k)
    except ArithmeticError as ex:
        rec.update({"pass": False, "error": f"contour: {ex}"})
        return rec
    zl_ = fa(zlo) * pt.OmegaP
    zh_ = fa(zhi) * pt.OmegaP
    rec["exit"] = exit_checks(pt, side, D, C, arb(zh_.upper()), eps_max, W_bound)
    rec["w1_check"] = w1_check(pt, side, D, C)
    tq = Fraction(design["target"]) if (design is not None and "target" in design) else \
        dyadic(max(rec.get("float_margin", 0.02) / 4, 1e-3), 20)
    rec["target"] = str(tq)
    target = fa(tq)
    # eps-boxes: halve on failure up to EPS_LEVELS
    leaves = []
    stack = [(Fraction(0), eps_max, 0)]
    while stack:
        el, eh, lev = stack.pop()
        r = certify_eps(pt, side, D, C, arb(zl_.lower()), arb(zh_.upper()), el, eh, target, max_depth)
        r["level"] = lev
        if not r["pass"] and lev < eps_levels:
            m = (el + eh) / 2
            stack.append((m, eh, lev + 1))
            stack.append((el, m, lev + 1))
            continue
        leaves.append(r)
    rec["eps_boxes"] = leaves
    ex = rec["exit"]
    rec["pass"] = bool(all(x["pass"] for x in leaves) and ex["exit_forbidden_ok"] and ex["exit_order_ok"]
                       and ex["phi_zero_at_exit"] and ex["eta_zero_at_exit"] and rec["w1_check"]["pass"])
    return rec


# ---------------------------------------------------------------------------------------------
# closing box (0, zeta_1]: near-horizon map
# ---------------------------------------------------------------------------------------------
NH_ORDER = 14


def nh_model(pt: Point, side: int, d_r: Fraction, N=NH_ORDER):
    """Near-horizon Taylor model.  With s = side (r - r_h) and rho = e^{+-2 kappa_h (x - c_h)} = s f(r),
    f = exp(side 2 kappa_h sum_{i != h} A_i log(1 + side s/(r_h - r_i))), G(rho) = U(r_h + side psi(rho))
    = sum_{n<=N} g_n rho^n + E_N, |E_N| <= M_G (|rho|/R)^{N+1}/(1 - |rho|/R) for |rho| <= R = m_f d_r
    (Rouche: psi(rho) stays in |r - r_h| <= d_r, where |f| >= m_f and |U| <= M_G).  Also M_W = sup|Wt|."""
    rh = pt.rp if side > 0 else pt.rc
    kap = pt.kp if side > 0 else pt.kc
    idx_h = 2 if side > 0 else 3
    prec = N + 2
    old_cap = ctx.cap
    ctx.cap = max(old_cap, prec + 2)
    try:
        L = acb_series([0], prec=prec)
        for i, (Ai, ri) in enumerate(zip(pt.Acoef, pt.roots)):
            if i == idx_h:
                continue
            L = L + Ai * acb_series([1, acb(side) / (rh - ri)], prec=prec).log()
        f_s = (side * 2 * kap * L).exp()
        rho_s = acb_series([0, 1], prec=prec) * f_s
        psi = rho_s.reversion()
        num = acb_series([acb(-pt.Lam / 3)], prec=prec)
        for i, ri in enumerate(pt.roots):
            if i != idx_h:
                num = num * acb_series([rh - ri, side], prec=prec)
        U_s = acb_series([0, side], prec=prec) * num / acb_series([rh, side], prec=prec) ** 4
        G = U_s(psi)
        g = G.coeffs()
    finally:
        ctx.cap = old_cap
    if len(g) < N + 1:
        raise ArithmeticError("near-horizon series truncated")
    g = [g[n] for n in range(N + 1)]
    # m_f on the box enclosing the disc: term-wise extremes of A_i log|r - r_i| (dist_bounds)
    d = fa(d_r)
    box = acb(rh + arb(0, 1) * d, arb(0, 1) * d)
    lo_sum = arb(0)
    for i, (Ai, ri) in enumerate(zip(pt.Acoef, pt.roots)):
        if i == idx_h:
            continue
        dlo, dhi = dist_bounds(box, ri)
        c = side * 2 * kap * Ai
        t1, t2 = c * (dlo / abs(rh - ri)).log(), c * (dhi / abs(rh - ri)).log()
        lo_sum = lo_sum + arb(t1.lower()).min(arb(t2.lower()))
    m_f = arb(arb(lo_sum.lower()).exp().lower())
    # sup |U|, sup |Wt| on the box by subdivision (n x n sub-boxes)
    nsub = 24
    M_G, M_W = arb(0), arb(0)
    for i in range(nsub):
        for j in range(nsub):
            cx = rh - d + (2 * i + 1) * d / nsub
            cy = -d + (2 * j + 1) * d / nsub
            sb = acb(cx + arb(0, 1) * (d / nsub), cy + arb(0, 1) * (d / nsub))
            U, W = pot_series(pt, sb)
            M_G = M_G.max(arb(abs(U).upper()))
            M_W = M_W.max(arb(abs(W).upper()))
    R = arb((m_f * d).lower())
    return {"g": g, "g1": g[1].real, "R": R, "M_G": M_G, "M_W": M_W, "m_f": m_f, "d_r": d_r, "N": N}


def nh_P(nh, rho):
    """P_N(rho) = sum_{n=1}^N g_n rho^(n-1) at a point/ball rho (acb)."""
    g = nh["g"]
    acc = acb(0)
    for n in range(nh["N"], 0, -1):
        acc = acc * rho + g[n]
    return acc


def nh_P_seg(nh, z):
    """Enclosure of {P_N(lam z) : lam in [0, 1]} as g_1 + sum_{n>=2} g_n z^(n-1) [0, 1] (each lam^(n-1) in
    [0, 1]); z an acb, or an acb_series [z, z'] in t, then [value, d/dt] (d/dt P(lam z(t)) =
    sum (n-1) g_n lam^(n-1) z^(n-2) z')."""
    g = nh["g"]
    U01 = arb(0).union(arb(1))
    if isinstance(z, acb_series):
        c = z.coeffs() + [acb(0)] * 2
        z0, z1 = c[0], c[1]
        val = acb(g[1])
        der = acb(0)
        zp = acb(1)
        for n in range(2, nh["N"] + 1):
            der += (n - 1) * g[n] * zp * z1 * U01
            zp = zp * z0
            val += g[n] * zp * U01
        return acb_series([val, der], prec=2)
    val = acb(g[1])
    zp = acb(1)
    for n in range(2, nh["N"] + 1):
        zp = zp * z
        val += g[n] * zp * U01
    return val


def nh_P_err(nh, rho_abs):
    """sup |E_N(rho)/rho| over |rho| <= rho_abs (< R)."""
    q = rho_abs / nh["R"]
    return arb(((nh["M_G"] / nh["R"]) * q ** nh["N"] / (1 - q)).upper())


def certify_closing(pt: Point, g: FGeo, side: int, W_bound: Fraction, design=None):
    """zeta in (0, zeta_1], eps in [0, 1/Omega_S]: contour x_T(zeta) + side Y(t)/kappa_h; along it
    rho = rho_T e^{2Y} and U~ = e^{2Y} P(rho_T e^{2Y})/P(rho_T), P = G/rho, rho_T in [0, rho_1]."""
    eps_max = eps_max_of(pt, Fraction(0))
    rec = {"zeta_over_OmegaP": ["0", str(ZETA_1)], "eps_max": str(eps_max), "closing": True}
    if design is None:
        best = select_design(g, side, 0.0, float(ZETA_1) * g.OmP, float(eps_max), float(W_bound), closing=True)
        D, fm = best
        rec["float_margin"] = fm["margin"]
    else:
        D = SDesign.from_key(design["design"])
    rec["design"] = D.key()
    kap = pt.kp if side > 0 else pt.kc
    nh = nh_model(pt, side, D_R[side])
    z1 = fa(ZETA_1) * pt.OmegaP
    g1 = nh["g1"]
    # rho_1: G(rho_1) >= zeta_1^2 with G increasing on [0, rho_1] (real); G >= rho (P_N - err)
    rho1 = arb(((z1 * z1) / g1).upper()) * arb("1.0001")
    rho1_ok = False
    for _ in range(40):
        Pl = nh_P(nh, acb(rho1)).real - nh_P_err(nh, rho1)
        if rho1 < nh["R"] and rho1 * Pl >= z1 * z1:
            rho1_ok = True
            break
        rho1 = rho1 * arb("1.05")
    ta = fa(D.ta)
    e2ta = (2 * ta).exp()
    rho_max = arb((rho1 * e2ta).upper())
    in_disc = bool(rho_max < nh["R"])
    rec["near_horizon"] = {"d_r": str(nh["d_r"]), "m_f": nh["m_f"].str(6), "R": nh["R"].str(6), "M_G": nh["M_G"].str(6),
                           "M_W": nh["M_W"].str(6), "order": nh["N"], "g1": g1.str(8), "g2": nh["g"][2].str(6),
                           "rho_1": rho1.str(6), "rho_max": rho_max.str(6), "rho_1_ok": rho1_ok, "in_disc": in_disc}
    if not (rho1_ok and in_disc):
        rec["pass"] = False
        rec["error"] = "near-horizon disc too small for this design"
        return rec
    errT = nh_P_err(nh, rho1)
    PT = nh_P_seg(nh, acb(rho1))                     # P_N(rho_T), rho_T = lam rho_1, lam in [0, 1]
    PT_abs_lo = arb((abs(PT).lower() - errT).lower())
    prm = Params(pt, None, None, Fraction(0), eps_max)
    k2 = kap * kap

    def model(tt, pb, pp):
        e1, _, _ = D.eta_derivs(tt, pb)
        eta = D.eta_value(tt, pb)
        e2Y = (2 * acb_series([acb(tt, eta), acb(1, e1)], prec=2)).exp()
        Ut = e2Y * nh_P_seg(nh, acb(rho1) * e2Y) / PT
        p = parts(pt, D, tt, None, None, side, kap, pb, pp, model=Ut)
        return p

    def err_cell(civ, p):
        """|Y'|^2 |U~_N| delta + eps^2 |Y'|^2 |Wt| over the cell."""
        e2t = (2 * civ).exp()
        rabs = arb((rho1 * e2t).upper())
        PN = nh_P_seg(nh, acb(rho1 * e2t))
        e1 = nh_P_err(nh, rabs) / arb((abs(PN).lower()).max(arb("1e-300")))
        e2 = errT / PT_abs_lo
        delta = (e1 + e2) / (1 - e2)
        Uabs = e2t * arb((abs(PN).upper() + nh_P_err(nh, rabs)).upper()) / PT_abs_lo
        q = rabs / nh["R"]
        Wabs = nh["M_W"] * q / (1 - q)
        return p["Yp2"][0] * (Uabs * delta + prm.E2.upper() * Wabs)

    def cell(a: Fraction, b: Fraction):
        pb = D.piece_bend(a, b)
        pp = D.piece_phase(a, b)
        a_, b_ = fa(a), fa(b)
        w = b_ - a_
        tm = (a_ + b_) / 2
        civ = tm + arb(0, 1) * (w / 2)
        pm = model(tm, pb, pp)
        pc = model(civ, pb, pp)
        qm = prm.combine(pm, 0, k2)
        dq = prm.combine(pc, 1, k2)
        err = err_cell(civ, pc)
        return arb((qm.lower() - abs(dq).upper() * (w / 2) - err.upper()).lower())

    # exit: U~(t_a) - 1 - eps_max^2 W_bound > 0 for all rho_T in [0, rho_1]
    pa = model(ta, 1, 1)
    Ua_lo = pa["P1"][0]                              # = Re(G U~) with G = 1 at the exit (phi = 0, Y' = 1)
    ea = err_cell(ta, pa)
    ex_margin = Ua_lo - ea - 1 - fa(eps_max) ** 2 * fa(W_bound)
    rh = pt.rp if side > 0 else pt.rc
    r_edge = rh + side * fa(nh["d_r"])
    x_edge = rstar(pt, acb(r_edge))
    x0re = sum((Ai * abs(ri).log() for Ai, ri in zip(pt.Acoef, pt.roots)), arb(0))
    strip = bool(x_edge < x0re) if side > 0 else bool(x_edge > 0)
    order = bool(r_edge < pt.rP) if side > 0 else bool(r_edge > pt.rP)
    rec["exit"] = {"Ut_minus_1_minus_epsmax2_Wbound": lohi(ex_margin), "exit_forbidden_ok": bool(ex_margin > 0),
                   "exit_order_ok": order, "phi_zero_at_exit": bool(D.p1 <= D.ta), "eta_zero_at_exit": bool(D.w <= D.ta)}
    rec["w1_check"] = {"note": "closing box: r(x(t)) = r_h + side psi(rho_T e^{2Y}) stays in |r - r_h| <= d_r "
                               "(Rouche); x there lies beyond the disc-edge value",
                       "x_at_disc_edge": x_edge.str(8), "Re_rstar_0": x0re.str(8),
                       "strip": "Re x < Re r*(0)" if side > 0 else "Re x > 0", "pass": strip}
    nodes = sorted({T_LEFT, D.ta} | {x for x in D.breakpoints() if T_LEFT < x < D.ta}
                   | {m * H_T for m in range(math.ceil(T_LEFT / H_T), math.floor(D.ta / H_T) + 1) if T_LEFT < m * H_T < D.ta})
    tq = Fraction(design["target"]) if (design is not None and "target" in design) else \
        dyadic(max(rec.get("float_margin", 0.1) / 4, 1e-3), 20)
    rec["target"] = str(tq)
    target = fa(tq)
    mu, worst, n_eval, fails = None, None, 0, []
    for a, b in itertools.pairwise(nodes):
        stack = [(a, b, 0)]
        while stack:
            a_, b_, dpt = stack.pop()
            low = cell(a_, b_)
            n_eval += 1
            if low < target and dpt < MAX_DEPTH:
                m = (a_ + b_) / 2
                stack.append((m, b_, dpt + 1))
                stack.append((a_, m, dpt + 1))
                continue
            if not (low > 0):
                fails.append([str(a_), str(b_), low.str(5)])
            if mu is None or low < mu:
                mu, worst = low, (float(a_), float(b_))
    # tail t <= t_L: |rho| <= rho_1 e^{2 t_L}; U~ = e^{2Y} P(rho)/P(rho_T), |e^{2Y}| = e^{2t} <= e^{2 t_L}
    phinf = D.phinf()
    e2L = (2 * fa(T_LEFT)).exp()
    rL = arb((rho1 * e2L).upper())
    UL = e2L * arb((abs(nh_P_seg(nh, acb(rL))).upper() + nh_P_err(nh, rL)).upper()) / PT_abs_lo
    WL = nh["M_W"] * (rL / nh["R"]) / (1 - rL / nh["R"])
    tail = (-2 * phinf.sin() + prm.EN * (-4 * phinf.cos()) + prm.E2N2 * (2 * phinf.sin())
            - 2 * UL - 2 * WL * prm.E2.upper())
    rec.update({"mu_cells": lohi(mu), "mu_tail": lohi(tail), "worst_cell": worst, "cells": n_eval,
                "fails": fails[:10], "n_fails": len(fails)})
    ex = rec["exit"]
    rec["pass"] = bool(mu > 0 and tail > 0 and not fails and ex["exit_forbidden_ok"] and ex["exit_order_ok"]
                       and ex["phi_zero_at_exit"] and ex["eta_zero_at_exit"] and strip)
    return rec


# ---------------------------------------------------------------------------------------------
# runs
# ---------------------------------------------------------------------------------------------
def zboxes():
    out = []
    z = ZETA_1
    while z < ZETA_SU:
        out.append((z, min(z + ZBOX, ZETA_SU)))
        z += ZBOX
    return out


def certify_zbox_refined(pt, g, side, zlo, zhi, W_bound):
    rec = certify_zbox(pt, g, side, zlo, zhi, W_bound)
    if rec["pass"] or (zhi - zlo) / 2 < ZBOX_MIN:
        return [rec]
    m = (zlo + zhi) / 2
    return certify_zbox_refined(pt, g, side, zlo, m, W_bound) + certify_zbox_refined(pt, g, side, m, zhi, W_bound)


def run_side(side_name: str, worker: int = 0, nworkers: int = 1):
    from .part1 import load_bulk
    side = 1 if side_name == "+" else -1
    W_bound, _ = load_bulk()
    tag = "plus" if side > 0 else "c"
    RES_DIR.mkdir(parents=True, exist_ok=True)
    out = RES_DIR / f"S_PB_{tag}.w{worker}.jsonl"
    done = set()
    if out.exists():
        for line in out.read_text().splitlines():
            if line.strip():
                done.add(tuple(json.loads(line)["key"]))
    pt = Point("PB", 256)
    g = FGeo(pt)
    jobs = [("closing", "0", str(ZETA_1))] + [("box", str(a), str(b)) for a, b in zboxes()]
    with precision(PREC):
        for j, key in enumerate(jobs):
            if j % nworkers != worker or key in done:
                continue
            with Timer() as tm:
                if key[0] == "closing":
                    leaves = [certify_closing(pt, g, side, W_bound)]
                else:
                    leaves = certify_zbox_refined(pt, g, side, Fraction(key[1]), Fraction(key[2]), W_bound)
            row = {"key": list(key), "pass": all(x["pass"] for x in leaves), "leaves": leaves, "commit": git_commit(),
                   "timing_s": {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)},
                   "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}
            with out.open("a") as f:
                f.write(json.dumps(row) + "\n")
            print(side_name, key, row["pass"], row["timing_s"], flush=True)


def assemble(side_name: str):
    side = 1 if side_name == "+" else -1
    tag = "plus" if side > 0 else "c"
    rows = {}
    for f in sorted(RES_DIR.glob(f"S_PB_{tag}.w*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[tuple(r["key"])] = r
    jobs = [("closing", "0", str(ZETA_1))] + [("box", str(a), str(b)) for a, b in zboxes()]
    missing = [j for j in jobs if j not in rows]
    leaves = [x for j in jobs if j in rows for x in rows[j]["leaves"]]
    # the leaves tile [0, zeta_SU]
    iv = sorted((Fraction(x["zeta_over_OmegaP"][0]), Fraction(x["zeta_over_OmegaP"][1])) for x in leaves)
    tiles = bool(iv) and iv[0][0] == 0 and iv[-1][1] == ZETA_SU and all(a[1] == b[0] for a, b in itertools.pairwise(iv))
    ok = (not missing) and tiles and all(x["pass"] for x in leaves)
    cert = {"regime": "S", "side": tag, "point": "PB", "h0": str(H0), "precision_bits": PREC,
            "design_grid": {"w": [str(x) for x in GRID_W], "y_a": [str(x) for x in GRID_YA],
                            "phi_inf": [f"-pi/{n}" for n in GRID_PHI], "phase_centre": [str(x) for x in GRID_PC],
                            "phase_width": [str(x) for x in GRID_PW]},
            "constants": {"Omega_S_over_km": str(OMEGA_S_KM), "zeta_SU_over_OmegaP": str(ZETA_SU),
                          "zeta_1_over_OmegaP": str(ZETA_1), "zeta_box": str(ZBOX), "zeta_box_min": str(ZBOX_MIN),
                          "eps_levels": EPS_LEVELS, "t_L": str(T_LEFT), "h_t": str(H_T), "taylor_order": N_TAYLOR,
                          "cell_depth": MAX_DEPTH, "closing_t_a_max": str(T_CLOSE_MAX)},
            "missing": missing, "zeta_leaves_tile_0_to_zeta_SU": tiles, "n_zeta_leaves": len(leaves),
            "w1_check": [{"zeta": x["zeta_over_OmegaP"], **x.get("w1_check", {})} for x in leaves],
            "cpu_s_total": round(sum(rows[j]["timing_s"]["cpu"] for j in jobs if j in rows), 1),
            "commits": sorted({rows[j]["commit"] for j in jobs if j in rows}), "pass": ok, "leaves": leaves}
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    (CERT_DIR / f"S_PB_{tag}.json").write_text(json.dumps(cert, indent=1) + "\n")
    return cert
