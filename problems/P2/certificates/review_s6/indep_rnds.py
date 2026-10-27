"""S6 review: an independent toolkit for the RNdS radial equation at P_B.

Shares no code with p2: python-flint (Arb) and the standard library only.  Everything is derived here
from FORMULATION v1.1 section 1:

    Lambda = 3/50, M = 1;  r_ext = smallest positive root of r^3 - 25 r + 25;  Q_max^2 = r_ext - r_ext^4/50;
    Q^2 = (999/1000)^2 Q_max^2;  Delta(r) = -(Lambda/3) r^4 + r^2 - 2 r + Q^2;  kappa_h = |Delta'(r_h)|/(2 r_h^2).

The phase-factored equations (derived independently in sympy_lemmaE.py):
    side +1 (Z_+ = e^{i w (x - c_+)} R_+):  Delta Z'' + (Delta' - 2 i w r^2) Z' - (2 i w r + L) Z = 0,
    side -1 (Z_c = e^{-i w (x - c_c)} R_c): Delta Z'' + (Delta' + 2 i w r^2) Z' + (2 i w r - L) Z = 0,
and W_l(w) = [Delta (Z_+ Z_c' - Z_+' Z_c) + 2 i w r^2 Z_+ Z_c]/(r_+ r_c) at any r in (r_+, r_c).

Validated series: Frobenius at r_+ and r_c, Taylor continuation at regular points, each with a
geometric majorant tail (the recurrence-coefficient bounds g_i hold for every index above the
truncation; q solves sum g_i q^-i <= 1, checked in Arb).
"""
from __future__ import annotations

import math
from fractions import Fraction

from flint import acb, arb, ctx, fmpq

I = acb(0, 1)


def A(x) -> arb:
    if isinstance(x, Fraction):
        return arb(fmpq(x.numerator, x.denominator))
    return arb(x)


def up(x: arb) -> arb:
    """an exact arb upper bound of x (x real ball)"""
    return arb(x.upper())


def newton_root(f, fp, x0: float, w: float = 1e-6, iters: int = 12) -> arb:
    """Interval Newton: returns a ball containing the unique root of f in [x0 - w, x0 + w]."""
    X = arb(x0, w)
    for k in range(iters):
        m = arb(X.mid())
        d = fp(X)
        if d.contains(0):
            raise ArithmeticError("f' contains 0 on the Newton box")
        N = m - f(m) / d
        if not (N.lower() > X.lower() and N.upper() < X.upper()):
            if k == 0:
                raise ArithmeticError(f"Newton box not contracted at x0 = {x0}")
            break
        X = N
    return X


class Point:
    """P_B from its definitions (FORMULATION section 1.3), as Arb balls at the current precision."""

    def __init__(self, qfrac=Fraction(999, 1000)):
        self.lam3 = A(Fraction(1, 50))                       # Lambda/3, Lambda = 3/50
        self.rext = newton_root(lambda r: r ** 3 - 25 * r + 25, lambda r: 3 * r ** 2 - 25, 1.0457442422066)
        self.Qmax2 = self.rext - self.rext ** 4 / 50
        self.Q2 = A(qfrac) ** 2 * self.Qmax2
        self.Dc = [self.Q2, arb(-2), arb(1), arb(0), -self.lam3]     # Delta, ascending powers of r
        D = lambda r: (((-self.lam3 * r) * r + 1) * r - 2) * r + self.Q2
        Dp = lambda r: ((-4 * self.lam3 * r) * r + 2) * r - 2
        self.Delta, self.Deltap = D, Dp
        guesses = {"n": -7.96033597, "m": 0.99736659, "p": 1.09434835, "c": 5.86862103}
        self.r = {k: newton_root(D, Dp, g, 1e-5) for k, g in guesses.items()}
        self.rn, self.rm, self.rp, self.rc = (self.r[k] for k in "nmpc")
        self.kappa = {k: abs(Dp(self.r[k])) / (2 * self.r[k] ** 2) for k in "mpc"}
        self.km, self.kp, self.kc = self.kappa["m"], self.kappa["p"], self.kappa["c"]
        # photon sphere r_P = (3 + sqrt(9 - 8 Q^2))/2 (M = 1), Omega_P = sqrt(F(r_P))/r_P, gamma0 = lambda_P/2
        self.rP = (3 + (9 - 8 * self.Q2).sqrt()) / 2
        F = lambda r: D(r) / r ** 2
        self.F = F
        self.OmegaP = F(self.rP).sqrt() / self.rP
        self.gamma0 = (F(self.rP) * (2 - 3 / self.rP)).sqrt() / self.rP / 2
        # sign checks of the root ordering
        assert bool(self.rn < 0) and bool(0 < self.rm) and bool(self.rm < self.rp) and bool(self.rp < self.rc)

    def summary(self) -> dict:
        s = lambda x: x.str(25, radius=True)
        return {"Q2": s(self.Q2), "r_n": s(self.rn), "r_-": s(self.rm), "r_+": s(self.rp), "r_c": s(self.rc),
                "kappa_-": s(self.km), "kappa_+": s(self.kp), "kappa_c": s(self.kc), "r_P": s(self.rP),
                "Omega_P": s(self.OmegaP), "gamma0": s(self.gamma0),
                "kappa_+/kappa_-": s(self.kp / self.km), "kappa_c/kappa_-": s(self.kc / self.km),
                "Omega_P/kappa_-": s(self.OmegaP / self.km), "gamma0/kappa_-": s(self.gamma0 / self.km)}


# ------------------------------------------------------------------------------------------------
# polynomials (ascending coefficient lists of acb)

def padd(p, q):
    n = max(len(p), len(q))
    return [(p[k] if k < len(p) else acb(0)) + (q[k] if k < len(q) else acb(0)) for k in range(n)]


def pscale(p, c):
    return [c * x for x in p]


def shift(p, x0):
    """coefficients of p(x0 + h) in h"""
    n = len(p)
    out = []
    for j in range(n):
        s = acb(0)
        for k in range(j, n):
            s += p[k] * math.comb(k, j) * (x0 ** (k - j))
        out.append(s)
    return out


def eq_polys(pt: Point, w: acb, L: int, side: int):
    """(Delta, B, C) as polynomials in r for the side-s phase-factored equation."""
    Dl = [acb(c) for c in pt.Dc]
    Dp = [Dl[k] * k for k in range(1, len(Dl))]
    r2 = [acb(0), acb(0), acb(1)]
    r1 = [acb(0), acb(1)]
    B = padd(Dp, pscale(r2, -2 * side * I * w))
    C = padd(pscale(r1, -2 * side * I * w), [acb(-L)])
    return Dl, B, C


def _solve_q(g: list) -> arb:
    """q > 0 with sum_i g_i q^-i <= 1 (i = 1..len(g)), verified in Arb."""
    gf = [float(x.upper()) for x in g]
    f = lambda q: sum(gi * q ** (-(i + 1)) for i, gi in enumerate(gf))
    lo, hi = 1e-12, 1.0
    while f(hi) > 1:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 1:
            lo = mid
        else:
            hi = mid
    q = arb(hi * (1 + 1e-9))
    tot = sum((g[i] * q ** (-(i + 1)) for i in range(len(g))), arb(0))
    if not bool(tot <= 1):
        raise ArithmeticError("majorant root not verified")
    return q


def _tails(K: arb, q: arb, h: arb, N: int):
    """sum_{n>N} K q^n h^n and sum_{n>N} n K q^n h^(n-1), with x = q h < 1."""
    x = q * h
    if not bool(x < 1):
        raise ArithmeticError("step outside the majorant disc")
    t0 = K * x ** (N + 1) / (1 - x)
    t1 = K * q * x ** N * ((N + 1) - N * x) / (1 - x) ** 2
    return up(t0), up(t1)


def _ball(r: arb) -> acb:
    return acb(arb(0, r), arb(0, r))


def frobenius(pt: Point, rh: arb, w: acb, L: int, side: int, a_init, t_eval: arb, N: int, b0_zero: bool = False):
    """Frobenius series Z = sum a_m t^m at the regular singular point rh (Delta(rh) = 0), t = r - rh,
    started from a_init = [a_0] (exponent 0, D_m = m (Delta_1 (m - 1) + B_0) != 0 for m >= 1) or
    a_init = [0, 1] with b0_zero (the exponent-1 solution at the E_+ point, where B_0 = 0 exactly).
    Returns (Z, dZ/dr) at t_eval (acb balls, rigorous tails) and the number of terms."""
    Dl, B, C = eq_polys(pt, w, L, side)
    Ds, Bs, Cs = shift(Dl, rh), shift(B, rh), shift(C, rh)
    Ds[0] = acb(0)                                          # Delta(r_h) = 0 exactly
    D1 = Ds[1]
    B0 = acb(0) if b0_zero else Bs[0]
    a = list(a_init)
    m0 = len(a)
    for m in range(m0, N + 1):
        acc = acb(0)
        for k in range(2, 5):                                # Delta_k (m - k + 1)(m - k) a_{m-k+1}
            j = m - k + 1
            if 0 <= j < len(a) and k < len(Ds):
                acc += Ds[k] * ((m - k + 1) * (m - k)) * a[j]
        for k in range(1, 4):                                # B_k (m - k) a_{m-k}
            j = m - k
            if 0 <= j < len(a) and k < len(Bs):
                acc += Bs[k] * (m - k) * a[j]
        for k in range(0, 2):                                # C_k a_{m-1-k}
            j = m - 1 - k
            if 0 <= j < len(a) and k < len(Cs):
                acc += Cs[k] * a[j]
        Dm = m * (D1 * (m - 1) + B0)
        a.append(-acc / Dm)
    # majorant: for m >= N + 1, |a_m| <= sum_{i=1..3} g_i |a_{m-i}|
    beta = abs(B0) / abs(D1)
    if not bool(N - beta > 1):
        raise ArithmeticError("N too small")
    g = []
    for i in range(1, 4):
        Di = abs(Ds[i + 1]) if i + 1 < len(Ds) else arb(0)
        Bi = abs(Bs[i]) if i < len(Bs) else arb(0)
        Ci = abs(Cs[i - 1]) if i - 1 < len(Cs) else arb(0)
        g.append(up((Di * (N + 1) + Bi + Ci / (N + 1)) / (abs(D1) * (N - beta))))
    q = _solve_q(g)
    K = max((up(abs(a[n]) / q ** n) for n in range(N - 2, N + 1)), key=lambda x: float(x.upper()))
    h = abs(t_eval)
    t0, t1 = _tails(K, q, h, N)
    Z = acb(0)
    dZ = acb(0)
    tt = acb(t_eval)
    for n in range(N, -1, -1):                               # Horner
        Z = Z * tt + a[n]
    for n in range(N, 0, -1):
        dZ = dZ * tt + a[n] * n
    return Z + _ball(t0), dZ + _ball(t1), {"terms": N, "q_times_h": float((q * h).upper())}


def taylor_step(pt: Point, x0: arb, Z0: acb, Z1: acb, w: acb, L: int, side: int, h: arb, N: int):
    """(Z, Z') at x0 + h from (Z, Z') at the regular point x0 (Delta(x0) != 0)."""
    Dl, B, C = eq_polys(pt, w, L, side)
    d, b, c = shift(Dl, acb(x0)), shift(B, acb(x0)), shift(C, acb(x0))
    d0 = d[0]
    if abs(d0).contains(0):
        raise ArithmeticError("Delta(x0) not separated from 0")
    z = [Z0, Z1]
    for m in range(2, N + 1):
        acc = acb(0)
        for i in range(1, 5):
            j = m - i
            if j < 0:
                continue
            t = acb(0)
            if i < len(d):
                t += d[i] * ((m - i) * (m - i - 1))
            if 0 <= i - 1 < len(b):
                t += b[i - 1] * (m - i)
            if 0 <= i - 2 < len(c):
                t += c[i - 2]
            acc += t * z[j]
        z.append(-acc / (d0 * (m * (m - 1))))
    g = []
    for i in range(1, 5):
        di = abs(d[i]) if i < len(d) else arb(0)
        bi = abs(b[i - 1]) if 0 <= i - 1 < len(b) else arb(0)
        ci = abs(c[i - 2]) if 0 <= i - 2 < len(c) else arb(0)
        g.append(up((di + bi / N + ci / (N * (N + 1))) / abs(d0)))
    q = _solve_q(g)
    K = max((up(abs(z[n]) / q ** n) for n in range(N - 3, N + 1)), key=lambda x: float(x.upper()))
    t0, t1 = _tails(K, q, abs(h), N)
    Z = acb(0)
    dZ = acb(0)
    hh = acb(h)
    for n in range(N, -1, -1):
        Z = Z * hh + z[n]
    for n in range(N, 0, -1):
        dZ = dZ * hh + z[n] * n
    return Z + _ball(t0), dZ + _ball(t1), float((q * abs(h)).upper())


def step_q(pt: Point, x0: arb, w: acb, L: int, side: int, N: int) -> float:
    """the majorant growth rate q at x0 (to choose a step with q h <= 1/2)"""
    Dl, B, C = eq_polys(pt, w, L, side)
    d, b, c = shift(Dl, acb(x0)), shift(B, acb(x0)), shift(C, acb(x0))
    g = []
    for i in range(1, 5):
        di = abs(d[i]) if i < len(d) else arb(0)
        bi = abs(b[i - 1]) if 0 <= i - 1 < len(b) else arb(0)
        ci = abs(c[i - 2]) if 0 <= i - 2 < len(c) else arb(0)
        g.append(up((di + bi / N + ci / (N * (N + 1))) / abs(d[0])))
    return float(_solve_q(g).upper())


def continue_to(pt: Point, x0: Fraction, Z0: acb, Z1: acb, x1: Fraction, w: acb, L: int, side: int,
                N: int = 120, ratio: float = 0.45):
    """validated continuation from x0 to x1 along the real axis (both in (r_+, r_c)), exact dyadic nodes"""
    x = x0
    steps = 0
    worst = 0.0
    sgn = 1 if x1 > x0 else -1
    while x != x1:
        q = step_q(pt, A(x), w, L, side, N)
        hmax = ratio / q
        hq = Fraction(math.floor(hmax * 2 ** 20), 2 ** 20)
        if hq <= 0:
            raise ArithmeticError("step underflow")
        nxt = x + sgn * hq
        if (nxt - x1) * sgn > 0:
            nxt = x1
        Z0, Z1, qh = taylor_step(pt, A(x), Z0, Z1, w, L, side, A(nxt - x), N)
        worst = max(worst, qh)
        x = nxt
        steps += 1
    return Z0, Z1, {"steps": steps, "max_q_h": worst}


def rat(x: arb, den: int = 2 ** 24) -> Fraction:
    return Fraction(round(float(x.mid()) * den), den)


def wronskian(pt: Point, w: acb, L: int, Zp, dZp, Zc, dZc, rm: arb) -> acb:
    D = acb(pt.Delta(rm))
    return (D * (Zp * dZc - dZp * Zc) + 2 * I * w * acb(rm * rm) * Zp * Zc) / acb(pt.rp * pt.rc)


def solve_W(pt: Point, w: acb, l: int, which: str = "Z+", N_frob: int = 160, N_tay: int = 120,
            t_frac: float = 0.2, s_frac: float = 0.15, rm_frac: Fraction = Fraction(1, 2)):
    """W_l(w) (which = 'Z+': the exponent-0 Jost solution a_0 = 1) or Wt_l (which = 'Z2': the exponent-1
    solution at the E_+ point, w must contain omega0 = -i kappa_+ exactly)."""
    L = l * (l + 1)
    rp, rc = pt.rp, pt.rc
    t0 = Fraction(round(t_frac * float((pt.rp - pt.rm).mid()) * 2 ** 20), 2 ** 20)      # inside r_+ - r_-
    s0 = Fraction(round(s_frac * float((pt.rc - pt.rp).mid()) * 2 ** 20), 2 ** 20)      # inside r_c - r_+
    xm = rat(rp + (rc - rp) * A(rm_frac), 2 ** 20)
    xp = rat(rp, 2 ** 40) + t0             # an exact dyadic near r_+ + t0
    Zp, dZp, i1 = (frobenius(pt, rp, w, L, +1, [acb(0), acb(1)] if which == "Z2" else [acb(1)],
                             A(xp) - rp, N_frob, b0_zero=(which == "Z2")))
    Zp, dZp, c1 = continue_to(pt, xp, Zp, dZp, xm, w, L, +1, N_tay)
    xc = rat(rc, 2 ** 40) - s0
    Zc, dZc, i2 = frobenius(pt, rc, w, L, -1, [acb(1)], A(xc) - rc, N_frob)
    Zc, dZc, c2 = continue_to(pt, xc, Zc, dZc, xm, w, L, -1, N_tay)
    W = wronskian(pt, w, L, Zp, dZp, Zc, dZc, A(xm))
    info = {"frob_plus": i1, "frob_c": i2, "cont_plus": c1, "cont_c": c2, "r_match": str(xm),
            "x_start_plus": str(xp), "x_start_c": str(xc)}
    return W, info
