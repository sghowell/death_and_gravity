"""S3 review: an independent floating-point W_l(omega) (reviewer's code; imports nothing from p2).

Formulation (different from both S1 arbw and S3 zsol): the unfactored radial equation as the
first-order system in r
    u = R,  v = Delta R',   u' = v/Delta,   v' = -(w^2 r^4/Delta - L) u,
so that Delta (R_+ R_c' - R_+' R_c) = R_+ v_c - v_+ R_c.  Start values from Frobenius series
R = (+-(r - r_h))^rho sum a_n (r - r_h)^n, rho = -i w/(2 kappa_h), computed from the equation
Delta^2 R'' + Delta Delta' R' + (w^2 r^4 - L Delta) R = 0 (vectorised over w, complex128), then
scipy DOP853 (rtol 1e-12) to the matching point.  W_l includes the frozen Jost factor
exp(-i w c_+ + i w c_c)/(r_+ r_c) of FORMULATION section 1.2 (c_+, c_c from A_i = r_i^2/Delta'(r_i)).
Parameters from the exact definitions of FORMULATION section 1.3 (mpmath, 50 digits).
Not rigorous; used only as an independent cross-check of the certified counts and increments.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

mp.mp.dps = 50


class Geo:
    def __init__(self, qfrac=Fraction(999, 1000)):
        cands = [mp.findroot(lambda r: r ** 3 - 25 * r + 25, x0) for x0 in (1.0, 1.2, 4.0)]
        rext = min(c for c in cands if c > 0)
        qmax2 = rext - rext ** 4 / 50
        Q2 = mp.mpf(qfrac.numerator) ** 2 / mp.mpf(qfrac.denominator) ** 2 * qmax2
        self.dco_mp = [Q2, mp.mpf(-2), mp.mpf(1), mp.mpf(0), -mp.mpf(1) / 50]   # ascending
        roots = sorted(mp.re(x) for x in mp.polyroots(self.dco_mp[::-1], maxsteps=200, extraprec=200))
        self.roots_mp = roots
        self.rn, self.rm_, self.rp, self.rc = [float(x) for x in roots]
        dD = lambda r: sum(k * self.dco_mp[k] * r ** (k - 1) for k in range(1, 5))
        self.A_mp = [r * r / dD(r) for r in roots]                 # r* = sum A_i log|r - r_i|
        rp, rc = roots[2], roots[3]
        self.km = float(abs(dD(roots[1])) / (2 * roots[1] ** 2))
        self.kp = float(abs(dD(rp)) / (2 * rp ** 2))
        self.kc = float(abs(dD(rc)) / (2 * rc ** 2))
        self.c_plus = float(sum(self.A_mp[i] * mp.log(abs(rp - roots[i])) for i in (0, 1, 3)))
        self.c_c = float(sum(self.A_mp[i] * mp.log(abs(rc - roots[i])) for i in (0, 1, 2)))
        self.dco = np.array([float(x) for x in self.dco_mp])
        rP = (3 + mp.sqrt(9 - 8 * Q2)) / 2
        self.rP = float(rP)

    def shifted(self, c: float) -> np.ndarray:
        """Coefficients d_k of Delta(c + x) = sum d_k x^k (k = 0..4), computed in mpmath."""
        cm = mp.mpf(c)
        out = []
        for k in range(5):
            out.append(float(sum(comb(i, k) * self.dco_mp[i] * cm ** (i - k) for i in range(k, 5))))
        return np.array(out)

    def Delta(self, r):
        return sum(self.dco[k] * r ** k for k in range(5))


def _polymul(a, b):
    out = np.zeros(len(a) + len(b) - 1, dtype=complex)
    for i, x in enumerate(a):
        out[i:i + len(b)] += x * b
    return out


def frobenius_start(g: Geo, which: str, w: np.ndarray, L: int, x: float, nmax=600, tol=1e-19):
    """(u, v) = (R, Delta R') at r = r_h + x for the frozen horizon solution, vectorised over w."""
    rh_mp = g.roots_mp[2] if which == "+" else g.roots_mp[3]
    rh = float(rh_mp)
    kap = g.kp if which == "+" else g.kc
    d = g.shifted(rh)          # d[0] = Delta(r_h) ~ 0: use exactly 0
    p = d[1:]                  # P(x) = Delta/x = p0 + p1 x + p2 x^2 + p3 x^3
    U = _polymul(p, p)                                   # P^2, degree 6
    Q = np.array([(k + 1) * p[k] for k in range(4)])
    V = _polymul(p, Q)                                   # P (P + x P'), degree 6
    r4 = np.array([comb(4, k) * rh ** (4 - k) for k in range(5)])
    xP = np.concatenate([[0.0], p])                      # x P, degree 4
    d1 = p[0]
    rho = -1j * w / (2 * kap)
    W2 = w * w
    nw = len(w)
    a = [np.ones(nw, dtype=complex)]
    S = np.ones(nw, dtype=complex)
    Sd = np.zeros(nw, dtype=complex)
    xp = 1.0
    small = 0
    for m in range(1, nmax):
        acc = np.zeros(nw, dtype=complex)
        for n in range(max(0, m - 6), m):
            j = m - n
            term = np.zeros(nw, dtype=complex)
            if j < len(U):
                term = term + U[j] * (n + rho) * (n + rho - 1)
            if j < len(V):
                term = term + V[j] * (n + rho)
            if j < 5:
                term = term + W2 * r4[j] - L * xP[j]
            acc = acc + a[n] * term
        am = -acc / (d1 * d1 * m * (m + 2 * rho))
        a.append(am)
        Sd = Sd + m * am * xp          # x^(m-1)
        xp = xp * x
        tv = am * xp
        S = S + tv
        if m > 20 and np.all(np.abs(tv) < tol * np.abs(S)) and np.all(np.abs(m * am * xp / x) < tol * (np.abs(Sd) + np.abs(S))):
            small += 1
            if small >= 6:
                break
        else:
            small = 0
    else:
        raise RuntimeError("Frobenius series did not converge")
    sx = x if which == "+" else -x                       # r - r_+ > 0 ; r_c - r > 0
    pref = np.exp(rho * np.log(sx))
    R = pref * S
    dR = pref * (rho * S / x + Sd)                       # d/dr of (+-x)^rho g
    Dl = sum(d[k] * x ** k for k in range(1, 5))
    return R, Dl * dR


def propagate(g: Geo, which: str, w: np.ndarray, L: int, rmatch: float, x0=None, rtol=1e-12):
    rh = g.rp if which == "+" else g.rc
    if x0 is None:
        x0 = 0.03 if which == "+" else -0.6
    u0, v0 = frobenius_start(g, which, w, L, x0)
    W2 = w * w
    n = len(w)

    def rhs(r, y):
        D = g.Delta(r)
        u, v = y[:n], y[n:]
        return np.concatenate([v / D, -(W2 * r ** 4 / D - L) * u])

    sol = solve_ivp(rhs, (rh + x0, rmatch), np.concatenate([u0, v0]), method="DOP853",
                    rtol=rtol, atol=1e-30)
    if not sol.success:
        raise RuntimeError(sol.message)
    y = sol.y[:, -1]
    return y[:n], y[n:]


def W_l(g: Geo, l: int, w, rmatch=None, rtol=1e-12):
    """Frozen W_l at the absolute frequencies w (array), FORMULATION section 1.2 normalisation."""
    w = np.atleast_1d(np.asarray(w, dtype=complex))
    L = l * (l + 1)
    rm = g.rP if rmatch is None else rmatch
    up, vp = propagate(g, "+", w, L, rm, rtol=rtol)
    uc, vc = propagate(g, "c", w, L, rm, rtol=rtol)
    W = (up * vc - vp * uc) * np.exp(-1j * w * g.c_plus + 1j * w * g.c_c) / (g.rp * g.rc)
    return W


def arg_sum(Wv: np.ndarray) -> float:
    """Sum of principal arguments of consecutive ratios (the continuous argument change)."""
    return float(np.sum(np.angle(Wv[1:] / Wv[:-1])))
