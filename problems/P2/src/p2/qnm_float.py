"""Float (mpmath) QNM solver for the massless scalar on RNdS -- evidence, never a certificate.

Method (as in src/s0_sanity.py, with adaptive radii so that large l, large |omega| and near-
extremal points do not lose digits to cancellation):
  * Frobenius series at r_+ and r_c with exponents -i w/(2 kappa) (FORMULATION section 1.2),
    evaluated at a distance t0 <= 0.25 d_min and t0 <= K d_min / |rho| from the horizon;
  * Taylor-series continuation of (Delta R')' + (w^2 r^4/Delta - L) R = 0 with step
    h <= frac * (distance to the nearest singular point) and h <= K / k_loc,
    k_loc = |w| r^2/|Delta| + sqrt(L)/r  (local oscillation/growth rate in r);
  * the Jost-normalised Wronskian W_l (FORMULATION section 1.2) at a matching point.
Zeros are found with the secant method from seeds.
"""
from __future__ import annotations

import mpmath as mp


class Geom:
    """RNdS with M = 1: Delta(r) = -(Lam/3) r^4 + r^2 - 2 r + Q^2 (mpmath, current dps)."""

    def __init__(self, Lam, Q2):
        self.Lam, self.Q2 = mp.mpf(Lam), mp.mpf(Q2)
        self.dco = [self.Q2, mp.mpf(-2), mp.mpf(1), mp.mpf(0), -self.Lam / 3]   # ascending
        rts = mp.polyroots(self.dco[::-1], maxsteps=800, extraprec=4 * mp.mp.prec)
        if any(abs(mp.im(z)) > mp.mpf(10) ** (-mp.mp.dps // 2) for z in rts):
            raise ValueError("Delta does not have four real roots (not subextremal)")
        self.rn, self.rm, self.rp, self.rc = sorted(mp.re(z) for z in rts)
        self.sing = [self.rn, self.rm, self.rp, self.rc]
        self.kap = {k: abs(self.dDelta(x)) / (2 * x * x) for k, x in zip("n-+c", self.sing)}
        self.km, self.kp, self.kc = self.kap["-"], self.kap["+"], self.kap["c"]
        self.A = [x * x / self.dDelta(x) for x in self.sing]
        self.cp = sum(self.A[i] * mp.log(abs(self.rp - self.sing[i])) for i in (0, 1, 3))
        self.cc = sum(self.A[i] * mp.log(abs(self.rc - self.sing[i])) for i in (0, 1, 2))
        disc = 9 - 8 * self.Q2
        self.rP = (3 + mp.sqrt(disc)) / 2
        FP = self.F(self.rP)
        self.OmegaP = mp.sqrt(FP) / self.rP
        self.lamP = mp.sqrt(FP * (2 - 3 / self.rP)) / self.rP
        self.gamma0 = self.lamP / 2

    # -- constructors ---------------------------------------------------------------------------
    @staticmethod
    def extremal_roots(Lam):
        """Positive roots of 2 Lam r^3/3 - r + 1 = 0 (double-root condition), ascending."""
        Lam = mp.mpf(Lam)
        cub = mp.polyroots([2 * Lam / 3, 0, -1, 1], maxsteps=800, extraprec=4 * mp.mp.prec)
        return sorted(mp.re(z) for z in cub if abs(mp.im(z)) < mp.mpf(10) ** (-mp.mp.dps // 2) and mp.re(z) > 0)

    @classmethod
    def from_qfrac(cls, Lam, qfrac):
        """Q = qfrac * Q_max(Lam) (FORMULATION section 1.3 for Lam = 3/50)."""
        Lam = mp.mpf(Lam)
        r_ext = cls.extremal_roots(Lam)[0]
        Qmax2 = r_ext - Lam * r_ext ** 4 / 3
        g = cls(Lam, mp.mpf(qfrac) ** 2 * Qmax2)
        g.Qmax2, g.qfrac = Qmax2, mp.mpf(qfrac)
        return g

    @classmethod
    def frozen(cls, name):
        """P_A / P_B from the exact rationals (qfrac 249/250, 999/1000; Lambda = 3/50)."""
        qf = {"PA": mp.mpf(249) / 250, "PB": mp.mpf(999) / 1000}[name]
        return cls.from_qfrac(mp.mpf(3) / 50, qf)

    # -- functions ------------------------------------------------------------------------------
    def Delta(self, r):
        return mp.polyval(self.dco[::-1], r)

    def dDelta(self, r):
        return -4 * self.Lam / 3 * r ** 3 + 2 * r - 2

    def F(self, r):
        return self.Delta(r) / r ** 2

    def V(self, r, l):
        Fr = (self.dDelta(r) * r - 2 * self.Delta(r)) / r ** 3
        return self.F(r) * (l * (l + 1) / r ** 2 + Fr / r)

    def shifted(self, c):
        """Coefficients of Delta(c + t) in t (ascending)."""
        out, pc = [], self.dco[::-1]
        for k in range(5):
            out.append(mp.polyval(pc, c) / mp.factorial(k))
            pc = [coef * (len(pc) - 1 - i) for i, coef in enumerate(pc[:-1])]
        return out


def _polymul(x, y):
    out = [mp.mpc(0)] * (len(x) + len(y) - 1)
    for i, xi in enumerate(x):
        if xi == 0:
            continue
        for j, yj in enumerate(y):
            out[i + j] += xi * yj
    return out


def ode_polys(g, c, w, L):
    """A = Delta^2, B = Delta Delta', C = w^2 r^4 - L Delta at r = c + t (ascending in t)."""
    d = g.shifted(c)
    dp = [d[k + 1] * (k + 1) for k in range(4)]
    A = _polymul(d, d)
    B = _polymul(d, dp)
    r4 = [c ** 4, 4 * c ** 3, 6 * c ** 2, 4 * c, mp.mpf(1)]
    C = [w * w * r4[k] - L * d[k] for k in range(5)]
    return A, B, C


def frobenius(g, s, rho, w, L, t, singular=True, nmax=20000):
    """R, dR/dt at r = s + t for the Frobenius solution |t|^rho sum a_k t^k (a_0 = 1)."""
    A, B, C = ode_polys(g, s, w, L)
    if singular:
        A[0] = A[1] = B[0] = mp.mpf(0)          # Delta(s) = 0 exactly
    a = [mp.mpc(1)]
    S, Sd = mp.mpc(1), rho / t
    tk = mp.mpf(1)
    eps = mp.mpf(10) ** (-mp.mp.dps + 2)
    small, peak = 0, mp.mpf(1)
    for n in range(1, nmax):
        acc = mp.mpc(0)
        for k in range(max(0, n - 6), n):
            j2, j1, j0 = n - k + 2, n - k + 1, n - k
            term = 0
            if j2 < len(A):
                term += A[j2] * (k + rho) * (k + rho - 1)
            if j1 < len(B):
                term += B[j1] * (k + rho)
            if j0 < len(C):
                term += C[j0]
            acc += a[k] * term
        ind = A[2] * n * (n + 2 * rho)               # = A2 (n+rho)(n+rho-1) + B1 (n+rho) + C0
        an = -acc / ind
        a.append(an)
        tk *= t
        tv = an * tk
        S += tv
        Sd += an * (n + rho) * tk / t
        peak = max(peak, abs(tv))
        if abs(tv) < eps * abs(S):
            small += 1
            if small > 8:
                break
        else:
            small = 0
    pref = abs(t) ** rho
    return pref * S, pref * Sd, peak


def taylor_step(g, c, R, dR, w, L, d, nmax=20000):
    """Continue (R, R') from the ordinary point c to c + d."""
    A, B, C = ode_polys(g, c, w, L)
    bco = [R, dR]
    S, Sd = R + dR * d, dR
    eps = mp.mpf(10) ** (-mp.mp.dps + 2)
    small = 0
    dk = d
    peak = abs(S)
    for n in range(nmax):
        acc = mp.mpc(0)
        for j in range(1, min(len(A), n + 3)):
            m = n - j + 2
            if m >= 0:
                acc += A[j] * m * (m - 1) * bco[m]
        for j in range(min(len(B), n + 2)):
            m = n - j + 1
            if m >= 0:
                acc += B[j] * m * bco[m]
        for j in range(min(len(C), n + 1)):
            acc += C[j] * bco[n - j]
        bn2 = -acc / (A[0] * (n + 2) * (n + 1))
        bco.append(bn2)
        dk *= d
        tv = bn2 * dk
        S += tv
        Sd += (n + 2) * bn2 * dk / d
        peak = max(peak, abs(tv))
        if abs(tv) < eps * (abs(S) + abs(Sd * d)):
            small += 1
            if small > 8:
                break
        else:
            small = 0
    return S, Sd, peak


class Solver:
    """Wronskian evaluator with step-size controls (K, frac); dps is the caller's mp.dps."""

    def __init__(self, g: Geom, K=6, frac="0.4", rmatch=None):
        self.g, self.K, self.frac = g, mp.mpf(K), mp.mpf(frac)
        self.rmatch = g.rP if rmatch is None else mp.mpf(rmatch)
        self.max_peak = mp.mpf(1)       # largest |term|/|sum| seen: digits lost to cancellation

    def kloc(self, r, w, L):
        return abs(w) * r * r / abs(self.g.Delta(r)) + mp.sqrt(L) / r + 1

    def continue_to(self, r0, R, dR, w, L, r1):
        g = self.g
        r = r0
        while abs(r1 - r) > mp.mpf(10) ** (-mp.mp.dps + 5):
            dist = min(abs(r - s) for s in g.sing)
            step = min(self.frac * dist, self.K / self.kloc(r, w, L))
            d = r1 - r
            if abs(d) > step:
                d = step if d > 0 else -step
            R, dR, pk = taylor_step(g, r, R, dR, w, L, d)
            self.max_peak = max(self.max_peak, pk / (abs(R) + abs(dR * d)))
            r = r + d
        return R, dR

    def wronskian(self, w, l, jost=True):
        g = self.g
        L = l * (l + 1)
        w = mp.mpc(w)
        rho_p, rho_c = -1j * w / (2 * g.kp), -1j * w / (2 * g.kc)
        dpm = min(g.rp - g.rm, g.rc - g.rp)
        tp = min(mp.mpf("0.25") * dpm, self.K * dpm / max(1, abs(rho_p)) / 4,
                 self.K / self.kloc(g.rp + mp.mpf("0.25") * dpm, w, L) / 4)
        tc = min(mp.mpf("0.25") * (g.rc - g.rp), self.K * (g.rc - g.rp) / max(1, abs(rho_c)) / 4)
        Rp, dRp, pk1 = frobenius(g, g.rp, rho_p, w, L, tp)
        Rp, dRp = self.continue_to(g.rp + tp, Rp, dRp, w, L, self.rmatch)
        Rc, dRc, pk2 = frobenius(g, g.rc, rho_c, w, L, -tc)
        Rc, dRc = self.continue_to(g.rc - tc, Rc, dRc, w, L, self.rmatch)
        self.max_peak = max(self.max_peak, pk1, pk2)
        W = g.Delta(self.rmatch) * (Rp * dRc - dRp * Rc)
        if jost:
            W = W * mp.exp(-1j * w * g.cp + 1j * w * g.cc) / (g.rp * g.rc)
        return W

    def find(self, l, seed, tol=None, maxsteps=60):
        """Secant iteration; converged when the relative step is below tol (scale-free)."""
        f = lambda x: self.wronskian(x, l)
        x0 = mp.mpc(seed)
        x1 = x0 * (1 + mp.mpf("1e-6")) + mp.mpc(0, "1e-9")
        tol = mp.mpf(10) ** (-mp.mp.dps + 12) if tol is None else tol
        f0, f1 = f(x0), f(x1)
        for _ in range(maxsteps):
            if f1 == f0:
                break
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
            x0, f0, x1 = x1, f1, x2
            f1 = f(x1)
            if abs(x1 - x0) <= tol * abs(x1):
                return x1
        raise ValueError(f"secant did not converge (last relative step {mp.nstr(abs(x1 - x0) / abs(x1), 3)})")


def qnm(g: Geom, l, seed_over_km, dps=50, K=6, check=True):
    """QNM omega/kappa_- from a seed (in kappa_- units), with an internal consistency check
    (recomputed at dps + 20 and K/2: returns (value, |difference|))."""
    with mp.workdps(dps):
        s = Solver(g, K=K)
        w1 = s.find(l, mp.mpc(seed_over_km) * g.km) / g.km
    if not check:
        return w1, None
    with mp.workdps(dps + 20):
        s2 = Solver(g, K=mp.mpf(K) / 2)
        w2 = s2.find(l, w1 * g.km) / g.km
    return w2, abs(w2 - w1)
