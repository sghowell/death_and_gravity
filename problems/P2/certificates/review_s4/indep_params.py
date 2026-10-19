"""Independent P_B parameters from the exact definitions (FORMULATION 1.3, D6), mpmath only.

No import from the p2 package.  Used by the reviewer's S4 checks.
"""
from __future__ import annotations

import mpmath as mp

mp.mp.dps = 50


class PB:
    def __init__(self, dps=50):
        with mp.workdps(dps):
            self.Lam = mp.mpf(3) / 50
            # r_ext: smallest positive root of r^3 - 25 r + 25 = 0
            rr = sorted([r.real for r in mp.polyroots([1, 0, -25, 25], maxsteps=200, extraprec=200) if abs(r.imag) < 1e-30])
            self.r_ext = [r for r in rr if r > 0][0]
            self.Qmax2 = self.r_ext - self.r_ext ** 4 / 50
            self.Q2 = (mp.mpf(999) / 1000) ** 2 * self.Qmax2
            co = [-self.Lam / 3, 0, 1, -2, self.Q2]          # Delta, descending
            roots = sorted([r.real for r in mp.polyroots(co, maxsteps=400, extraprec=400)])
            self.rn, self.rm, self.rp, self.rc = roots
            self.roots = roots
            self.dco = co
            self.A = [r ** 2 / self.dDelta(r) for r in roots]
            self.km = abs(self.dDelta(self.rm)) / (2 * self.rm ** 2)
            self.kp = abs(self.dDelta(self.rp)) / (2 * self.rp ** 2)
            self.kc = abs(self.dDelta(self.rc)) / (2 * self.rc ** 2)
            self.rP = (3 + mp.sqrt(9 - 8 * self.Q2)) / 2
            FP = self.F(self.rP)
            self.OmegaP = mp.sqrt(FP) / self.rP
            self.gamma0 = mp.sqrt(FP * (2 - 3 / self.rP)) / self.rP / 2
            self.d = (mp.mpf(1) / 2 + mp.mpf(1) / 20) * self.km      # d_delta
            self.cp = sum(a * mp.log(abs(self.rp - r)) for a, r in zip(self.A, roots) if r != self.rp)
            self.cc = sum(a * mp.log(abs(self.rc - r)) for a, r in zip(self.A, roots) if r != self.rc)
            self.xP = self.rstar(self.rP)
            self.X0 = sum(a * mp.log(abs(r)) for a, r in zip(self.A, roots))   # Re r*(0)
            self.astar = 2 * self.OmegaP * self.gamma0

    def Delta(self, r):
        return -self.Lam / 3 * r ** 4 + r ** 2 - 2 * r + self.Q2

    def dDelta(self, r):
        return -4 * self.Lam / 3 * r ** 3 + 2 * r - 2

    def F(self, r):
        return self.Delta(r) / r ** 2

    def U(self, r):
        return self.Delta(r) / r ** 4

    def Wt(self, r):
        D, Dp = self.Delta(r), self.dDelta(r)
        F = D / r ** 2
        Fr = Dp / r ** 2 - 2 * D / r ** 3
        return F * Fr / r - F / r ** 2 / 4

    def V(self, r, l):
        F = self.F(r)
        Fr = self.dDelta(r) / r ** 2 - 2 * self.Delta(r) / r ** 3
        return F * (l * (l + 1) / r ** 2 + Fr / r)

    def rstar(self, r):
        """Re r*(r) = sum A_i log|r - r_i|."""
        return sum(a * mp.log(abs(r - ri)) for a, ri in zip(self.A, self.roots))

    def r_of_x(self, X):
        """Real r in (r_+, r_c) with r*(r) = X (bisection + Newton, mpmath)."""
        lo, hi = self.rp + mp.mpf(10) ** (-40), self.rc - mp.mpf(10) ** (-40)
        for _ in range(400):
            m = (lo + hi) / 2
            if self.rstar(m) > X:
                hi = m
            else:
                lo = m
            if hi - lo < mp.mpf(10) ** (-mp.mp.dps + 5):
                break
        r = (lo + hi) / 2
        for _ in range(5):
            r = r - (self.rstar(r) - X) * self.F(r)
        return r


if __name__ == "__main__":
    p = PB()
    for k in ("rn", "rm", "rp", "rc", "km", "kp", "kc", "rP", "OmegaP", "gamma0", "d", "cp", "cc", "xP", "X0", "astar"):
        print(k, mp.nstr(getattr(p, k), 20))
    print("sum A", mp.nstr(sum(p.A), 5))
    print("Omega_P/km", mp.nstr(p.OmegaP / p.km, 15), "gamma0/km", mp.nstr(p.gamma0 / p.km, 15))
