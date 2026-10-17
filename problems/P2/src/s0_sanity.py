"""P2 S0 sanity numerics -- floating-point / mpmath only, NOTHING here is certified.

Purpose (FORMULATION.md sections 1.3 and 5; notes/s0-gate.md sections 1 and 5):
  * compute the geometric constants of the frozen parameter points
    (horizon radii, surface gravities, photon sphere, trapping gap gamma_0);
  * reproduce a few of the published quasinormal modes of
    Cardoso-Costa-Destounis-Hintz-Jansen, arXiv:1711.10502, Table I
    (Lambda M^2 = 0.06, Q/Q_max = 0.996) by an independent method:
    Frobenius series at r_+ and r_c, Taylor-series analytic continuation
    of the radial ODE (polynomial coefficients, D-finite), zeros of the
    Jost-normalised Wronskian;
  * sample the Wronskian on the bottom edge of the target strip to size
    the validated argument-principle computation (part (iii)).

Radial equation (massless minimally coupled scalar, phi = e^{-i w t} R(r) Y_lm):
    (Delta R')' + (w^2 r^4 / Delta - L) R = 0,
    Delta(r) = r^2 F(r) = -(Lam/3) r^4 + r^2 - 2 M r + Q^2,  L = l(l+1).
QNM boundary conditions: R ~ (r - r_+)^{-i w/(2 k_+)} at r_+,
R ~ (r_c - r)^{-i w/(2 k_c)} at r_c (both analytic factors otherwise).

Run:  .venv/bin/python problems/P2/src/s0_sanity.py [all|geom|modes|rate|winding]
      (geom <1 s, modes ~20 s, rate ~15 s, winding ~95 s on one core)
"""
from __future__ import annotations

import itertools
import sys
import time

import mpmath as mp

mp.mp.dps = 40


# ---------------------------------------------------------------- geometry
class RNdS:
    """Subextremal RNdS with M = 1 units unless M is given."""

    def __init__(self, Lam, qfrac=None, Q=None, M=1):
        self.M = mp.mpf(M)
        self.Lam = mp.mpf(Lam)
        # extremal charge: double root r0 of Delta solves 2 Lam r^3/3 - r + M = 0
        # (smallest positive root, continuous with r0 = M at Lam = 0)
        cub = mp.polyroots([2 * self.Lam / 3, 0, -1, self.M], maxsteps=200, extraprec=100)
        pos = sorted(mp.re(z) for z in cub if abs(mp.im(z)) < mp.mpf(10) ** (-30) and mp.re(z) > 0)
        self.r_ext = pos[0]
        self.Qmax = mp.sqrt(self.M * self.r_ext - self.Lam * self.r_ext ** 4 / 3)
        if Q is None:
            Q = mp.mpf(qfrac) * self.Qmax
        self.Q = mp.mpf(Q)
        self.qfrac = self.Q / self.Qmax
        # Delta coefficients, ascending powers of r
        self.dco = [self.Q ** 2, -2 * self.M, mp.mpf(1), mp.mpf(0), -self.Lam / 3]
        roots = mp.polyroots(self.dco[::-1], maxsteps=400, extraprec=200)
        rr = sorted(mp.re(z) for z in roots)
        self.rn, self.rm, self.rp, self.rc = rr
        self.sing = rr
        self.kap = {name: abs(self.dDelta(r)) / (2 * r ** 2)
                    for name, r in zip(("n", "-", "+", "c"), rr)}
        # partial fractions 1/F = r^2/Delta = sum A_i/(r - r_i)
        self.A = [r ** 2 / self.dDelta(r) for r in rr]
        # photon sphere (Lambda independent): r^2 - 3 M r + 2 Q^2 = 0
        self.rP = (3 * self.M + mp.sqrt(9 * self.M ** 2 - 8 * self.Q ** 2)) / 2
        FP = self.F(self.rP)
        self.OmegaP = mp.sqrt(FP) / self.rP                     # eikonal Re w / (l + 1/2)
        self.lamP = mp.sqrt(FP * (2 - 3 * self.M / self.rP)) / self.rP   # Lyapunov exponent
        self.gamma0 = self.lamP / 2                              # Hintz-Vasy essential gap

    def Delta(self, r):
        return mp.polyval(self.dco[::-1], r)

    def dDelta(self, r):
        return -4 * self.Lam / 3 * r ** 3 + 2 * r - 2 * self.M

    def F(self, r):
        return self.Delta(r) / r ** 2

    def shifted(self, c):
        """Coefficients (ascending) of Delta(c + t) in t."""
        out = []
        p = self.dco[::-1]
        for k in range(5):
            out.append(mp.polyval(p, c) / mp.factorial(k))
            p = [coef * (len(p) - 1 - i) for i, coef in enumerate(p[:-1])]
        return out

    def report(self):
        k = self.kap
        km = k["-"]
        rows = [
            ("Lambda M^2", self.Lam * self.M ** 2),
            ("Q/Q_max", self.qfrac),
            ("Q/M", self.Q / self.M),
            ("Q_max/M", self.Qmax / self.M),
            ("r_n", self.rn), ("r_-", self.rm), ("r_+", self.rp), ("r_c", self.rc),
            ("kappa_-", km), ("kappa_+", k["+"]), ("kappa_c", k["c"]),
            ("kappa_c^dS = sqrt(Lam/3)", mp.sqrt(self.Lam / 3)),
            ("kappa_+/kappa_-", k["+"] / km), ("kappa_c/kappa_-", k["c"] / km),
            ("r_P", self.rP), ("Omega_P", self.OmegaP), ("Omega_P/kappa_-", self.OmegaP / km),
            ("lambda_P (Lyapunov)", self.lamP), ("gamma0 = lambda_P/2", self.gamma0),
            ("gamma0/kappa_- (eikonal beta_PS)", self.gamma0 / km),
        ]
        for name, val in rows:
            print(f"  {name:34s} {mp.nstr(val, 12)}")


def polymul(a, b):
    out = [mp.mpc(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def ode_polys(bh, c, w, L):
    """A = Delta^2, B = Delta Delta', C = w^2 r^4 - L Delta, expanded at r = c + t."""
    d = bh.shifted(c)
    dp = [d[k + 1] * (k + 1) for k in range(4)]
    A = polymul(d, d)
    B = polymul(d, dp)
    r4 = [c ** 4, 4 * c ** 3, 6 * c ** 2, 4 * c, mp.mpf(1)]
    C = [w ** 2 * r4[k] - L * d[k] for k in range(5)]
    return A, B, C


def frobenius(bh, s, rho, w, L, t, nmax=4000):
    """R and dR/dt at r = s + t for the exponent-rho Frobenius solution at singular s.
    R = |t|^rho * sum a_k t^k (a_0 = 1); the |t|^rho convention handles t < 0."""
    A, B, C = ode_polys(bh, s, w, L)
    a = [mp.mpc(1)]
    S = mp.mpc(1)
    Sd = rho / t
    tk = mp.mpf(1)
    small = 0
    eps = mp.mpf(10) ** (-mp.mp.dps + 3)
    for n in range(1, nmax):
        acc = mp.mpc(0)
        for k in range(max(0, n - 8), n):
            j2, j1, j0 = n - k + 2, n - k + 1, n - k
            term = 0
            if j2 < len(A):
                term += A[j2] * (k + rho) * (k + rho - 1)
            if j1 < len(B):
                term += B[j1] * (k + rho)
            if j0 < len(C):
                term += C[j0]
            acc += a[k] * term
        ind = A[2] * (n + rho) * (n + rho - 1) + B[1] * (n + rho) + C[0]
        an = -acc / ind
        a.append(an)
        tk *= t
        term_v = an * tk
        S += term_v
        Sd += an * (n + rho) * tk / t
        if abs(term_v) < eps * abs(S):
            small += 1
            if small > 6:
                break
        else:
            small = 0
    pref = abs(t) ** rho
    return pref * S, pref * Sd


def taylor_step(bh, c, R, dR, w, L, d, nmax=4000):
    """Continue (R, R') from ordinary point c to c + d by the Taylor recurrence."""
    A, B, C = ode_polys(bh, c, w, L)
    b = [R, dR]
    S = R + dR * d
    Sd = dR
    eps = mp.mpf(10) ** (-mp.mp.dps + 3)
    small = 0
    dk = d  # d^(n+1) bookkeeping below
    for n in range(nmax):
        # coefficient of t^n: sum_j A_j (n-j+2)(n-j+1) b_{n-j+2} + B_j (n-j+1) b_{n-j+1} + C_j b_{n-j}
        acc = mp.mpc(0)
        for j in range(1, len(A)):
            m = n - j + 2
            if 0 <= m < len(b):
                acc += A[j] * m * (m - 1) * b[m]
        for j in range(len(B)):
            m = n - j + 1
            if 0 <= m < len(b):
                acc += B[j] * m * b[m]
        for j in range(len(C)):
            m = n - j
            if 0 <= m < len(b):
                acc += C[j] * b[m]
        bn2 = -acc / (A[0] * (n + 2) * (n + 1))
        b.append(bn2)
        dk = dk * d  # d^(n+2)
        tv = bn2 * dk
        S += tv
        Sd += (n + 2) * bn2 * dk / d
        if abs(tv) < eps * abs(S):
            small += 1
            if small > 6:
                break
        else:
            small = 0
    return S, Sd


def continue_to(bh, r0, R, dR, w, L, r1, frac=None):
    frac = mp.mpf("0.5") if frac is None else frac
    r = r0
    while abs(r1 - r) > mp.mpf(10) ** (-30):
        dist = min(abs(r - s) for s in bh.sing)
        step = frac * dist
        d = r1 - r
        if abs(d) > step:
            d = step if d > 0 else -step
        R, dR = taylor_step(bh, r, R, dR, w, L, d)
        r = r + d
    return R, dR


def wronskian(bh, w, l, rmatch=None, jost=True):
    """Jost-normalised Wronskian W = Delta (R_+ R_c' - R_+' R_c) at rmatch."""
    L = l * (l + 1)
    w = mp.mpc(w)
    kp, kc = bh.kap["+"], bh.kap["c"]
    rho_p = -1j * w / (2 * kp)
    rho_c = -1j * w / (2 * kc)
    if rmatch is None:
        rmatch = bh.rP
    hp = mp.mpf("0.4") * min(bh.rp - bh.rm, bh.rc - bh.rp)
    hc = mp.mpf("0.4") * (bh.rc - bh.rp)
    Rp, dRp = frobenius(bh, bh.rp, rho_p, w, L, hp)
    Rp, dRp = continue_to(bh, bh.rp + hp, Rp, dRp, w, L, rmatch)
    Rc, dRc = frobenius(bh, bh.rc, rho_c, w, L, -hc)
    Rc, dRc = continue_to(bh, bh.rc - hc, Rc, dRc, w, L, rmatch)
    W = bh.Delta(rmatch) * (Rp * dRc - dRp * Rc)
    if jost:
        rr = bh.sing
        cp = sum(bh.A[i] * mp.log(abs(bh.rp - rr[i])) for i in range(4) if i != 2)
        cc = sum(bh.A[i] * mp.log(abs(bh.rc - rr[i])) for i in range(4) if i != 3)
        W = W * mp.exp(-1j * w * cp) * mp.exp(1j * w * cc) / (bh.rp * bh.rc)
    return W


def find_qnm(bh, l, guess):
    f = lambda w: wronskian(bh, w, l)
    g = mp.mpc(guess)
    return mp.findroot(f, (g, g * (1 + mp.mpf("1e-3")) + mp.mpc(0, "1e-6")), solver="secant",
                       tol=mp.mpf(10) ** (-28), maxsteps=80)


# ---------------------------------------------------------------- drivers
TABLE_I = [  # CCDHJ arXiv:1711.10502 Table I, omega_0 / kappa_-, Lambda M^2 = 0.06, Q/Q_max = 0.996
    (0, mp.mpc(0, "-0.8539013779"), "NE l=0"),
    (1, mp.mpc(0, "-1.5003853731"), "dS l=1"),
    (1, mp.mpc("3.2426164126", "-0.7958326323"), "PS l=1"),
    (2, mp.mpc("5.4796067815", "-0.7754179185"), "PS l=2"),
    (5, mp.mpc("12.1322536641", "-0.7665731858"), "PS l=5"),
    (10, mp.mpc("23.1897597770", "-0.7649242108"), "PS l=10"),
]


# Frozen points (FORMULATION.md section 1.3).  M = 1, Lambda = 3/50.
POINT_A = {"Lam": "0.06", "qfrac": "0.996"}   # P_A: CCDHJ Table I point (second target)
POINT_B = {"Lam": "0.06", "qfrac": "0.999"}   # P_B: near-extremal window of CCDHJ Fig. 3 (primary)

# P_B float spectrum probe: (l, seed in units of kappa_- unless flagged, label)
SEEDS_B = [
    (0, mp.mpc(0, "-0.9"), "NE l=0 n=0"),
    (0, mp.mpc(0, "-2.0"), "NE l=0 n=1"),
    (1, mp.mpc(0, "-1.9"), "NE l=1 n=0"),
    (2, mp.mpc(0, "-3.2"), "NE l=2 n=0"),
    (1, mp.mpc(0, "-3.3"), "dS l=1"),
    (1, mp.mpc("7.2", "-1.75"), "PS l=1"),
    (2, mp.mpc("12.2", "-1.71"), "PS l=2"),
    (10, mp.mpc("51.6", "-1.68"), "PS l=10"),
]


def run_geom():
    for name, pt in (("A", POINT_A), ("B", POINT_B)):
        bh = RNdS(**pt)
        print(f"point {name}: {pt}")
        bh.report()
        # Hintz 2504.01734 Thm 1.1 near-horizon prediction (r_e, r_c frozen at current radii)
        re_, rc_ = bh.rp, bh.rc
        vk = (rc_ ** 2 + 2 * re_ * rc_ - 3 * re_ ** 2) / (rc_ ** 2 + 2 * re_ * rc_ + 3 * re_ ** 2)
        lam = [(1 + mp.sqrt(1 + 4 * l * (l + 1) / vk)) / 2 for l in range(3)]
        print(f"  varkappa (Hintz Thm 1.1)           {mp.nstr(vk, 10)}; lambda_l^+ (l=0,1,2) = "
              + ", ".join(mp.nstr(x, 8) for x in lam))
        for delta, tag in ((mp.mpf("0.01"), "1/100"), (mp.mpf("0.05"), "1/20")):
            print(f"  strip edge (1/2+{tag}) kappa_- as fraction of gamma0: "
                  f"{mp.nstr((mp.mpf('0.5') + delta) * bh.kap['-'] / bh.gamma0, 4)}")


def run_modes():
    bh = RNdS(**POINT_A)
    km = bh.kap["-"]
    print("point A modes (w/kappa_-), independent Frobenius + Taylor continuation vs CCDHJ Table I:")
    for l, ref, label in TABLE_I:
        t0 = time.time()
        w = find_qnm(bh, l, ref * km)
        val = w / km
        print(f"  {label:8s} computed {mp.nstr(val, 14):>38s}   CCDHJ {mp.nstr(ref, 11):>30s}"
              f"   |diff| {mp.nstr(abs(val - ref), 3)}   ({time.time() - t0:.1f}s)", flush=True)
    ref100 = mp.mpc("222.0602900249", "-0.7643094446")
    print(f"  eikonal (l+1/2)Omega_P/kappa_- at l=100: {mp.nstr(mp.mpf('100.5') * bh.OmegaP / km, 10)}"
          f"  vs Re {mp.nstr(mp.re(ref100), 13)};  gamma0/kappa_- {mp.nstr(bh.gamma0 / km, 10)}"
          f"  vs -Im {mp.nstr(-mp.im(ref100), 10)}")
    bhB = RNdS(**POINT_B)
    kB = bhB.kap["-"]
    print("point B float spectrum (w/kappa_-), same method (no published table at this point):")
    for l, seed, label in SEEDS_B:
        t0 = time.time()
        w = find_qnm(bhB, l, seed * kB)
        print(f"  {label:11s} l={l:2d}: {mp.nstr(w / kB, 12):>34s}   ({time.time() - t0:.1f}s)", flush=True)


def run_rate(ls=(2, 10, 20), nsamp=12):
    """Local phase rate |d arg(W/w)/d Re w| [rad per kappa_-] on Im w = -0.55 kappa_-."""
    mp.mp.dps = 30
    for name, pt in (("A", POINT_A), ("B", POINT_B)):
        bh = RNdS(**pt)
        km = bh.kap["-"]
        for l in ls:
            top = (mp.mpf("1.25") * (l + mp.mpf("0.5")) * bh.OmegaP / km) + 2
            worst = 0
            for i in range(nsamp):
                x = top * (i + mp.mpf("0.5")) / nsamp
                h = mp.mpf("1e-4")
                w0 = mp.mpc(x, "-0.55") * km
                w1 = mp.mpc(x + h, "-0.55") * km
                rate = mp.im(mp.log((wronskian(bh, w1, l) / w1) / (wronskian(bh, w0, l) / w0))) / h
                worst = max(worst, abs(rate))
            print(f"  point {name} l={l:3d}: max |d arg(W/w)/d Re w| over {nsamp} samples = "
                  f"{mp.nstr(worst, 3)} rad per kappa_-", flush=True)
    mp.mp.dps = 40


def _winding(bh, l, xl, xr, ylo, yhi, h0):
    """Float winding number of W_l around [xl,xr] x [ylo,yhi] (kappa_- units), adaptive arg steps."""
    km = bh.kap["-"]
    corners = [mp.mpc(xl, ylo), mp.mpc(xr, ylo), mp.mpc(xr, yhi), mp.mpc(xl, yhi), mp.mpc(xl, ylo)]
    f = lambda t: wronskian(bh, t * km, l)
    total, nev, minabs = mp.mpf(0), 0, mp.inf
    for a, b in itertools.pairwise(corners):
        n = max(2, int(abs(b - a) / h0) + 1)
        pts = [a + (b - a) * k / n for k in range(n + 1)]
        vals = [f(p) for p in pts]
        nev += len(pts)
        k = 0
        while k < len(pts) - 1:
            d = mp.im(mp.log(vals[k + 1] / vals[k]))
            if abs(d) > 0.6 and abs(pts[k + 1] - pts[k]) > 1e-6:
                mid = (pts[k] + pts[k + 1]) / 2
                pts.insert(k + 1, mid)
                vals.insert(k + 1, f(mid))
                nev += 1
                continue
            total += d
            if pts[k] != 0:
                minabs = min(minabs, abs(vals[k] / (pts[k] * km)))
            k += 1
    return total / (2 * mp.pi), nev, minabs


def run_winding(ls=(0, 1, 2, 5)):
    """Float (non-rigorous) zero counts in the part-(iii) boxes: expect 1 for l=0 (w=0), else 0."""
    mp.mp.dps = 25
    for name, pt, h0 in (("B", POINT_B, 0.25), ("A", POINT_A, 0.15)):
        bh = RNdS(**pt)
        km = bh.kap["-"]
        for l in ls:
            t0 = time.time()
            xr = mp.mpf("1.25") * (l + mp.mpf("0.5")) * bh.OmegaP / km + 2
            xl = -xr if l == 0 else mp.mpf(0)
            wn, nev, mn = _winding(bh, l, xl, xr, mp.mpf("-0.55"), mp.mpf("0.1"), h0)
            print(f"  point {name} l={l}: box [{mp.nstr(xl, 4)},{mp.nstr(xr, 4)}] x [-0.55,0.1] kappa_-: "
                  f"winding {mp.nstr(wn, 4)} ({nev} evals, min|W/w| {mp.nstr(mn, 3)}, "
                  f"{time.time() - t0:.0f}s)", flush=True)
    mp.mp.dps = 40


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    t = time.time()
    if what in ("all", "geom"):
        run_geom()
    if what in ("all", "modes"):
        run_modes()
    if what in ("all", "rate"):
        run_rate()
    if what in ("all", "winding"):
        run_winding()
    print(f"total {time.time() - t:.1f}s")
