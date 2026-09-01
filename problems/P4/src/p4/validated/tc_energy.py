"""Theorem C, S4-2 (route b1): certified direct-energy exclusion  "no eigenvalues with
Re kappa >= c0"  for the reduced (chi, eta) problem of tc_reduce (see notes/s4-energy.md).

Mechanism (GHJS 2509.12435, Lemmas 3.7-3.9, transferred).  The scalar chi-equation
chi'' - T chi' - U chi = 0 is put in Liouville normal form v'' = Q(x; kappa) v via
v = chi/f, 2 f'/f = T; then for an eigenfunction (class of s2-theorem-b.md section 3.8)
the energy identity  int |v'|^2 + int Q |v|^2 = 0  holds on (-inf, 0) (boundary terms
vanish: notes/s4-energy.md section 2), and
    Im Q = tau (2 eps q2 + q1) + Im R,   Q(eps) = q2 eps^2 + q1 eps + q0 + R(eps),
with Q = QN/QD exactly rational in kappa (real coefficients), q_j the kappa-quadratic
quotient of the division QN = (q2 k^2 + q1 k + q0) QD + rem, R = rem/QD.  The certificate
checks, on boxes covering (-inf, 0) (A3/A2 tube + certified sonic/centre series):
    (a) q2 > 0;  (b) 2 c0 q2 + q1 - B1 > 0;  (c) c0^2 q2 + c0 q1 + q0 - B0 > 0;
    (d) all poles of Q(x; .) have Re kappa <= rbar < c0 (QD-factor roots: A(x), -F_N(x)
        and the Bc-cubic via Routh-Hurwitz),
where B0 >= sup_{Re kappa >= c0} |R|, B1 >= sup |dR/dkappa| (maximum modulus: sup on the
line Re kappa = c0, tau-boxes + coefficient tail; B1 also by the Cauchy shift).  (a)-(d)
for all x  =>  Im Q/tau > 0 pointwise (kills Im kappa != 0) and Q(eps) > 0 for real
eps >= c0 (kills the rest); chi = 0 then forces p = const * g and the sonic gauge
N_p(0) = 0 forces const = 0 (notes/s4-energy.md section 1).

Exact algebra used at the sonic end (all verified in fmpq_mpoly at import):
    QN = Dpoly^10 QN_red,  QD = SD^3 D2^2 Bc^2 = Dpoly^13 QD_red,  TN1 = Dpoly^5 TN1_red,
    SD = Dpoly SD_red, D2 = Dpoly^3 D2_red, Bc = Dpoly^2 Bc_red,
and the field identity QN_red(u0) = 0 in Q(V0, sqrt3) at the closed-form sonic point
(checked once, scripts record it; re-run with mode ``field``).  Hence on [-dc, 0]:
    x^2 Q = QQ / (E_D^3 QD_red),   QN_red = x*QQ (mean value, exact zero at x = 0),
    E_D = D(x)/x = hull(D'), all factors nonvanishing balls: interval conditions hold
    up to and including x = 0.  Known answer: x^2 Q -> (sigma(kappa)^2 - 1)/4.

Certified region: x in [-4.5, 0] (tube steps + sonic ladder + closing box at the
sonic point); the deep tail x < -4.5 is Lemma D of notes/s4-energy.md section 4
(analytic, float-validated; the identified closing computation is exact reduction
modulo the centre fixed-point ideal (2NV+1) -- S4-3).  Working threshold c0 = 6.

Run:  PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_energy all <tube.json> [c0]
      (modes: tube | sonic | centre | all | field;  writes results/theorem_c/tc_energy.json)
"""
from __future__ import annotations

import json
import os
import sys
import time

from flint import acb, arb, fmpq

from .arbseries import precision, to_arb
from .tc_reduce import Reduction

KD = 12                      # kappa-degree of QN
RBAR = 2.0                   # certified upper bound for Re(poles); needs c0 > RBAR
DELC = 1.2                   # Cauchy shift for B1; needs RBAR < c0 - DELC
TTAIL = 64.0                 # line-sup switches to the coefficient tail at |tau| = TTAIL


# ------------------------------------------------------------------ exact atom algebra
def _dpow(f, Dp, expect=None):
    m = 0
    while True:
        q, r = divmod(f, Dp)
        if not r.is_zero():
            break
        f, m = q, m + 1
    if expect is not None and m != expect:
        raise RuntimeError(f"Dpoly-power {m} != expected {expect}")
    return f, m


class Atoms:
    """Exact mpoly data + compiled term lists.  full=True adds the QN_red/TN1_red layer."""

    def __init__(self, full=True):
        red = Reduction()
        rp = red.reduce_pair()
        Ac, Bc = rp["chi"][0], rp["chi"][1]
        Ae, Be = rp["eta"][0], rp["eta"][1]
        D2, SD, Dp = rp["D2"], red.SD, red.Dpoly
        D1 = red.D1
        mp = dict(Dp=Dp, dDp=D1(Dp))
        mp["ddDp"] = D1(mp["dDp"])
        mp["SDr"], _ = _dpow(SD, Dp, 1)
        mp["D2r"], _ = _dpow(D2, Dp, 3)
        mp["Bcr"], _ = _dpow(Bc, Dp, 2)
        mp["Acr"], _ = _dpow(Ac, Dp, 2)
        mp["Aer"], _ = _dpow(Ae, Dp, 1)
        mp["Ber"], _ = _dpow(Be, Dp, 1)
        mp["SD"], mp["dSD"] = SD, D1(SD)
        for n in ("SDr", "Acr", "Bcr", "Aer", "Ber", "D2r"):
            mp["d" + n] = D1(mp[n])
        for n in ("SDr", "Acr", "Bcr", "Aer", "Ber", "D2r"):
            mp["dd" + n] = D1(mp["d" + n])
        for n in ("Dp", "SDr", "Acr", "Bcr", "Aer", "Ber", "D2r"):
            mp["ddd" + n] = D1(mp["dd" + n])
        self.red, self.mp = red, mp
        if full:
            dAc, dBc, dBe, dD2, dSD = D1(Ac), D1(Bc), D1(Be), D1(D2), D1(SD)
            ddBc, ddD2 = D1(dBc), D1(dD2)
            TN1 = SD * Bc * (Ac + Be) + dBc * D2 - Bc * dD2
            TD1 = dD2 * Bc + D2 * dBc
            UN3c = dAc * D2 + SD * (Bc * Ae - Ac * Be)
            U1 = Bc * UN3c - Ac * dBc * D2
            TN1d = SD * SD * (dBc * (Ac + Be) + Bc * (dAc + dBe)) \
                + SD * (ddBc * D2 - Bc * ddD2) \
                - dSD * (dBc * D2 - Bc * dD2)
            QNm = U1 * SD * SD * Bc + TN1 * TN1 * SD * fmpq(1, 4) \
                - (TN1d * D2 * Bc - SD * TN1 * TD1) * fmpq(1, 2)
            mp["QNr"], _ = _dpow(QNm, Dp, 10)
            mp["TN1r"], _ = _dpow(TN1, Dp, 5)
            mp["dTN1r"] = D1(mp["TN1r"])
            mp["ddTN1r"] = D1(mp["dTN1r"])
            mp["dddTN1r"] = D1(mp["ddTN1r"])
            # exact reduced form of U1 = Dp^6 U6:  U6 = Bcr UN3r - Acr EBc D2r,
            # UN3r = EAc D2r + SDr (Bcr Aer - Acr Ber), EAc/EBc = (D1 of Dp^2 f)/Dp
            EAc = 2 * mp["dDp"] * mp["Acr"] + Dp * mp["dAcr"]
            EBc = 2 * mp["dDp"] * mp["Bcr"] + Dp * mp["dBcr"]
            UN3r = EAc * mp["D2r"] + mp["SDr"] * (mp["Bcr"] * mp["Aer"] - mp["Acr"] * mp["Ber"])
            U6 = mp["Bcr"] * UN3r - mp["Acr"] * EBc * mp["D2r"]
            if not (U1 - Dp ** 6 * U6).is_zero():
                raise RuntimeError("U6 reduced identity failed")
            for i, nm in ((2, "N"), (3, "W"), (4, "V")):     # d/du_i for mean value / E_D
                mp["gQNr_" + nm] = mp["QNr"].derivative(i)
                mp["gDp_" + nm] = Dp.derivative(i)
            mp["gQNr_A"], mp["gDp_A"] = mp["QNr"].derivative(1), Dp.derivative(1)
        self.cmp = {n: compile_mp(f) for n, f in mp.items()}
        self.dmax = [max(c[1][i] for c in self.cmp.values()) for i in range(4)]


def compile_mp(f):
    """[(kdeg, fmpq coef, eA, eN, eW, eV)] and the per-variable max degrees."""
    terms, dmax = [], [0, 0, 0, 0]
    for exps, coef in f.terms():
        e = list(map(int, exps))
        terms.append((e[9], coef, e[1], e[2], e[3], e[4]))
        for i in range(4):
            dmax[i] = max(dmax[i], e[1 + i])
    return terms, dmax


def powlist(x, n):
    out = [arb(1) if isinstance(x, arb) else acb(1)]
    for _ in range(n):
        out.append(out[-1] * x)
    return out


def kval(cmpf, pw, kd=None):
    """kappa-poly value (ascending coefficient list) of a compiled mpoly at power tables."""
    terms, dmax = cmpf
    kd = max(t[0] for t in terms) if kd is None else kd
    out = [arb(0) for _ in range(kd + 1)]
    for j, coef, a, b, c, e in terms:
        out[j] += arb(coef) * pw[0][a] * pw[1][b] * pw[2][c] * pw[3][e]
    return out


def u_to_pw(u, dmax):
    return [powlist(u[i], dmax[i]) for i in range(4)]


# ------------------------------------------------------------------ kappa-poly helpers
def pmul(a, b):
    out = [arb(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def padd(*ps):
    out = [arb(0)] * max(len(p) for p in ps)
    for p in ps:
        for i, x in enumerate(p):
            out[i] += x
    return out


def pscale(p, c):
    return [x * c for x in p]


def pder(p):
    return [p[i] * i for i in range(1, len(p))]


def phorner(p, k):
    acc = acb(0)
    for c in reversed(p):
        acc = acc * k + c
    return acc


def kdiv(QN, QD):
    """QN = (q2 k^2 + q1 k + q0) QD + rem (ball long division); QN deg<=12, QD deg 10."""
    R = list(QN) + [arb(0)] * (KD + 1 - len(QN))
    n = len(QD) - 1
    lead = QD[n]
    if not (lead != 0):
        raise ArithmeticError("QD leading coefficient contains 0")
    q = [arb(0), arb(0), arb(0)]
    for j in (KD, KD - 1, KD - 2):
        c = R[j] / lead
        q[j - n] = c
        for i in range(n + 1):
            R[i + j - n] -= c * QD[i]
        R[j] = arb(0)
    return q[2], q[1], q[0], R[:n]


# ------------------------------------------------------------------ half-plane bounds
def line_sup(rem, QD, c0, dens, nbox=48):
    """B0 >= sup_{Re kappa >= c0} |R|, B1 >= sup |dR/dkappa|, R = rem/QD (arb uppers).
    Maximum modulus: R is analytic on Re kappa >= c0 - DELC (pole condition (d)) and
    vanishes at infinity, so sup_{Re >= s} |R| is attained on Re kappa = s; the line
    is covered by tau in [0, TTAIL] boxes (conjugation symmetry) plus a coefficient
    tail; |QD| is lower-bounded factor-wise (``dens``), and the numerator uses the
    midpoint value + width * |derivative| + width^2/2 * |curvature| (recentring; the
    coefficient balls are tight, the box width is the enemy).  B1 is the smaller of
    the direct bound (numerator N2 = rem' QD - rem QD', |QD|^2 below) and the Cauchy
    shift sup_{Re >= c0 - DELC} |R| / DELC."""

    def sup_line(s0, num, dpow):
        nump = pder(num)
        numpp = pder(nump)

        def box(a, b, depth=0):
            kap = acb(arb(s0), arb((a + b) / 2, (b - a) / 2))
            lo = arb(1)
            for d in dens:
                fl = arb(abs(phorner(d, kap)).lower())
                if not (fl > 0):
                    lo = None
                    break
                lo = lo * fl
            if lo is not None:
                kap0 = acb(arb(s0), arb((a + b) / 2))
                wd = arb((b - a) / 2)
                up = arb(abs(phorner(num, kap0)).upper()) \
                    + wd * arb(abs(phorner(nump, kap0)).upper()) \
                    + wd * wd / 2 * arb(abs(phorner(numpp, kap)).upper())
                return arb((up / lo ** dpow).upper())
            if depth >= 9:
                raise ArithmeticError("denominator lower bound stays 0 on a sup line")
            mid = (a + b) / 2
            return box(a, mid, depth + 1).max(box(mid, b, depth + 1))

        edges = [12.0 * i / nbox for i in range(nbox)] \
            + [12.0 * (TTAIL / 12.0) ** (i / 16) for i in range(17)]
        B = arb(0)
        for a, b in zip(edges, edges[1:]):
            B = B.max(box(a, b))
        # coefficient tail, valid for every kappa with |kappa| >= TTAIL
        T = arb(TTAIL)
        cQD, nQD = arb(1), 0
        for d in dens:
            nd = len(d) - 1
            cd = abs(d[nd])
            for jj in range(nd):
                cd -= abs(d[jj]) / T ** (nd - jj)
            if not (cd.lower() > 0):
                raise ArithmeticError("den-factor tail coefficient bound not positive")
            cQD, nQD = cQD * cd, nQD + nd
        gap = dpow * nQD - len(num) + 1          # decay order of num/QD^dpow at infinity
        Cnum = sum((abs(c) / T ** (len(num) - 1 - jj) for jj, c in enumerate(num)), arb(0))
        B = B.max((Cnum / (cQD ** dpow * T ** gap)).upper())
        if not B.is_finite():
            raise ArithmeticError("line-sup bound not finite")
        return arb(B.upper())

    B0 = sup_line(c0, rem, 1)
    N2 = padd(pmul(pder(rem), QD), pscale(pmul(rem, pder(QD)), arb(-1)))
    B1 = sup_line(c0, N2, 2).min(arb((sup_line(c0 - DELC, rem, 1) / DELC).upper()))
    return B0, arb(B1.upper())


def hurwitz_shift(bc, rbar):
    """True iff every root of the cubic bc (ball coeffs, ascending) has Re < rbar."""
    if len(bc) < 4 or not (bc[3] != 0):
        return False
    r = arb(rbar)
    a0 = bc[0] + bc[1] * r + bc[2] * r ** 2 + bc[3] * r ** 3
    a1 = bc[1] + 2 * bc[2] * r + 3 * bc[3] * r ** 2
    a2 = bc[2] + 3 * bc[3] * r
    a3 = bc[3]
    if bool(a3 < 0):
        a0, a1, a2, a3 = -a0, -a1, -a2, -a3
    elif not bool(a3 > 0):
        return False
    return bool(a2 > 0) and bool(a0 > 0) and bool(a1 > 0) and bool(a2 * a1 > a3 * a0)


# ------------------------------------------------------------------ the box conditions
def scaled_to_u(Z):
    """(n, w, v, e^x) balls -> (A, N, W, V) via the background constraint."""
    t = Z[3]
    N, W, V = Z[0] / t, Z[1] * t * t, Z[2] * t
    A = 1 + 2 * W * (1 + V * V / 3 + 4 * N * V / 3) / (1 - V * V)
    return [A, N, W, V]


class J2:
    """Order-2 jet of a kappa-poly-valued function of x on a box [m - w, m + w]:
    (v, d) = value and d/dx at the midpoint (tight), (vb, db, cb) = enclosures of the
    value, d/dx and d^2/dx^2 over the box.  enclose() gives the Taylor--Lagrange ball
    f(m) + (x - m) f'(m) + (x - m)^2/2 f''(xi)."""

    __slots__ = ("v", "d", "vb", "db", "cb")

    def __init__(self, v, d, vb, db, cb):
        self.v, self.d, self.vb, self.db, self.cb = v, d, vb, db, cb


def j2add(*xs):
    return J2(*(padd(*(getattr(x, sl) for x in xs)) for sl in J2.__slots__))


def j2scale(x, c):
    c = arb(c) if not isinstance(c, arb) else c
    return J2(*(pscale(getattr(x, sl), c) for sl in J2.__slots__))


def j2mul(x, y):
    return J2(pmul(x.v, y.v), padd(pmul(x.v, y.d), pmul(x.d, y.v)),
              pmul(x.vb, y.vb), padd(pmul(x.vb, y.db), pmul(x.db, y.vb)),
              padd(pmul(x.cb, y.vb), pscale(pmul(x.db, y.db), arb(2)), pmul(x.vb, y.cb)))


def j2inv0(x):
    """Reciprocal of a kappa-degree-0 jet (scalar); needs the value balls nonzero."""
    v, vb = x.v[0], x.vb[0]
    if not (abs(v).lower() > 0) or not (abs(vb).lower() > 0):
        raise ArithmeticError("jet reciprocal of a ball containing 0")
    return J2([1 / v], [-x.d[0] / v ** 2], [1 / vb], [-x.db[0] / vb ** 2],
              [(2 * x.db[0] ** 2 / vb - x.cb[0]) / vb ** 2])


def jcoef(x, i):
    z = arb(0)
    return J2(*([getattr(x, sl)[i] if i < len(getattr(x, sl)) else z]
                for sl in J2.__slots__))


def j2horner(x, z):
    out = []
    for sl in J2.__slots__:
        acc = acb(0)
        for c in reversed(getattr(x, sl)):
            acc = acc * z + c
        out.append([acc])
    return J2(*out)


def j2kder(x):
    return J2(*(pder(getattr(x, sl)) for sl in J2.__slots__))


def kdiv_jet(QNj, QDj):
    """Quotient scalar-jets (q0, q1, q2) of the kappa-division QN/QD (deg 12 by 10),
    slot-synchronized so that enclosures of the quotient stay correlated."""
    leadinv = j2inv0(jcoef(QDj, 10))
    R = [jcoef(QNj, i) for i in range(KD + 1)]
    q = [None, None, None]
    for j in (KD, KD - 1, KD - 2):
        c = j2mul(R[j], leadinv)
        q[j - 10] = c
        for i in range(11):
            R[i + j - 10] = j2add(R[i + j - 10], j2scale(j2mul(c, jcoef(QDj, i)), -1))
    return q, R[:10]


class BoxEval:
    """Atom jets on an x-box anchored at its midpoint (all through the exact reduced
    mpolys; derivatives are D1-atom values divided by the S*Delta value)."""

    def __init__(self, at, u_m, u_box, w):
        self.at, self.dx = at, arb(0, w)
        self.dx2h = self.dx * self.dx / 2
        dm = at.dmax
        self.pm, self.pb = u_to_pw(u_m, dm), u_to_pw(u_box, dm)
        self.sd_m = kval(at.cmp["SD"], self.pm, 0)[0]
        self.sd_b = kval(at.cmp["SD"], self.pb, 0)[0]
        self.dsd_b = kval(at.cmp["dSD"], self.pb, 0)[0]
        if not (self.sd_b != 0) or not (self.sd_m != 0):
            raise ArithmeticError("S*Delta contains 0 on a box")

    def jet(self, n):
        a = kval(self.at.cmp[n], self.pm)
        kd = len(a) - 1
        d = pscale(kval(self.at.cmp["d" + n], self.pm, kd), 1 / self.sd_m)
        vb = kval(self.at.cmp[n], self.pb, kd)
        dbraw = kval(self.at.cmp["d" + n], self.pb, kd)
        db = pscale(dbraw, 1 / self.sd_b)
        ddraw = kval(self.at.cmp["dd" + n], self.pb, kd)
        cb = [y / self.sd_b ** 2 - x * self.dsd_b / self.sd_b ** 3
              for x, y in zip(dbraw, ddraw)]
        return J2(a, d, vb, db, cb)

    def enclose(self, x):
        return [v + self.dx * d + self.dx2h * c for v, d, c in
                zip(x.v, padd(x.d, [arb(0)] * (len(x.v) - len(x.d))),
                    padd(x.cb, [arb(0)] * (len(x.v) - len(x.cb))))]

    def enclose1(self, x):
        """First-order (mean-value) enclosure v(m) + dx * hull(d/dx over the box);
        better than enclose() for division outputs whose curvature slot is inflated."""
        return [v + self.dx * d for v, d in
                zip(x.v, padd(x.db, [arb(0)] * (len(x.v) - len(x.db))))]


def kval_s(cmpf, pw, pwt, kd=None):
    """kappa-poly value of a compiled mpoly under N = n/t, W = w t^2, V = v t: the
    t-powers are explicit per term (pwt maps integer p -> t^p ball), so the deep-x
    t-power cancellations of the unscaled variables never happen."""
    terms, dmax = cmpf
    kd = max(t[0] for t in terms) if kd is None else kd
    out = [arb(0) for _ in range(kd + 1)]
    for j, coef, a, b, c, e in terms:
        out[j] += arb(coef) * pw[0][a] * pw[1][b] * pw[2][c] * pw[3][e] * pwt[-b + 2 * c + e]
    return out


class BoxEvalScaled:
    """BoxEval in the tube's scaled variables (n, w, v, t = e^x): same jet interface,
    with every atom evaluated through kval_s.  Z balls are (n, w, v, t)."""

    def __init__(self, at, Zm, Zb, w):
        self.at, self.dx = at, arb(0, w)
        self.dx2h = self.dx * self.dx / 2
        dm = at.dmax

        def tables(Z):
            t = Z[3]
            S = 1 - (Z[2] * t) ** 2
            A = 1 + 2 * (Z[1] * t * t) * (1 + (Z[2] * t) ** 2 / 3 + 4 * Z[0] * Z[2] / 3) / S
            pw = [powlist(A, dm[0]), powlist(Z[0], dm[1]), powlist(Z[1], dm[2]),
                  powlist(Z[2], dm[3])]
            pmax = dm[1] + 2 * dm[2] + dm[3] + 2
            tp, ti = powlist(t, pmax), powlist(1 / t, pmax)
            pwt = {i: tp[i] for i in range(pmax + 1)}
            pwt.update({-i: ti[i] for i in range(1, pmax + 1)})
            return pw, pwt, A

        self.pm_, self.pmt, self.Am = tables(Zm)
        self.pb_, self.pbt, self.Ab = tables(Zb)
        self.sd_m = kval_s(at.cmp["SD"], self.pm_, self.pmt, 0)[0]
        self.sd_b = kval_s(at.cmp["SD"], self.pb_, self.pbt, 0)[0]
        self.dsd_b = kval_s(at.cmp["dSD"], self.pb_, self.pbt, 0)[0]
        if not (self.sd_b != 0) or not (self.sd_m != 0):
            raise ArithmeticError("S*Delta contains 0 on a box")

    def jet(self, n):
        a = kval_s(self.at.cmp[n], self.pm_, self.pmt)
        kd = len(a) - 1
        d = pscale(kval_s(self.at.cmp["d" + n], self.pm_, self.pmt, kd), 1 / self.sd_m)
        vb = kval_s(self.at.cmp[n], self.pb_, self.pbt, kd)
        dbraw = kval_s(self.at.cmp["d" + n], self.pb_, self.pbt, kd)
        db = pscale(dbraw, 1 / self.sd_b)
        ddraw = kval_s(self.at.cmp["dd" + n], self.pb_, self.pbt, kd)
        cb = [y / self.sd_b ** 2 - x * self.dsd_b / self.sd_b ** 3
              for x, y in zip(dbraw, ddraw)]
        return J2(a, d, vb, db, cb)

    enclose = BoxEval.enclose
    enclose1 = BoxEval.enclose1


def qn_qd(at, ev):
    """QN, QD kappa-polys at a BoxEval box via order-2 jets of the reduced atoms:
        T = TN1r/(P W1), W1 = SDr D2r Bcr, U = U6/(P^3 SDr D2r^2 Bcr),
        U6 = Bcr UN3r - Acr EBc D2r (exact),  QN_red-form of Q = U + T^2/4 - T'/2:
        QN = U6 SDr^2 Bcr + P SDr TN1r^2/4 - (P/2)[dTN1r W1/sd... ] (see the note),
        QD = P^3 SDr^3 D2r^2 Bcr^2."""
    Dp, dDp = ev.jet("Dp"), ev.jet("dDp")
    SDr, dSDr = ev.jet("SDr"), ev.jet("dSDr")
    Acr, dAcr = ev.jet("Acr"), ev.jet("dAcr")
    Bcr, dBcr = ev.jet("Bcr"), ev.jet("dBcr")
    Aer, Ber, dBer = ev.jet("Aer"), ev.jet("Ber"), ev.jet("dBer")
    D2r, dD2r = ev.jet("D2r"), ev.jet("dD2r")
    TN1r, dTN1r = ev.jet("TN1r"), ev.jet("dTN1r")
    sdinv = j2inv0(j2mul(Dp, SDr))
    Pp = j2mul(dDp, sdinv)
    EAc = j2add(j2scale(j2mul(dDp, Acr), 2), j2mul(Dp, dAcr))
    EBc = j2add(j2scale(j2mul(dDp, Bcr), 2), j2mul(Dp, dBcr))
    UN3r = j2add(j2mul(EAc, D2r),
                 j2mul(SDr, j2add(j2mul(Bcr, Aer), j2scale(j2mul(Acr, Ber), -1))))
    U6 = j2add(j2mul(Bcr, UN3r), j2scale(j2mul(j2mul(Acr, EBc), D2r), -1))
    DB = j2mul(D2r, Bcr)
    W1 = j2mul(SDr, DB)
    W1p = j2mul(sdinv, j2add(j2mul(dSDr, DB),
                             j2mul(SDr, j2add(j2mul(dD2r, Bcr), j2mul(D2r, dBcr)))))
    br = j2add(j2mul(dTN1r, W1),
               j2scale(j2mul(SDr, j2mul(TN1r, j2add(j2mul(Pp, W1), j2mul(Dp, W1p)))), -1))
    QNj = j2add(j2mul(U6, j2mul(j2mul(SDr, SDr), Bcr)),
                j2scale(j2mul(j2mul(TN1r, TN1r), j2mul(Dp, SDr)), arb(fmpq(1, 4))),
                j2scale(j2mul(br, Dp), arb(fmpq(-1, 2))))
    S3 = j2mul(SDr, j2mul(SDr, SDr))
    P3 = j2mul(Dp, j2mul(Dp, Dp))
    QDrj = j2mul(S3, j2mul(j2mul(D2r, D2r), j2mul(Bcr, Bcr)))
    return QNj, QDrj, dict(Bcr=ev.enclose(Bcr), D2r=ev.enclose(D2r), P3=P3,
                           fac=(j2mul(P3, S3), D2r, Bcr))


def conditions_jet(QNj, QDj, u, bc_poly, c0, ev, fac=None, p3j=None):
    """Margins (arb) of (a)-(d) for one box; every margin must be > 0.  q_j and the
    remainder are computed in jet arithmetic (division and subtraction stay
    correlated; the point-anchored slot carries the true value).  Pole condition (d):
    Routh-Hurwitz on the shifted Bc-cubic plus the A(x), -F_N(x) balls."""
    qj, remj = kdiv_jet(QNj, QDj)          # QDj here is the P-free QDr
    p3 = ev.enclose(fac[0] if p3j is None else p3j)[0]
    if not (p3 != 0):
        raise ArithmeticError("P^3 enclosure contains 0")
    q0, q1, q2 = (ev.enclose(x)[0] / p3 for x in qj)
    rem = [ev.enclose(r)[0] for r in remj]
    A, W = u[0], u[2]
    FN = A - 2 - 2 * W / 3
    if not hurwitz_shift(bc_poly, RBAR):
        raise ArithmeticError("Bc-cubic roots not certified left of RBAR (Hurwitz)")
    dens = [[ev.enclose(jcoef(fac[0], 0))[0]], ev.enclose(fac[1]), ev.enclose(fac[1]),
            ev.enclose(fac[2]), ev.enclose(fac[2])]
    B0, B1 = line_sup(rem, pscale(ev.enclose(QDj), p3), c0, dens=dens)
    out = dict(q2=q2, mB=2 * c0 * q2 + q1 - B1, mC=q2 * c0 ** 2 + q1 * c0 + q0 - B0,
               poleA=arb(RBAR) - A, poleFN=arb(RBAR) + FN)
    ok = all(bool(m > 0) for m in out.values())
    return ok, out, (B0, B1)


def conditions(QN, QD, u, bc_poly, c0, dens=None):
    """Plain-poly version of the margins (used by the sonic closing box)."""
    q2, q1, q0, rem = kdiv(QN, QD)
    A, W = u[0], u[2]
    FN = A - 2 - 2 * W / 3
    if not hurwitz_shift(bc_poly, RBAR):
        raise ArithmeticError("Bc-cubic roots not certified left of RBAR (Hurwitz)")
    B0, B1 = line_sup(rem, QD, c0, dens=dens or [QD])
    out = dict(q2=q2, mB=2 * c0 * q2 + q1 - B1, mC=q2 * c0 ** 2 + q1 * c0 + q0 - B0,
               poleA=arb(RBAR) - A, poleFN=arb(RBAR) + FN)
    ok = all(bool(m > 0) for m in out.values())
    return ok, out, (B0, B1)


def kdiv_gen(num, den, qdeg):
    """num = q * den + rem with deg q = qdeg (ball long division); returns (q, rem)."""
    n = len(den) - 1
    lead = den[n]
    if not (lead != 0):
        raise ArithmeticError("leading denominator coefficient contains 0")
    R = list(num) + [arb(0)] * (n + qdeg + 1 - len(num))
    q = [arb(0)] * (qdeg + 1)
    for j in range(n + qdeg, n - 1, -1):
        c = R[j] / lead
        q[j - n] = c
        for i in range(n + 1):
            R[i + j - n] -= c * den[i]
        R[j] = arb(0)
    return q, R[:n]


# ------------------------------------------------------------------ region: tube steps
def run_tube(at, tube, c0, x_min=-4.5, wfrac=1500.0, prec=256):
    """Conditions (a)-(d) on every tube step with x - h >= x_min (the deep tail
    x < x_min is covered by the analytic Lemma D of notes/s4-energy.md; the interval
    margins genuinely die around x ~ -5.5, where the t^2-relational cancellation of
    the kappa-rational remainder crosses the enclosure floor).  Each step is split
    into sub-boxes of half-width ~ |x|/wfrac (adaptive doubling on failure)."""
    from . import shootsys as ss
    mins, fails, t0, nboxes = {}, [], time.time(), 0
    with precision(prec):
        for sd in tube.steps:
            if float(sd.x) - sd.h < x_min:
                continue
            h = arb(sd.h)
            nsub = max(1, min(512, int(sd.h * wfrac / (2 * abs(float(sd.x))))))
            ok = False
            while not ok and nsub <= 512:
                ok, worst = True, {}
                for j in range(nsub):
                    w = sd.h / (2 * nsub)
                    sm = (-h * (2 * j + 1)) / (2 * nsub)
                    Zm = [zi + arb(0, sd.eps_z) for zi in ss.horner_vec(sd.co, sm)]
                    Zb = [zi + arb(0, sd.eps_z)
                          for zi in ss.horner_vec(sd.co, sm + arb(0, w))]
                    try:
                        ev = BoxEvalScaled(at, Zm, Zb, w)
                        QNj, QDj, v = qn_qd(at, ev)
                        u = [ev.Ab, Zb[0] / Zb[3], Zb[1] * Zb[3] ** 2, Zb[2] * Zb[3]]
                        okb, m, _ = conditions_jet(QNj, QDj, u, v["Bcr"], c0, ev,
                                                   fac=v["fac"], p3j=v["P3"])
                    except (ArithmeticError, RuntimeError):
                        okb, m = False, {}
                    nboxes += 1
                    ok = ok and okb
                    if not okb:
                        break
                    for k, val in m.items():
                        worst[k] = val if k not in worst else worst[k].min(val)
                if not ok:
                    nsub *= 2
            for k, val in worst.items():
                mins[k] = val if k not in mins else mins[k].min(val)
            if not ok:
                fails.append(float(sd.x))
    return dict(ok=not fails, fails=fails, x_min=x_min, boxes=nboxes,
                time=time.time() - t0,
                mins={k: float(v.lower()) for k, v in mins.items()})


# ------------------------------------------------------------------ region: sonic end
V0C, V0W = "0.1124394013880983", 1e-16
SIGMA0 = (0.090105970507920, 5e-16)          # certified ball (linsonic4.sigma, A1 bg)


def sonic_expansion_cached(prec=256):
    from . import sonic
    with precision(prec):
        ex = sonic.sonic_expansion(V0C, K=41, width=V0W, m=5)
        cert = ex.certify()
        if not cert.ok:
            raise RuntimeError("sonic background certificate failed")
    return ex


def u_prime_balls(ex, xball, prec=256):
    """(A', N', W', V') enclosures over the x-ball: series derivative + Cauchy tail."""
    with precision(prec):
        rho = arb(ex.cert.nu) * arb(0.75)
        xmax = abs(xball).upper()
        if not (rho > arb(xmax)):
            raise RuntimeError("Cauchy radius does not cover the x-ball")
        tb = ex.cert.tail_bound(rho) / (rho - arb(xmax))
        out = []
        for i in range(4):
            acc = arb(0)
            for n in reversed(range(1, ex.K + 1)):
                acc = acc * xball + n * ex.balls[n][i]
            out.append(acc + arb(0, tb.upper()))
    return out


def _sonic_box(at, ex, a, b, c0, mins, stats, depth=0):
    m, w = (a + b) / 2, (b - a) / 2
    okb = None
    if w <= abs(m) / 1500:              # heuristic: don't even try wide boxes
        try:
            ev = BoxEval(at, ex.eval(arb(m)), ex.eval(arb(m, w)), w)
            QNj, QDj, v = qn_qd(at, ev)
            okb, mg, _ = conditions_jet(QNj, QDj, ex.eval(arb(m, w)), v["Bcr"], c0, ev,
                                        fac=v["fac"], p3j=v["P3"])
        except (ArithmeticError, RuntimeError):
            okb, mg = False, None
    if okb:
        stats[0] += 1
        for k, val in mg.items():
            mins[k] = val if k not in mins else mins[k].min(val)
        return True
    if depth >= 52:
        stats[1] += 1
        return False
    mid = (a + b) / 2
    okl = _sonic_box(at, ex, a, mid, c0, mins, stats, depth + 1)
    return _sonic_box(at, ex, mid, b, c0, mins, stats, depth + 1) and okl


def closing_box(at, ex, c0, dc=1e-8, prec=256):
    """The closing box [-dc, 0]: x^2 Q = QQ/(E_D^3 QD_red) with QN_red = x QQ (mean
    value; the exact zero at x = 0 is the field identity QN_red(u0) = 0 in
    Q(V0, sqrt3)).  All factors are nonvanishing balls, so the conditions hold up to
    and including the endpoint.  Also returns the sigma known-answer containments and
    the xT = 1 + sigma(kappa) exponent data for the boundary lemma."""
    with precision(prec):
        xb = arb(-dc / 2, dc / 2)
        u = ex.eval(xb)
        up = u_prime_balls(ex, xb, prec)
        pw = u_to_pw(u, at.dmax)
        QQ = [arb(0)] * (KD + 1)
        for nm, d in (("A", up[0]), ("N", up[1]), ("W", up[2]), ("V", up[3])):
            g = kval(at.cmp["gQNr_" + nm], pw, KD)
            for j in range(KD + 1):
                QQ[j] += g[j] * d
        ED = sum((kval(at.cmp["gDp_" + nm], pw, 0)[0] * d for nm, d in
                  (("A", up[0]), ("N", up[1]), ("W", up[2]), ("V", up[3]))), arb(0))
        if not (ED > 0):
            raise ArithmeticError("E_D = D(x)/x not certified positive")
        SDr = kval(at.cmp["SDr"], pw, 0)[0]
        D2r, Bcr = kval(at.cmp["D2r"], pw), kval(at.cmp["Bcr"], pw)
        QDr = pmul(pmul([SDr ** 3], pmul(D2r, D2r)), pmul(Bcr, Bcr))
        den = pscale(QDr, ED ** 3)
        okc, m, _ = conditions(QQ, den, u, Bcr, c0,
                               dens=[[ED ** 3 * SDr ** 3], D2r, D2r, Bcr, Bcr])
        # known answer: x^2 Q -> (sigma(kappa)^2 - 1)/4 at x = 0
        s0 = arb(SIGMA0[0], SIGMA0[1])
        L = [(s0 * s0 - 1) / 4, s0 * (1 - s0) / 2, (1 - s0) ** 2 / 4]
        q2, q1, q0, _rem = kdiv(QQ, den)
        ka = [bool(q.contains(Lj)) for q, Lj in zip((q0, q1, q2), L)]
        # xT enclosure: exponent data 1 + sigma(kappa) for the boundary lemma
        TN1r = kval(at.cmp["TN1r"], pw)
        denT = pscale(pmul(D2r, Bcr), ED * SDr)
        qT, remT = kdiv_gen(TN1r, denT, len(TN1r) - len(denT))
        BT0, _ = line_sup(remT, denT, c0, dens=[[ED * SDr], D2r, Bcr])
        return okc, m, dict(known_answer_sigma=ka,
                            ED=[float(ED.lower()), float(ED.upper())],
                            xT_affine=[float(qT[0].mid()),
                                       float(qT[1].mid()) if len(qT) > 1 else 0.0],
                            xT_rem_bound=float(BT0.upper()), dc=dc)


def run_sonic(at, c0, dc=1e-8, prec=256):
    """Ladder [-0.05, -dc]: adaptive bisection with the jet pipeline; then the closing
    box [-dc, 0] (``closing_box``)."""
    ex, t0 = sonic_expansion_cached(prec), time.time()
    with precision(prec):
        mins, stats = {}, [0, 0]
        ok_all = _sonic_box(at, ex, -0.05, -dc, c0, mins, stats)
        okc, m, extra = closing_box(at, ex, c0, dc, prec)
        ok_all = ok_all and okc
        for k, val in m.items():
            mins[k] = val if k not in mins else mins[k].min(val)
        out = dict(ok=ok_all, mins={k: float(v.lower()) for k, v in mins.items()},
                   boxes=stats[0], box_fails=stats[1], time=time.time() - t0, **extra)
    return out


# --------------------------------------------------- the exact field identity (sonic)
def field_identity(at=None):
    """QN_red(u0) = 0 exactly in Q(V0, sqrt3) at the KHA closed-form sonic point,
    identically in V0 and kappa (the order-11 vanishing of QN = Dpoly^10 QN_red).
    ~2 min of exact fmpq_poly arithmetic; returns True iff the identity holds."""
    from flint import fmpq_poly as P
    at = at or Atoms(full=True)
    fm = lambda x, y: (x[0] * y[0] + 3 * (x[1] * y[1]), x[0] * y[1] + x[1] * y[0])
    ONE, Sp = (P([1]), P([0])), (P([0]), P([1]))
    msV = (P([1]), P([0, -1]))                                   # 1 - sqrt3 V
    qF = fm(fm((P([0]), P([4])), (P([1, 0, -1]), P([0]))), msV)  # 4 s (1-V^2)(1-sV)
    pA = fm(fm(Sp, msV), (P([7, 0, -3]), P([0, 2])))
    pN = fm((P([0]), P([4])), fm((P([1, 0, -1]), P([0])), (P([0, -1]), P([1]))))
    pW = fm((P([fmpq(3, 2)]), P([0])), fm(msV, (P([0, -2]), P([1, 0, -1]))))
    pV = fm((P([0, 1]), P([0])), qF)
    terms, dmax = at.cmp["QNr"]
    DM = sum(dmax)
    pows = lambda p, n: [ONE] + [None] * n
    PA, PN, PW, PV, PQ = (pows(p, d) for p, d in
                          ((pA, dmax[0]), (pN, dmax[1]), (pW, dmax[2]), (pV, dmax[3]), (qF, DM)))
    for tab, base in ((PA, pA), (PN, pN), (PW, pW), (PV, pV), (PQ, qF)):
        for i in range(1, len(tab)):
            tab[i] = fm(tab[i - 1], base)
    acc = {}
    for j, coef, a, b, c, e in terms:
        t = fm(fm(fm(PA[a], PN[b]), fm(PW[c], PV[e])), PQ[DM - a - b - c - e])
        cp = P([coef])
        cur = acc.get(j, (P([0]), P([0])))
        acc[j] = (cur[0] + cp * t[0], cur[1] + cp * t[1])
    return all(v[0].degree() < 0 and v[1].degree() < 0 for v in acc.values())


# ------------------------------------------------------------------------------ main
def main(argv):
    mode = argv[0] if argv else "all"
    c0 = float(argv[2]) if len(argv) > 2 else 6.0
    if c0 <= RBAR:
        raise SystemExit(f"c0 = {c0} must exceed the pole bound RBAR = {RBAR}")
    t0 = time.time()
    at = Atoms(full=True)
    res = dict(c0=c0, rbar=RBAR, ttail=TTAIL, atoms_time=time.time() - t0)
    if mode == "field":
        res["field_identity"] = field_identity(at)
        print("QN_red(u0) == 0 exactly in Q(V0, sqrt3):", res["field_identity"])
        return res
    if mode in ("sonic", "all"):
        res["sonic"] = run_sonic(at, c0)
        print(f"sonic:  ok={res['sonic']['ok']}  mins={res['sonic']['mins']}  "
              f"sigma-known-answer={res['sonic']['known_answer_sigma']}")
    if mode == "centre":
        raise SystemExit("the deep tail x < -4.5 is Lemma D (notes/s4-energy.md "
                         "section 4); its closing computation is an S4-3 item")
    if mode in ("tube", "all"):
        from .lintube import Tube
        path = argv[1] if len(argv) > 1 else os.environ.get("P4_TUBE_CACHE", "")
        if not path or not os.path.exists(path):
            raise SystemExit("tube mode needs the certified tube JSON (arg or P4_TUBE_CACHE)")
        tube = Tube.load(path)
        res["tube_info"] = {k: str(v) for k, v in tube.info.items()}
        res["tube"] = run_tube(at, tube, c0)
        print(f"tube:   ok={res['tube']['ok']}  steps={len(tube.steps)}  "
              f"mins={res['tube']['mins']}  fails={res['tube']['fails'][:5]}")
    ok = all(res[r]["ok"] for r in ("sonic", "centre", "tube") if r in res)
    res["all_ok"], res["total_time"] = ok, time.time() - t0
    out = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                       "results", "theorem_c", "tc_energy.json")
    if mode == "all":
        json.dump(res, open(out, "w"), indent=1)
        print(f"certificate {'PASSES' if ok else 'FAILS'} at c0 = {c0}; written {out}")
    return res


if __name__ == "__main__":
    main(sys.argv[1:])
