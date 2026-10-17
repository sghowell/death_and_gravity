"""Theorem C, S4-3: (1) Lemma D closure (deep tail of the direct-energy
certificate, making P4-8 unconditional) and (2) the commutator exclusion on the
strip  Re kappa in [0, 6], |Im kappa| >= b1  (route (b2), GHJS-style).

Both certificates run in Taylor-model (TM) arithmetic in the local variable of
each region (tc_tm): y = e^x on the centre region, sigma = x_k - x on tube
steps, xi = -x on the sonic region.  Nothing here modifies tc_energy's
certified [-4.5, 0] run; Lemma D covers x <= -4.5 in one TM computation.

(1) Lemma D.  QN and QD (tc_energy's kappa-polynomials with Q = QN/QD) are
assembled as kappa-polys of Laurent TM coefficients on y in (0, e^{-4.5}];
the *exact layer identity*  (verified at import, `check_layer_identity`):
      Z := QNr - (9/4) Dp^3 SDr^3 D2r^2 Bcr^2   has empty t-layers < 16
      in the centre grading  p = -eN + 2eW + eV  with A -> 1  (and no odd
      layers at all),
together with A - 1 = O(t^2), gives val(Z o u) >= 16 while val(QD o u) = 14,
so  Q - 9/4 = Z/QD  evaluates junk-free:  q2, q1, q0 - 9/4 and the remainder
carry the honest t^2 scale and conditions (a)-(d) of tc_energy hold on the
whole open tail with hat-margins ~ q2/t^2 = O(1).

(2) Strip exclusion (m = 1 commutation, c = -1/32).  For an eigenfunction chi
(class of s2-theorem-b.md 3.8) the scalar reduction gives chi'' + A0 chi' +
B0 chi = 0 with A0 = -T, B0 = -U exact rational in kappa (tc_reduce atoms).
Where B0 != 0, chi1 := chi' solves chi1'' + A1 chi1' + B1 chi1 = 0 with
    A1 = A0 - (log B0)',   B1 = B0 + A0' - A0 (log B0)'.
With tt := the kappa-lead of T (real), c < 0, weight g = exp int(Re A1 - c tt),
testing with g conj(chi1) and combining Re + Im parts (toy-verified identity):
    0 = -int g|chi1'|^2 + int g K |chi1|^2 - int g nu Im(conj(chi1) chi1'),
    K = P - (tau/c) Q - (1/2)(f' + rho f),   P + iQ = B1,   nu = Im A1 + tau tt,
    rho = Re A1 - c tt,   f = (c + tau^2/c) tt - (tau/c) nu,
boundary terms vanishing per the endpoint lemma (note s4-commutator.md).
Young: |nu Im(conj chi1 chi1')| <= |chi1'|^2 + nu^2 |chi1|^2 / 4, so
    (K + nu^2/4)(x) < 0 for all x in (-inf, 0)   =>   chi1 == 0  =>  chi == 0
(B0 != 0 unwinds), => p = const * g => const = 0: no eigenvalue.  The
certificate checks K + nu^2/4 < 0 as a tau-quadratic with explicit
remainder bounds, uniformly for eps in [0, 6], tau >= b1 (tau > 0 suffices:
real coefficients pair kappa with conj kappa).

Certificate JSON: results/theorem_c/tc_commutator.json.  Run:
  PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_commutator all \
      <tube.json> [b1]
Modes: lemma_d | centre | tube | sonic | all | checks.
"""
from __future__ import annotations

import json
import os
import sys
import time

from flint import arb, fmpq, fmpq_mpoly_ctx

from .arbseries import precision
from .tc_energy import Atoms, SIGMA0, V0C, V0W, hurwitz_shift
from .tc_tm import (AtomKP, LC, TM, centre_components, jlc_add, jlc_inv,
                    jlc_mul, jlc_neg, kpj_add, kpj_deriv, kpj_depth, kpj_kmul,
                    kpj_mul, kpj_neg, kpj_scale, lc_raise, lc_zero,
                    sonic_components, step_components, tm_zero)

MC = 1                      # commutation level m
CNUM, CDEN = -1, 32         # c = -1/32 (exact dyadic)
B1_DEFAULT = 14.0           # target |Im kappa| threshold
EPS_MAX = 6.0               # strip Re kappa in [0, EPS_MAX] (= P4-8's c0)
C0_LEMMA_D = 6.0            # Lemma D threshold (tc_energy's c0)
YCUT = 0.011108996538242306   # e^{-4.5}: centre region y in (0, YCUT]
XI_MAX = 0.05               # sonic region xi = -x in [0, XI_MAX]


def cval():
    return arb(CNUM) / CDEN


def check_layer_identity(at, pshift=2):
    """(E2) Z = QNr - (9/4) Dp^3 SDr^3 D2r^2 Bcr^2 has empty layers below
    pmin+2 in the centre grading with A -> 1, and no odd-parity terms.
    Returns (pmin, ok).  ~20 s of exact arithmetic."""
    QQD = at.mp["Dp"] ** 3 * at.mp["SDr"] ** 3 * at.mp["D2r"] ** 2 * at.mp["Bcr"] ** 2
    Z = at.mp["QNr"] - fmpq(9, 4) * QQD
    ctx3 = fmpq_mpoly_ctx.get(("n", "w", "v", "k"))
    pmin, odd = None, 0
    lay = {}
    for e, cf in Z.terms():
        eN, eW, eV, ek = int(e[2]), int(e[3]), int(e[4]), int(e[9])
        p = -eN + 2 * eW + eV
        pmin = p if pmin is None else min(pmin, p)
        if (eN + eV) % 2:
            odd += 1
    for e, cf in Z.terms():
        eN, eW, eV, ek = int(e[2]), int(e[3]), int(e[4]), int(e[9])
        p = -eN + 2 * eW + eV
        if p < pmin + pshift:
            lay.setdefault(p, ctx3.constant(0))
            lay[p] += ctx3.from_dict({(eN, eW, eV, ek): cf})
    ok = odd == 0 and all(l.is_zero() for l in lay.values())
    return pmin, ok


def check_u6_top_identity(at):
    """(E1) U6_9 * (-Dp) = (3 - V^2) * [kappa^7 coefficient of denU], and
    Dp^2 | U6_9, where U6_9 is the kappa^9 coefficient of U6 and denU =
    Dp^3 SDr D2r^2 Bcr.  This is the exact form of  w * det C = -(3 - V^2)
    (s4-reduction.md 1.3) and justifies the val-2 Laurent raise of the kdiv
    leads at the sonic point (Dp^2 factor, D(u0) = 0 exact)."""
    red = at.red
    ctx = red.ctx
    _t, A, N, W, V, *_ = ctx.gens()

    def kcoef(f, j):
        out = ctx.constant(0)
        for e, cf in f.terms():
            if int(e[9]) == j:
                ee = list(map(int, e))
                ee[9] = 0
                out += ctx.from_dict({tuple(ee): cf})
        return out

    mp = at.mp
    EAc = 2 * mp["dDp"] * mp["Acr"] + mp["Dp"] * mp["dAcr"]
    EBc = 2 * mp["dDp"] * mp["Bcr"] + mp["Dp"] * mp["dBcr"]
    UN3r = EAc * mp["D2r"] + mp["SDr"] * (mp["Bcr"] * mp["Aer"]
                                          - mp["Acr"] * mp["Ber"])
    U6 = mp["Bcr"] * UN3r - mp["Acr"] * EBc * mp["D2r"]
    U6_9 = kcoef(U6, 9)
    den_top = kcoef(mp["Dp"] ** 3 * mp["SDr"] * mp["D2r"] ** 2 * mp["Bcr"], 7)
    ok1 = (U6_9 * (-mp["Dp"]) - (3 - V ** 2) * den_top).is_zero()
    q, r = divmod(U6_9, mp["Dp"])
    ok2 = r.is_zero()
    if ok2:
        _, r2 = divmod(q, mp["Dp"])
        ok2 = r2.is_zero()
    return ok1 and ok2


def check_sonic_zeros(at):
    """(E3) D(u0) = 0, PW(u0) = 0, PV(u0) = 0 exactly in Q(V0, sqrt3) at the
    KHA closed-form sonic point: the cleared flow du and Delta~ vanish at x=0,
    justifying the val-1 Laurent raises of SD, Dp and the D1-atom values there.
    Uses the field arithmetic of tc_energy.field_identity (small polys)."""
    from flint import fmpq_poly as P
    red = at.red
    fm = lambda x, y: (x[0] * y[0] + 3 * (x[1] * y[1]), x[0] * y[1] + x[1] * y[0])
    ONE, Sp = (P([1]), P([0])), (P([0]), P([1]))
    msV = (P([1]), P([0, -1]))
    qF = fm(fm((P([0]), P([4])), (P([1, 0, -1]), P([0]))), msV)
    pA = fm(fm(Sp, msV), (P([7, 0, -3]), P([0, 2])))
    pN = fm((P([0]), P([4])), fm((P([1, 0, -1]), P([0])), (P([0, -1]), P([1]))))
    pW = fm((P([fmpq(3, 2)]), P([0])), fm(msV, (P([0, -2]), P([1, 0, -1]))))
    pV = fm((P([0, 1]), P([0])), qF)

    def field_eval(f):
        acc = (P([0]), P([0]))
        dmax = sum(f.degrees()[i] for i in (1, 2, 3, 4))
        for e, cf in f.terms():
            a, b, c, ee = int(e[1]), int(e[2]), int(e[3]), int(e[4])
            t = (P([1]), P([0]))
            for base, n in ((pA, a), (pN, b), (pW, c), (pV, ee)):
                for _ in range(n):
                    t = fm(t, base)
            for _ in range(dmax - a - b - c - ee):
                t = fm(t, qF)
            acc = (acc[0] + P([cf]) * t[0], acc[1] + P([cf]) * t[1])
        return acc[0].degree() < 0 and acc[1].degree() < 0

    P22, P33, P23, P32 = red.P[2][2], red.P[3][3], red.P[2][3], red.P[3][2]
    PW = red.Q[2] * P33 - red.Q[3] * P23
    PV = P22 * red.Q[3] - P32 * red.Q[2]
    return field_eval(red.Dpoly), field_eval(PW), field_eval(PV)


def check_sonic_vanishing(at):
    """(E5) every kappa-coefficient of U6 vanishes at the sonic point u0 in
    Q(V0, sqrt3) (U6 = U * denU with denU ~ xi^3 and U ~ xi^-2), and
    Aer_3(u0) = Ber_3(u0) = 0.  Justifies raising ALL nB kdiv inputs by one
    Laurent order (the top by two, with E1) in sonic_raise_inputs."""
    from flint import fmpq_poly as P
    red, mp, ctx = at.red, at.mp, at.red.ctx
    fm = lambda x, y: (x[0] * y[0] + 3 * (x[1] * y[1]), x[0] * y[1] + x[1] * y[0])
    Sp = (P([0]), P([1]))
    msV = (P([1]), P([0, -1]))
    qF = fm(fm((P([0]), P([4])), (P([1, 0, -1]), P([0]))), msV)
    pA = fm(fm(Sp, msV), (P([7, 0, -3]), P([0, 2])))
    pN = fm((P([0]), P([4])), fm((P([1, 0, -1]), P([0])), (P([0, -1]), P([1]))))
    pW = fm((P([fmpq(3, 2)]), P([0])), fm(msV, (P([0, -2]), P([1, 0, -1]))))
    pV = fm((P([0, 1]), P([0])), qF)

    def field_zero(f):
        if f.is_zero():
            return True
        acc = (P([0]), P([0]))
        dmax = sum(f.degrees()[i] for i in (1, 2, 3, 4))
        for e, cf in f.terms():
            a, b, c_, ee = int(e[1]), int(e[2]), int(e[3]), int(e[4])
            t = (P([1]), P([0]))
            for base, n in ((pA, a), (pN, b), (pW, c_), (pV, ee)):
                for _ in range(n):
                    t = fm(t, base)
            for _ in range(dmax - a - b - c_ - ee):
                t = fm(t, qF)
            acc = (acc[0] + P([cf]) * t[0], acc[1] + P([cf]) * t[1])
        return acc[0].degree() < 0 and acc[1].degree() < 0

    def kcoef(f, j):
        out = ctx.constant(0)
        for e, cf in f.terms():
            if int(e[9]) == j:
                ee = list(map(int, e))
                ee[9] = 0
                out += ctx.from_dict({tuple(ee): cf})
        return out

    EAc = 2 * mp["dDp"] * mp["Acr"] + mp["Dp"] * mp["dAcr"]
    EBc = 2 * mp["dDp"] * mp["Bcr"] + mp["Dp"] * mp["dBcr"]
    UN3r = EAc * mp["D2r"] + mp["SDr"] * (mp["Bcr"] * mp["Aer"]
                                          - mp["Acr"] * mp["Ber"])
    U6 = mp["Bcr"] * UN3r - mp["Acr"] * EBc * mp["D2r"]
    u6_all = all(field_zero(kcoef(U6, j)) for j in range(10))
    aer3 = field_zero(kcoef(mp["Aer"], 3))
    ber3 = field_zero(kcoef(mp["Ber"], 3))
    return u6_all, aer3, ber3


def _prod_jets(ak, names, depth, shifts=None):
    out = None
    for n in names:
        kp = ak.jet_kp(n, depth, num_shift=(shifts or {}).get(n, 0) or ak.sd_shift)
        out = kp if out is None else kpj_mul(out, kp, depth)
    return out


def build_pairs(ak, depth=2):
    """A0 = (-TN1r, denT), B0 = (-U6, denU) as jet kappa-polys (depth 2);
    denT = Dp SDr D2r Bcr, denU = Dp^3 SDr D2r^2 Bcr."""
    d = depth
    sh = ak.sd_shift
    TN = ak.jet_kp("TN1r", d, num_shift=sh)
    Dp = ak.jet_kp("Dp", d, num_shift=sh)
    SDr = ak.jet_kp("SDr", d, num_shift=sh)
    D2r = ak.jet_kp("D2r", d, num_shift=sh)
    Bcr = ak.jet_kp("Bcr", d, num_shift=sh)
    Acr = ak.jet_kp("Acr", d, num_shift=sh)
    Aer = ak.jet_kp("Aer", d, num_shift=sh)
    Ber = ak.jet_kp("Ber", d, num_shift=sh)
    dDp = ak.jet_kp("dDp", d, num_shift=sh)
    dAcr = ak.jet_kp("dAcr", d, num_shift=sh)
    dBcr = ak.jet_kp("dBcr", d, num_shift=sh)
    denT = kpj_mul(kpj_mul(Dp, SDr, d), kpj_mul(D2r, Bcr, d), d)
    EAc = kpj_add(kpj_scale(kpj_mul(dDp, Acr, d), arb(2)), kpj_mul(Dp, dAcr, d))
    EBc = kpj_add(kpj_scale(kpj_mul(dDp, Bcr, d), arb(2)), kpj_mul(Dp, dBcr, d))
    UN3r = kpj_add(kpj_mul(EAc, D2r, d),
                   kpj_mul(SDr, kpj_add(kpj_mul(Bcr, Aer, d),
                                        kpj_neg(kpj_mul(Acr, Ber, d))), d))
    U6 = kpj_add(kpj_mul(Bcr, UN3r, d), kpj_neg(kpj_mul(kpj_mul(Acr, EBc, d), D2r, d)))
    denU = kpj_mul(kpj_mul(kpj_mul(Dp, Dp, d), kpj_mul(Dp, SDr, d), d),
                   kpj_mul(kpj_mul(D2r, D2r, d), Bcr, d), d)
    return dict(nT=kpj_neg(TN), dT=denT, nB=kpj_neg(U6), dB=denU, TN=TN,
                _Dp0=Dp[0][0], _SDr0=SDr[0][0], _D2r2=D2r[2][0],
                _Bcr3=Bcr[3][0], _Bcr_kp=[x[0] for x in Bcr],
                _comps=(ak.A, ak.N, ak.W, ak.V),
                _akp=dict([(n, [x[0] for x in kp]) for n, kp in
                           (("Dp", Dp), ("SDr", SDr), ("D2r", D2r), ("Bcr", Bcr),
                            ("Acr", Acr), ("Aer", Aer), ("Ber", Ber), ("dDp", dDp),
                            ("dAcr", dAcr), ("dBcr", dBcr))]
                          + [("_A", ak.A), ("_W", ak.W)]),
                _ajets=dict((n, kp) for n, kp in
                            (("Dp", Dp), ("SDr", SDr), ("D2r", D2r), ("Bcr", Bcr),
                             ("Acr", Acr), ("Aer", Aer), ("Ber", Ber),
                             ("dDp", dDp), ("dAcr", dAcr), ("dBcr", dBcr))))


def commute(pr, dTtop_inv):
    """One commutation level: returns dict with A1 = (nA1, dA1) [depth 1],
    B1 = (nB1, dB1) [depth 0], tt jet [depth 1], and W = A1 + kappa*tt as
    (nW [top degree dropped: the exact cancellation tt*dT_top = TN1r_top], dA1).
    dTtop_inv: LC inverse of denT's top coefficient (offset-raised by the
    caller at the sonic point, justified by D(u0) = 0)."""
    nT, dT, nB, dB = pr["nT"], pr["dT"], pr["nB"], pr["dB"]
    d1 = lambda p: kpj_depth(p, 1)
    d0 = lambda p: kpj_depth(p, 0)
    # dlogB0 = (nB' dB - nB dB') / (nB dB), depth 1
    ln_num = kpj_add(kpj_mul(kpj_deriv(nB), d1(dB), 1),
                     kpj_neg(kpj_mul(d1(nB), kpj_deriv(dB), 1)))
    ln_den = kpj_mul(d1(nB), d1(dB), 1)
    # A1 = A0 - dlogB0 = (nT*ln_den - ln_num*dT) / (dT*ln_den), depth 1
    nA1 = kpj_add(kpj_mul(d1(nT), ln_den, 1), kpj_neg(kpj_mul(ln_num, d1(dT), 1)))
    dA1 = kpj_mul(d1(dT), ln_den, 1)
    # B1 = B0 + A0' - A0*dlogB0 over dT^2 * nB * dB, depth 0
    A0p_num = kpj_add(kpj_mul(kpj_deriv(nT), d0(dT), 0),
                      kpj_neg(kpj_mul(d0(nT), kpj_deriv(dT), 0)))
    nB1 = kpj_add(kpj_mul(kpj_mul(d0(nB), d0(nB), 0), kpj_mul(d0(dT), d0(dT), 0), 0),
                  kpj_mul(A0p_num, kpj_mul(d0(nB), d0(dB), 0), 0),
                  kpj_neg(kpj_mul(kpj_mul(d0(nT), d0(ln_num), 0), d0(dT), 0)))
    dB1 = kpj_mul(kpj_mul(d0(dT), d0(dT), 0), kpj_mul(d0(nB), d0(dB), 0), 0)
    # tt = kappa-lead of T = TN1r_top / denT_top (real), jet depth 1, built
    # with the raised inverse (T ~ tt kappa with tt ~ 1/xi at the sonic point)
    TNt, dTt = pr["TN"][-1], dT[-1]
    ti = dTtop_inv
    tt0 = TNt[0] * ti
    tt1 = (TNt[1] * dTt[0] - TNt[0] * dTt[1]) * (ti * ti)
    tt = [tt0, tt1]
    # W = A1 + kappa*tt: the top coefficient cancels EXACTLY (numA1's top is
    # -TN1r_top*(nB dB)_top and kappa*tt*dA1's top is +tt*dT_top*(nB dB)_top
    # with tt*dT_top == TN1r_top by the definition of tt), so it is omitted.
    ttden = [jlc_mul(tt, x[:2], 1) for x in dA1]
    nW = kpj_add(d1(nA1), kpj_kmul(kpj_depth(ttden, 1)))
    nW = nW[:-1]                                     # exact top cancellation
    # W' = (nW' dA1 - nW dA1') / dA1^2, depth 0 (for the nu' bound)
    nWp = kpj_add(kpj_mul(kpj_deriv(nW), d0(dA1), 0),
                  kpj_neg(kpj_mul(d0(nW), kpj_deriv(dA1), 0)))
    dA1sq = kpj_mul(d0(dA1), d0(dA1), 0)
    return dict(nA1=nA1, dA1=dA1, nB1=nB1, dB1=dB1, tt=tt, nW=nW, nWp=nWp,
                dA1sq=dA1sq, ln_num=ln_num, ln_den=ln_den)


def kdiv_lc(num, den, lead_inv, gap):
    """num = q(kappa)*den + rem at the LC level; q of degree = gap; lead_inv
    must be an LC inverse of den[-1] (possibly offset-raised by the caller,
    justified by the exact vanishing identities).  Returns (q asc, rem)."""
    nd = len(den) - 1
    R = list(num) + [lc_zero(den[0].tm.R)] * (nd + gap + 1 - len(num))
    q = [None] * (gap + 1)
    for j in range(nd + gap, nd - 1, -1):
        c = R[j] * lead_inv
        q[j - nd] = c
        for i in range(nd):
            R[i + j - nd] = R[i + j - nd] - c * den[i]
        R[j] = lc_zero(den[0].tm.R)
    return q, R[:nd]


def slot0(kp):
    return [x[0] for x in kp]


def raise_jet(jlc, k0):
    """Raise the Laurent offsets of a jet: slot i by max(0, k0 - i) (each x-
    derivative loses one vanishing order).  Sound via drop_lead's exact-zero /
    ve checks; the caller certifies the true orders (exact identities)."""
    return [lc_raise(x.trim(), max(0, k0 - i)) if max(0, k0 - i) else x
            for i, x in enumerate(jlc)]


def sonic_raise_inputs(pr):
    """Certified Laurent orders of the level-0 pair data at the sonic point
    (xi = -x -> 0): every dB = Dp^3 SDr D2r^2 Bcr coefficient vanishes to
    order 3 and every dT = Dp SDr D2r Bcr coefficient to order 1 (the exact
    Dp factor, D(u0) = 0: E3); nB = -U6 has orders >= (0,...,0,1,1,2) in
    kappa: U6_9 = Dp^2 * unit (E1) and U6_8(u0) = U6_7(u0) = 0 in Q(V0,sqrt3)
    (E5).  Raising the inputs keeps every later kdiv aligned and junk-free."""
    pr = dict(pr)
    pr["dB"] = [raise_jet(x, 3) for x in pr["dB"]]
    pr["dT"] = [raise_jet(x, 1) for x in pr["dT"]]
    # every U6_j(u0) = 0 in Q(V0, sqrt3) (checked exactly: U6 = U * denU with
    # denU ~ xi^3 and U ~ xi^{-2}); U6_9 additionally = Dp^2 * unit (E1)
    pr["nB"] = [raise_jet(x, 2 if j == 9 else 1)
                for j, x in enumerate(pr["nB"])]
    return pr


def slot1(kp):
    return [x[1] for x in kp]


def kdiv_jlc(num, den, lead_inv_j, gap, depth=1):
    """kdiv at jet depth (<= 2): num, den kappa-polys of JLCs; lead_inv_j =
    jet of 1/den_top.  Returns (quotient JLCs asc, remainder JLCs)."""
    R0 = den[0][0].tm.R
    d1 = depth + 1
    zj = [lc_zero(R0)] * d1
    nd = len(den) - 1
    R = [x[:d1] for x in num] + [zj] * (nd + gap + 1 - len(num))
    q = [None] * (gap + 1)
    for j in range(nd + gap, nd - 1, -1):
        c = jlc_mul(R[j], lead_inv_j[:d1], depth)
        q[j - nd] = c
        for i in range(nd):
            t = jlc_mul(c, den[i][:d1], depth)
            R[i + j - nd] = [R[i + j - nd][k] - t[k] for k in range(d1)]
        R[j] = zj
    return q, R[:nd]


class Bnd:
    """Upper bound of a nonnegative quantity as sum of val * ell^off terms."""

    def __init__(self, terms=()):
        self.t = list(terms)

    @classmethod
    def of(cls, val, off=0):
        return cls([(arb(val) if not isinstance(val, arb) else val, off)])


    def __add__(self, o):
        return Bnd(self.t + o.t)

    def scale(self, c, doff=0):
        return Bnd([(v * c, o + doff) for v, o in self.t])

    def mul(self, o):
        return Bnd([(v1 * v2, o1 + o2) for v1, o1 in self.t for v2, o2 in o.t])

    def eval_up(self, lo, hi, scale_off=0):
        acc = arb(0)
        el = arb((lo + hi) / 2, (hi - lo) / 2)
        for v, o in self.t:
            oo = o + scale_off
            if oo < 0 and lo <= 0:
                raise ArithmeticError("negative bound offset at ell = 0 "
                                      "(o=%d scale=%d v=%s terms=%d)"
                                      % (o, scale_off, str(v)[:16], len(self.t)))
            acc += arb(v) * abs(el) ** oo
        return acc.upper()

    def min_off(self):
        return min((o for v, o in self.t if not (v == 0)), default=0)


def lc_abs(lc, lo, hi):
    cont, junk = lc.split_junk()
    out = Bnd.of(abs(cont.tm.eval(lo, hi)).upper(), cont.off)
    if junk is not None:
        for k, x in enumerate(junk.tm.c):
            if not x.is_zero():
                out = out + Bnd.of(abs(x).upper(), junk.off + k)
    return out


def kp_dk(kp):
    """kappa-derivative of a kappa-poly of LCs (ascending)."""
    return [j * kp[j] for j in range(1, len(kp))]


def kp_mul(a, b):
    """product of two kappa-polys of LCs (ascending)."""
    out = [None] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            t = x * y
            out[i + j] = t if out[i + j] is None else out[i + j] + t
    return [x.trim() for x in out]


def kp_sub(a, b):
    n = max(len(a), len(b))
    out = []
    for j in range(n):
        if j >= len(b):
            out.append(a[j])
        elif j >= len(a):
            out.append(-1 * b[j])
        else:
            out.append(a[j] + (-1 * b[j]))
    return [x.trim() for x in out]


def lc_split_lead(lc, k):
    """Exact split of an LC into (first k slots, rest) -- both enclose their
    own parts; the sum is the original."""
    from .tc_tm import TM, LC
    lead = LC(TM([x for x in lc.tm.c[:k]], lc.tm.R), lc.off)
    rest = LC(TM([x for x in lc.tm.c[k:]], lc.tm.R), lc.off + k)
    return lead, rest


def lc_hat(lc):
    """(content LC, junk Bnd or None) of a quotient/coefficient LC."""
    cont, junk = lc.split_junk()
    jb = None
    if junk is not None:
        for k, x in enumerate(junk.tm.c):
            if not x.is_zero():
                b = Bnd.of(abs(x).upper(), junk.off + k)
                jb = b if jb is None else jb + b
    return cont, jb


class DenLower2:
    """|den(kappa)(ell)| >= C * ell^doff * tau^g on the piece, |Im kappa| = tau
    >= klo, Re kappa in [0, EPS_MAX]."""

    def __init__(self, C=None, doff=0, g=0):
        self.C = C if C is not None else arb(1)
        self.doff, self.g = doff, g

    def mul_lc(self, lc, lo, hi, power=1, tau_order=0):
        v = abs(lc.tm.eval(lo, hi)).lower()
        if not (v > 0):
            raise ArithmeticError("factor TM lower bound not positive")
        return DenLower2(self.C * arb(v) ** power, self.doff + lc.off * power,
                         self.g + tau_order * power)

    def scale(self, c):
        c = arb(c)
        if not (c > 0):
            raise ArithmeticError("nonpositive scale in DenLower2")
        return DenLower2(self.C * c, self.doff, self.g)

    def mul(self, o):
        return DenLower2(self.C * o.C, self.doff + o.doff, self.g + o.g)


def quad_lower_mixed(q2, q1, q0, klo, lo, hi):
    """DenLower2 for the B0-quadratic  q(k) = q2 k^2 + q1 k + q0  on the strip
    tail (k = eps + i tau, eps in [0, EPS_MAX], tau >= klo):
        |q| >= (|Re q| + |Im q|)/2 >= [|q2| tau^2 + |q1| tau - X]/2,
        X := |q2| eps^2 + |q1| eps + |q0|   (so |Re q| >= |q2| tau^2 - X),
        |Im q| = tau |2 eps q2 + q1| >= tau |q1|   (q2, q1 of one sign, eps>=0).
    If X <= (|q2| klo^2)/2 + |q1| klo pointwise, then |q| >= |q2| tau^2 / 4."""
    (q2, j2), (q1, j1), (q0, j0) = lc_hat(q2), lc_hat(q1), lc_hat(q0)
    s2, s1 = q2.tm.eval(lo, hi), q1.tm.eval(lo, hi)
    if not ((s2 < 0 and s1 < 0) or (s2 > 0 and s1 > 0)):
        raise ArithmeticError("mixed quad-lower: q2, q1 not of one sign")
    emax, kl = arb(EPS_MAX), arb(klo)
    X = (lc_abs(q2, lo, hi).scale(emax * emax) + lc_abs(q1, lo, hi).scale(emax)
         + lc_abs(q0, lo, hi))
    # near-zero junk parts of the quotient coefficients degrade |q| by at most
    # J tau^2 with J = (sum of junk hulls) * (1 + emax^2/klo^2); subtract from C
    J = arb(0)
    for jb in (j2, j1, j0):
        if jb is not None:
            J += arb(jb.eval_up(max(lo, 0.0), hi, 0)) if all(
                o >= 0 for _, o in jb.t) else arb(jb.eval_up(lo, hi, 0))
    J = J * (1 + emax * emax / (kl * kl))
    omin = min(q2.off, q1.off, X.min_off())
    el = arb((lo + hi) / 2, (hi - lo) / 2)
    Xup = X.eval_up(lo, hi, -omin)
    p2 = abs(el) ** (q2.off - omin)
    p1 = abs(el) ** (q1.off - omin)
    rhs = (abs(s2).lower() * (p2.lower() if q2.off > omin else 1) * kl * kl / 2
           + abs(s1).lower() * (p1.lower() if q1.off > omin else 1) * kl)
    if not (arb(Xup) < rhs):
        raise ArithmeticError("mixed quad-lower: X exceeds |q2| klo^2/2 + |q1| klo")
    el2 = abs(arb((lo + hi) / 2, (hi - lo) / 2)) ** q2.off
    C = (abs(s2).lower() / 4 - J / el2.lower() if q2.off else abs(s2).lower() / 4 - J)
    C = arb(C).lower()
    if not (C > 0):
        raise ArithmeticError("mixed quad-lower: |q2|/4 - junk not positive")
    return DenLower2(arb(C), q2.off, 2)


def rem_ratio_bound(rem, dlow, klo, lo, hi):
    """(Bnd, decay): pointwise bound |rem(kappa)(ell)| / |den(kappa)(ell)| <=
    Bnd(ell) * (klo/tau)^decay for ell in [lo, hi], tau >= klo (strip tail),
    |den| >= dlow.C ell^dlow.doff tau^dlow.g; per-coefficient offsets are kept
    (junk balls are split off so they do not drag the offsets down).
    decay = dlow.g - deg(rem) >= 0 required."""
    jmax = len(rem) - 1
    while jmax > 0 and all(x.is_zero() for x in rem[jmax].tm.c):
        jmax -= 1
    decay = dlow.g - jmax
    if decay < 0:
        raise ArithmeticError("remainder degree exceeds denominator growth")
    kl, emax = arb(klo), arb(EPS_MAX)
    out = None
    for j, r in enumerate(rem[:jmax + 1]):
        fac = (1 + emax * emax / (kl * kl)).sqrt() ** j / (kl ** (dlow.g - j) * dlow.C)
        cont, jb = lc_hat(r)
        hull = abs(cont.tm.eval(lo, hi))
        if hull.upper() != 0:
            t = Bnd.of((hull * fac).upper(), cont.off - dlow.doff)
            out = t if out is None else out + t
        if jb is not None:
            t = jb.scale(fac.upper(), -dlow.doff)
            out = t if out is None else out + t
    return (out if out is not None else Bnd.of(arb(0), 0)), decay


def _hat_kp(kp, slo, shi):
    """(hulls, minoff): per-coefficient TM hulls times ell^{off-minoff}
    (monotone endpoint powers), stripping the common ell^minoff."""
    kt = [c.trim() for c in kp]
    offs = [c.off for c in kt if not all(x.is_zero() for x in c.tm.c)]
    m = min(offs) if offs else 0
    out = []
    ehi, elo = arb(max(shi, 0.0)), arb(max(slo, 0.0))
    for c in kt:
        d = c.off - m
        if d < 0:
            raise ArithmeticError("hat offset below the minimum")
        pw = arb(1) if d == 0 else (elo ** d).union(ehi ** d)
        out.append(c.tm.eval(slo, shi) * pw)
    return out, m


def AFN_hulls(akp, slo, shi):
    """(A, F_N) hulls on the strip from the component LCs stashed in akp."""
    Alc, Wlc = akp["_A"], akp["_W"]
    Ah = Alc.eval(slo, shi)
    FNh = Ah - 2 + (arb(-2) / 3) * Wlc.eval(slo, shi)
    return Ah, FNh


def strip_list(lo, hi, nstrip):
    if lo <= 0:
        return ([(0.0, hi * 0.5 ** (nstrip - 1))] +
                [(hi * 0.5 ** (k + 1), hi * 0.5 ** k)
                 for k in reversed(range(nstrip - 1))])
    rr = (hi / lo) ** (1.0 / nstrip)
    return [(lo * rr ** k, lo * rr ** (k + 1)) for k in range(nstrip)]


def sup_ratio_boxes(rems, pr, uq, klo, lo, hi, thU=None, nbox=26, nstrip=10,
                    neps=6, diag=None, dp_raise=0, net_offs=None):
    """net_offs: name -> integer o; the returned bounds then certify
    |rem/den| <= D * ell^o (and <= Ddec ell^o / tau), allowing negative true
    Laurent orders at a closing piece where the conditions are ell^2-scaled."""
    """rems: name -> (rem kp, (pT, pU, pB)): pointwise bounds over the strip
    tail (Re kappa in [0, EPS_MAX], tau >= klo) and ell in [lo, hi] of
    |rem(kappa)(ell)| / (|dT|^pT |U6|^pU |dB|^pB)(kappa)(ell) <= D, <= Ddec/tau.
    U6 is lowered through the exact B0-quadratic dominance |U6| >= (1 - thU)
    |quad| |dB| (quad = kdiv quotient of B0, thU certified first with pU = 0
    families).  dT, dB, quad are box-evaluated from tight atom hats (D2r via
    its exact linear factors); kappa-boxes are geometric in tau ending in one
    pseudo-box w in [0, 1/T] that covers tau = infinity."""
    from flint import acb
    akp = pr["_akp"]
    prep = {}
    for name, (rem, pows) in rems.items():
        remt = [x.trim() for x in rem]
        pows = tuple(pows) + (0,) * (4 - len(pows))
        dendeg = 5 * pows[0] + 9 * pows[1] + 7 * pows[2] + 2 * pows[3]
        gap = dendeg - (len(remt) - 1)
        if gap < 0:
            raise ArithmeticError("negative decay reserve for " + name)
        if pows[1] and thU is None:
            raise ArithmeticError("U6 power needs a certified thU")
        prep[name] = (remt, pows, gap, dendeg)
    net_off = dict((n, 0) for n in rems) if net_offs is None else dict(net_offs)
    edges = [klo * 1.15 ** k for k in range(14)] +             [klo * 1.15 ** 13 * 2.0 ** (k + 1) for k in range(nbox)]
    out = {name: [arb(0), arb(0), arb(0)] for name in rems}
    one_m_th = arb(1) if thU is None else (1 - arb(thU))
    if not (one_m_th > 0):
        raise ArithmeticError("thU >= 1")
    for slo, shi in strip_list(lo, hi, nstrip):
        H, M = {}, {}
        for n, kp in akp.items():
            if n.startswith("_"):
                continue
            if n == "Dp" and dp_raise:
                # Dp = -xi E_D at the sonic point (D(u0) = 0 exact, E3):
                # evaluate the raised hat and carry the offset
                kpr = [lc_raise(x.trim(), dp_raise) for x in kp]
                H[n], M[n] = _hat_kp(kpr, slo, shi)
            else:
                H[n], M[n] = _hat_kp(kp, slo, shi)
        Ah, FNh = AFN_hulls(akp, slo, shi)
        uqc, uqj_ = [], []
        for q in uq:
            cq, jq = lc_hat(q)
            uqc.append(cq)
            uqj_.append(jq)
        oQ = min(q.off for q in uqc)

        def elp(d):
            if d == 0:
                return arb(1)
            return (arb(max(slo, 0.0)) ** d).union(arb(max(shi, 0.0)) ** d)

        ncs = {name: _hat_kp(remt, slo, shi) for name, (remt, _, _, _) in prep.items()}
        nBh, moffB = _hat_kp(slot0(pr["nB"]), slo, shi)
        for i in range(neps):
            e0, e1 = EPS_MAX * i / neps, EPS_MAX * (i + 1) / neps
            boxes = list(zip(edges, edges[1:])) + [(edges[-1], None)]
            for t0, t1 in boxes:
                if t1 is None:
                    wv = acb(arb(0, 1.0 / t0), arb(0, 1.0 / t0))
                    wlo = wilo = arb(0)
                else:
                    z = acb(arb((e0 + e1) / 2, (e1 - e0) / 2),
                            arb((t0 + t1) / 2, (t1 - t0) / 2))
                    wv = acb(1) / z
                    wlo = arb(1) / arb(abs(z).upper())
                    wilo = arb(t0) / (arb(e1) ** 2 + arb(t1) ** 2)
                # factored den hats
                Dp, SDr, dDp = H["Dp"][0], H["SDr"][0], H["dDp"][0]
                D2rf = H["D2r"][2] * (1 - Ah * wv) * (1 + FNh * wv)
                Bcr = None
                for x in H["Bcr"]:
                    Bcr = x + 0 * wv if Bcr is None else Bcr * wv + x
                lT = abs(Dp * SDr * D2rf * Bcr).lower()
                lB = abs(Dp * Dp * Dp * SDr * D2rf * D2rf * Bcr).lower()
                t2 = uqc[2].tm.eval(slo, shi) * elp(uqc[2].off - oQ)
                t1_ = (uqc[1].tm.eval(slo, shi) * elp(uqc[1].off - oQ)) * wv
                t0_ = (uqc[0].tm.eval(slo, shi) * elp(uqc[0].off - oQ)) * wv * wv
                jup = arb(0)
                for k_, jq in enumerate(uqj_):
                    if jq is not None:
                        jup += arb(jq.eval_up(slo, shi, -oQ)) / arb(t0) ** (2 - k_)
                quad = t2 + t1_ + t0_
                lQ = (abs(quad).lower() - jup).lower()
                # fallback where the ell-hull of the q2 term contains 0:
                # |quad| >= |q1 w| - |q2 t| - |q0 w^2| - junk
                # |q1 w| = |q1| |w|: product of abs-lowers survives wide
                # kappa boxes where the acb hull of t1_ self-cancels
                l1 = (abs(uqc[1].tm.eval(slo, shi)
                          * elp(uqc[1].off - oQ)).lower() * wlo)
                alt = (l1 - abs(t2).upper() - abs(t0_).upper() - jup).lower()
                # the q-hats are real, so t2 never touches Im(quad):
                # |quad| >= |q1| Im(w)_lo - |t0_| - junk
                alt2 = (abs(uqc[1].tm.eval(slo, shi)
                            * elp(uqc[1].off - oQ)).lower() * wilo
                        - abs(t0_).upper() - jup).lower()
                lQ = max(lQ, alt, alt2)
                if not (lT > 0 and lB > 0):
                    raise ArithmeticError(
                        "den box lower not positive (strip=[%g,%g] t0=%g "
                        "lT=%s lB=%s)" % (slo, shi, t0, lT, lB))
                # branches: (lQ lower, extra ell-powers, extra kappa-decay)
                branches = [(arb(lQ), 0, 0)] if lQ > 0 else []
                if lQ <= 0 and t1 is None and not (jup > 0):
                    # pseudo-box over an ell = 0 strip: both hulls contain 0,
                    # so no constant lower for |quad| exists.  Signs align
                    # (Re w >= 0, Im w <= 0, q-hats real of equal sign), so
                    # |quadN| >= (l2/sqrt2) ell + cI |w| holds pointwise;
                    # each summand is certified as its own branch below.
                    s2v = uqc[2].tm.eval(slo, shi)
                    s1v = uqc[1].tm.eval(slo, shi) * elp(uqc[1].off - oQ)
                    s0v = uqc[0].tm.eval(slo, shi) * elp(uqc[0].off - oQ)
                    if (s2v * s1v) > 0:
                        rt2 = arb(2).sqrt()
                        fI = (1 - arb(EPS_MAX) ** 2 / arb(t0) ** 2).sqrt()
                        cI = (abs(s1v).lower() * fI.lower() / rt2
                              - abs(s0v).upper() / arb(t0)).lower()
                        lQR = (abs(s2v).lower() / rt2).lower()
                        if lQR > 0:
                            branches.append((arb(lQR), 1, 0))
                        if cI > 0:
                            branches.append((arb(cI), 0, 1))
                if not branches:
                    raise ArithmeticError(
                        "quad box lower not positive (strip=[%g,%g] t0=%g "
                        "t2=%s t1_=%s jup=%s)" % (slo, shi, t0, t2, t1_, jup))
                # direct |U6| lower via top dominance: |nB(kappa)| >=
                # ell^moffB t0^9 (|top| - sum_j |c_j|/t0^{9-j}); no 1/tau loss
                # where quad-hat ~ u1 w (the ell = 0 strip at huge tau)
                nBtop = len(nBh) - 1
                sub6 = arb(0)
                for j in range(nBtop):
                    sub6 += abs(nBh[j]).upper() / arb(t0) ** (nBtop - j)
                lU6 = (abs(nBh[nBtop]).lower() - sub6).lower()
                if not (lU6 > 0):
                    lU6 = None
                oT = M["Dp"] + M["SDr"] + M["D2r"] + M["Bcr"]
                oB = 3 * M["Dp"] + M["SDr"] + 2 * M["D2r"] + M["Bcr"]
                for name, (nc, moff) in ncs.items():
                    remt, pows, gap, dendeg = prep[name]
                    nvw = acb(0)
                    for x in nc:
                        nvw = nvw * wv + x
                    nq = pows[1] + pows[3]
                    rD, rDec, rG = None, None, None
                    brs = [(lQb, dell, dgap, False) for lQb, dell, dgap
                           in branches]
                    if lU6 is not None and pows[1] and not pows[3]:
                        brs.append((None, 0, 0, True))
                    for lQb, dell, dgap, u6 in brs:
                        if u6:
                            dl = (arb(lT) ** pows[0] * arb(lU6) ** pows[1]
                                  * arb(lB) ** pows[2]).lower()
                            net = moff - (oT * pows[0] + moffB * pows[1]
                                          + oB * pows[2]) - net_off[name]
                            gap_ = gap
                        else:
                            dl = (arb(lT) ** pows[0]
                                  * (one_m_th * lQb * arb(lB)) ** pows[1]
                                  * arb(lB) ** pows[2] * lQb ** pows[3]).lower()
                            # branch-R replaces lQ ell^oQ by lQR ell^{oQ+1}:
                            # the den GAINS ell-powers, so net DECREASES
                            net = moff - (oT * pows[0] + (oQ + oB) * pows[1]
                                          + oB * pows[2] + oQ * pows[3]) \
                                - net_off[name] - dell * nq
                            gap_ = gap - dgap * nq
                        if not (dl > 0) or (net < 0 and slo <= 0) or gap_ < 0:
                            continue
                        nel = abs(elp(net)).upper() if net >= 0 else \
                            (arb(slo) ** net).upper()
                        r = abs(nvw).upper() * nel / arb(dl)
                        rD_ = (arb(r) / arb(t0) ** gap_).upper()
                        if t1 is None:
                            rDec_ = (arb(r) / arb(t0) ** (gap_ - 1)).upper() \
                                if gap_ >= 1 else None
                        else:
                            rDec_ = (arb(r) * t1 / arb(t0) ** gap_).upper()
                        if rD is None or float(rD_) < float(rD):
                            rD = rD_
                        if rDec_ is not None and (rDec is None
                                                  or float(rDec_) < float(rDec)):
                            rDec = rDec_
                        # sup of tau^gap * ratio over this box (valid only
                        # where the branch's tau-power deficit is coverable)
                        if gap_ == gap:
                            rG_ = arb(r)
                        elif t1 is not None:
                            rG_ = arb(r) * arb(t1) ** (gap - gap_)
                        else:
                            rG_ = None
                        if rG_ is not None and (rG is None
                                                or float(rG_) < float(rG)):
                            rG = rG_
                    if rD is None:
                        raise ArithmeticError(
                            "no valid quad branch for %s (strip=[%g,%g] "
                            "t0=%g moff=%d gap=%d)"
                            % (name, slo, shi, t0, moff, gap))
                    if diag is not None and rDec is not None and \
                            float(rDec) > max((d[3] for d in diag.get(name, [])),
                                              default=0.0):
                        diag.setdefault(name, []).append(
                            (slo, t0, float(rD), float(rDec)))
                    out[name][0] = out[name][0].max(arb(rD))
                    if rDec is not None:
                        out[name][1] = out[name][1].max(arb(rDec))
                    # sup of tau^gap * ratio (for Im tail-integration)
                    if rG is None:
                        rG = arb(float("inf"))
                    out[name][2] = out[name][2].max(arb(rG))
    return {name: (v[0].upper(), v[1].upper(), v[2].upper())
            for name, v in out.items()}


def cubic_root_bound(bcr_kp, lo, hi):
    """Cauchy: every root z of the Bcr cubic has |z| <= 1 + max_j |c_j/c_3|,
    uniformly on the subrange (needs matching Laurent offsets or lower ones on
    the lead; compares at the lead's offset)."""
    bcr_kp = [x.trim() for x in bcr_kp]
    c3 = bcr_kp[3]
    l3 = abs(c3.tm.eval(lo, hi)).lower()
    if not (l3 > 0):
        raise ArithmeticError("Bcr leading coefficient not bounded away from 0")
    el = arb((lo + hi) / 2, (hi - lo) / 2)
    m = arb(0)
    for j in range(3):
        cj = bcr_kp[j]
        if cj.off < c3.off and lo <= 0:
            raise ArithmeticError("Bcr coefficient offset below the lead at 0")
        ratio = abs(cj.tm.eval(lo, hi)) * abs(el) ** (cj.off - c3.off) / l3
        m = m.max(ratio)
    return (1 + m).upper()


def analyze_piece(pr, cm, lo, hi, b1, eps_cells, scale_off, raises,
                  dp_raise=None, tq=True):
    from .tc_tm import SPLIT_REL
    SPLIT_REL[0] = 1e-4 if (lo <= 0 and scale_off) else 1e-12
    """Certified strip-exclusion conditions on one local subrange [lo, hi]:
    B0 != 0 and Bcr != 0 on the tail (commutation + eta-elimination legal), and
    sup_{tau >= b1} [K2 tau^2 + R1 tau + K0 + R0] < 0 per eps cell.  raises =
    (r_dT, r_nB, r_dB): Laurent raises of the kdiv leads at the sonic point,
    justified by D(u0) = 0 (E3) and the U6_9 factorization (E1)."""
    c, cl = cval(), abs(cval())
    b1a = arb(b1)
    nT, dT, nB, dB = (slot0(pr[k]) for k in ("nT", "dT", "nB", "dB"))
    nA1, dA1 = slot0(cm["nA1"]), slot0(cm["dA1"])
    nB1, dB1 = slot0(cm["nB1"]), slot0(cm["dB1"])
    nW = slot0(cm["nW"])
    r_dT, r_nB, r_dB = raises
    closing0 = bool(scale_off)          # Laurent-sonic treatment
    inv_of = lambda lc, k: (lc_raise(lc.trim(), k) if k else lc).inv()
    dBi = inv_of(dB[-1], r_dB)
    dA1i = inv_of(dA1[-1], r_dT + r_nB + r_dB)
    dB1i = inv_of(dB1[-1], 2 * r_dT + r_nB + r_dB)
    trim = lambda xs: [x.trim() for x in xs]
    dBtop_j = pr["dB"][-1][:3]
    dBi1 = -1 * (dBtop_j[1] * (dBi * dBi))
    dBi_j = [dBi, dBi1,
             2 * (dBtop_j[1] * dBtop_j[1] * (dBi * dBi * dBi))
             - dBtop_j[2] * (dBi * dBi)]
    uqj, remUj = kdiv_jlc(pr["nB"], pr["dB"], dBi_j, 2, depth=2)
    uq = [x[0].trim() for x in uqj]
    remU = [x[0].trim() for x in remUj]
    bq, remB1 = map(trim, kdiv_lc(nB1, dB1, dB1i, 2))
    Bkp = remB1
    if closing0 and bq[0].off <= -3:
        # b0 carries a genuine ell^{-3} hat (-tt*lambda1 in the kappa
        # asymptotics) that the remainder cancels at fixed kappa.  Split it
        # off and keep the cancellation:  b0l + remB1/dB1 = G/dB1 with
        # G := b0l*dB1 + remB1, whose ell^3 layers cancel exactly (the
        # per-slot raise is containment-checked); Im(G/dB1) = Im(remB1/dB1).
        b0l, b0h = lc_split_lead(bq[0], -2 - bq[0].off)
        tgt = b0l.off + dB1[-1].off
        Bkp = [(b0l * x + r).trim() for x, r in zip(dB1[:-1], remB1)] \
            + [(b0l * dB1[-1]).trim()]
        Bkp = [lc_raise(x, max(0, tgt - x.off)) for x in Bkp]
        bq = [b0h.trim(), bq[1], bq[2]]
    dA1top_j = cm["dA1"][-1][:2]
    dA1i_j = [dA1i, -1 * (dA1top_j[1].trim() * (dA1i * dA1i))]
    nWk = cm["nW"]
    if closing0:
        # A1 has a simple pole at xi = 0 and dA1-top = (dT nB dB)-tops has
        # val exactly 5 (pure atom product), so nW-top -- the w0 numerator --
        # has val 4 and its derivative val 3; the sub-val slots are
        # cancellation noise (drop_lead checks 0-containment)
        nWk = list(cm["nW"])
        nWk[-1] = [lc_raise(x.trim(), 1) for x in nWk[-1]]
    wqj, remWj = kdiv_jlc(nWk, cm["dA1"], dA1i_j, 0)
    wq = [x[0].trim() for x in wqj]
    remW = [x[0].trim() for x in remWj]
    remWp = [x[1].trim() for x in remWj]
    u0_, u1_, u2_ = uq
    b0_, b1_, b2_ = bq
    w0_ = wq[0]
    tt, ttp = cm["tt"][0].trim(), cm["tt"][1].trim()
    # ---- eta-elimination root bound and real-pole distances
    Rbar = cubic_root_bound(pr["_Bcr_kp"], lo, hi)
    b1a = arb(b1)
    if not (arb(Rbar) < b1a):
        raise ArithmeticError("Bcr root bound reaches the strip tail")
    Alc, Wlc = pr["_akp"]["_A"], pr["_akp"]["_W"]
    Aup = abs(Alc.eval(lo, hi)).upper()
    FNup = abs((Alc - 2 + (arb(-2) / 3) * Wlc).eval(lo, hi)).upper()
    if not (arb(Aup) < b1a and arb(FNup) < b1a):
        raise ArithmeticError("real pole bound reaches the strip tail")
    # ---- remainder bounds for tau >= tstar (thU first: B0 quad dominance)
    tstar = TSTAR_FAC * b1
    dpr = r_dT if dp_raise is None else dp_raise
    nof = dict(U=0)
    if not tq:
        out = dict(Rbar=float(Rbar), cells=[], ok=True,
                   tau_quad="open (tau > tstar not certified on this piece)")
        for e0, e1 in eps_cells:
            try:
                okd, wd, dfails = direct_S_boxes(pr, cm, uqj, remUj, lo, hi,
                                                 b1, tstar, [(e0, e1)],
                                                 scale_off)
            except ArithmeticError as e:
                okd, wd, dfails = False, None, []
            m = dict(direct=wd)
            if dfails:
                m["direct_fails"] = dfails
            out["cells"].append(dict(eps=[e0, e1], ok=okd, margin=wd, **m))
            out["ok"] = out["ok"] and okd
        return out
    thU_b = sup_ratio_boxes(dict(U=(remU, (0, 0, 1, 1))), pr, uq, tstar, lo, hi,
                            dp_raise=dpr, net_offs=nof)
    thU = thU_b["U"][0]
    if not (arb(thU) < 1):
        raise ArithmeticError("B0 quadratic domination fails (thU >= 1)")
    diag = {}
    # per-name Laurent allowances at the closing piece (true ratio orders:
    # nu-rem ~ ell^-1, nu'-rem ~ ell^-2, G/dB1 ~ ell^-2, log-derivs ~ ell^-1)
    noW, noWp, noB, LLo = (-1, -2, -2, -1) if closing0 else (0, 0, 0, 0)
    bnds = sup_ratio_boxes(
        dict(W=(remW, (1, 1, 1)), Wp=(remWp, (1, 1, 1)), B=(Bkp, (2, 1, 1)),
             LT=(slot1(pr["dT"]), (1, 0, 0)), LU=(slot1(pr["nB"]), (0, 1, 0)),
             LB=(slot1(pr["dB"]), (0, 0, 1))),
        pr, uq, tstar, lo, hi, thU=thU, diag=diag, dp_raise=dpr,
        net_offs=dict(W=noW, Wp=noWp, B=noB, LT=LLo, LU=LLo, LB=LLo))
    DW, DWdec = bnds["W"][:2]
    DWp, DWpdec = bnds["Wp"][:2]
    DB, DBdec = bnds["B"][:2]
    DWKI = DWPK = None
    if closing0:
        # Im-parts of the real-coefficient ratios RW = remW/dA1 and its
        # xi-derivative decay like 1/tau even where |RW| does not (the ell = 0
        # hat):  Im F(eps + i tau) = -int_tau^oo d/dt Im F dt, |d/dt F| =
        # |F'(kappa)|, so  |Im F| <= (sup tau^2 |F'|) ell^o / tau  via the
        # gap-2 kappa-derivative numerators (real-rational tail lemma).
        dA1k = slot0(cm["dA1"])
        NW1 = kp_sub(kp_mul(kp_dk(remW), dA1k), kp_mul(remW, kp_dk(dA1k)))
        dA1p = [x[1].trim() for x in cm["dA1"]]
        NWp = kp_sub(kp_mul(remWp, dA1k), kp_mul(remW, dA1p))
        NWp2 = kp_sub(kp_mul(kp_dk(NWp), dA1k),
                      [2 * x for x in kp_mul(NWp, kp_dk(dA1k))])
        kb = sup_ratio_boxes(dict(KW=(NW1, (2, 2, 2)), KWp=(NWp2, (3, 3, 3))),
                             pr, uq, tstar, lo, hi, thU=thU, dp_raise=dpr,
                             net_offs=dict(KW=noW, KWp=noWp))
        DWKI, DWPK = kb["KW"][2], kb["KWp"][2]
        if not (arb(DWKI) < arb(float("inf"))
                and arb(DWPK) < arb(float("inf"))):
            raise ArithmeticError("Im tail bound not finite")
    LL = (arb(bnds["LT"][0]) + arb(bnds["LU"][0]) + arb(bnds["LB"][0])).upper()
    tt_b, ttp_b = lc_abs(tt, lo, hi), lc_abs(ttp, lo, hi)
    out = dict(Rbar=float(Rbar), DW=float(DW), DWdec=float(DWdec),
               DWp=float(DWp), DB=float(DB), DBdec=float(DBdec), LL=float(LL),
               thU=float(thU), cells=[], ok=True,
               diag={k: sorted(v, key=lambda d: -d[3])[:2] for k, v in diag.items()})
    cl_ = cl
    closing = closing0
    for e0, e1 in eps_cells:
        eb = arb((e0 + e1) / 2, (e1 - e0) / 2)
        rhoR = w0_ - (eb + c) * tt                     # explicit part of rho
        G = ttp + rhoR * tt
        K2 = -b2_ - (1 / c) * ((2 * eb) * b2_ + b1_) - (1 / (2 * c)) * G
        K0 = (eb * eb) * b2_ + eb * b1_ + b0_ - (c / 2) * G
        rho_b = lc_abs(w0_, lo, hi) + tt_b.scale(arb(EPS_MAX) + cl_)
        # slacks of the tau-quadratic (nu = Im(remW/denA1) etc.; all constants
        # are tau-uniform sups, the "dec" variants carry one 1/tau decay):
        #   K += -(1/2)[ d(rho) f + rho d(f) + d(f') ] + nu^2/4  where
        #   d(f) = -(tau/c) nu, d(f') = -(tau/c) nu',
        #   nu' = (remW/denA1)' bounded by DWp + DW*LL.
        R0 = (Bnd.of(arb(DB) + (0 if closing else arb(DBdec) / cl_),
                     noB)                                  # (1+tau/|c|) remB1
              + Bnd.of(arb(DW) ** 2 / 4, 2 * noW)          # nu^2/4 Young
              + tt_b.scale(cl_ * arb(DW) / 2, noW))        # (c/2) tt d-rho
        if not closing:
            R0 = R0 + Bnd.of(arb(DW) * arb(DWdec) / (2 * cl_ * arb(tstar)),
                             2 * noW)                      # d-rho d-f
        if closing:
            # nu-rem and nu'-rem via the Im tail bounds: |nu-rem| <= DWKI
            # ell^noW / tau, |nu'-rem| <= DWPK ell^noWp / tau; the tau-linear
            # and tau^2 slack consumers then land in R0/R1 with structural
            # Laurent offsets
            R0 = R0 + Bnd.of(arb(DWPK) / (2 * cl_), noWp) \
                + Bnd.of(arb(DW) * arb(DWKI) / (2 * cl_), 2 * noW)
            R0 = R0 + rho_b.scale(arb(DWKI) / (2 * cl_), noW)
            # (tau/|c|) |Im(G/dB1)| <= (tau/|c|) DG ell^{-2}: tau-linear slot
            R1 = tt_b.scale(arb(DWKI) / (2 * cl_), noW) \
                + Bnd.of(arb(DB) / cl_, noB)
        else:
            R1 = (Bnd.of(arb(DWp) / (2 * cl_), noWp)
                  + Bnd.of(arb(DW) * arb(LL) / (2 * cl_), noW + LLo)  # nu'
                  + rho_b.scale(arb(DW) / (2 * cl_), noW)  # (tau/2c) rho_R nu
                  + tt_b.scale(arb(DWdec) / (2 * cl_), noW))  # d-rho f tau^2
        okc, m = tau_quad_check(K2, R1, K0, R0, tstar, lo, hi, scale_off)
        try:
            okd, wd, dfails = direct_S_boxes(pr, cm, uqj, remUj, lo, hi, b1,
                                             tstar, [(e0, e1)], scale_off)
        except ArithmeticError as e:
            okd, wd, dfails = False, None, []
            m["direct_err"] = str(e)
        m["direct"] = wd
        if dfails:
            m["direct_fails"] = dfails
        okc = okc and okd
        out["cells"].append(dict(eps=[e0, e1], ok=okc, **m))
        out["ok"] = out["ok"] and okc
    return out


def tau_quad_check(K2, R1, K0, R0, b1, lo, hi, scale_off):
    """sup_{tau >= b1} [K2 tau^2 + R1 tau + K0 + R0] < 0 pointwise on [lo,hi]
    (all quantities multiplied by ell^scale_off).  If lo = 0 and K2's leading
    behaviour vanishes there (positive offset), uses the hat-split
    Q <= ell^h [K2^ tau^2 + R1^ tau] + (K0 + R0)."""
    b1a = arb(b1)
    K2, K0 = K2.trim(), K0.trim()
    A = K2.tm.eval(lo, hi)                  # hat of K2 at its own offset
    hK2 = K2.off + scale_off
    hR1 = R1.min_off() + scale_off
    C = (LC(K0.tm, K0.off + scale_off).eval(lo, hi) if K0.off + scale_off >= 0
         or lo > 0 else None)
    if C is None:
        return False, dict(err="K0 offset negative at 0")
    Cup = C.upper() + R0.eval_up(lo, hi, scale_off)
    m = dict(K2hat=float(A.upper()), K0R0=float(Cup))
    if not (A < 0):
        return False, dict(err="K2 hat not negative", **m)
    if hK2 < 0 and lo <= 0:
        return False, dict(err="K2 offset negative at 0", **m)
    h = min(hK2, hR1)
    if h < 0 and lo <= 0:
        return False, dict(err="R1 offset negative at 0", **m)
    # pointwise ell^h [K2h tau^2 + R1h tau] with K2h <= A*ellpow, R1h <= B
    el = arb((lo + hi) / 2, (hi - lo) / 2)
    Ah = (A * abs(el) ** (hK2 - h)).upper()
    Bh = arb(0)
    for v, o in R1.t:
        Bh += arb(v) * abs(el) ** (o + scale_off - h)
    Bh = Bh.upper()
    Aup = arb(Ah)
    if not (Aup < 0):
        # ell^{hK2-h} hull contains 0: fall back to requiring B-part tiny vs C
        return False, dict(err="K2 hat lost sign in hull", **m)
    tstar = arb(Bh) / (2 * abs(Aup))
    Qb1 = Aup * b1a ** 2 + arb(Bh) * b1a
    sup_in = arb(Bh) ** 2 / (4 * abs(Aup))
    sup_quad = Qb1.upper() if bool(tstar <= b1a) else max(Qb1.upper(), sup_in.upper())
    elh = abs(el) ** h
    tot = arb(sup_quad).max(arb(0)) * elh.upper() + Cup if h > 0 else \
        arb(sup_quad) * 1 + Cup
    if h > 0:
        tot = (arb(sup_quad).max(arb(0)) * elh).upper() + Cup
    else:
        tot = arb(sup_quad).upper() + Cup
    m["margin"] = float(tot if isinstance(tot, float) else arb(tot).upper())
    return bool(arb(tot) < 0), m


def _centre_expansion(prec=256):
    from . import centre
    root = os.path.join(os.path.dirname(__file__), "..", "..", "..")
    aj = json.load(open(os.path.join(root, "results", "a3_midpoint.json")))
    ball = lambda s: arb(s[1:-1].split(" +/-")[0]) + arb(0, float(s[1:-1].split(" +/- ")[1]))
    with precision(prec):
        a, mu = ball(aj["a"]), ball(aj["mu"])
        ce = centre.centre_expansion(mu * (2 * a).exp(), nhat=(-a).exp(), K=30)
        cert = ce.certify()
        if not cert.ok:
            raise RuntimeError("centre certificate failed")
    return ce


def eps_cells_default(n=None):
    e = [0.0, 0.03, 0.07, 0.12, 0.2, 0.35, 0.6, 1.0, 1.6, 2.4, 3.4, 4.6, 6.0]
    return list(zip(e, e[1:]))


def y_subranges(y0=1e-8, ratio=1.6):
    out, hi = [], YCUT
    while hi > y0 * ratio:
        lo = hi / ratio
        out.append((lo, hi))
        hi = lo
    out.append((0.0, hi))
    return out[::-1]


def centre_piece(at, ce, lo, hi, b1, eps_cells, prec):
    """Build and check one centre piece on its own TM domain [0, hi]."""
    with precision(prec):
        comps, _ = centre_components(ce, hi)
        ak = AtomKP(at, comps, sd_shift=0)
        pr = build_pairs(ak)
        cm = commute(pr, pr["dT"][-1][0].trim().inv())
        return analyze_piece(pr, cm, lo, hi, b1, eps_cells, 0, (0, 0, 0))


def run_centre(at, b1, prec=320, y0=1e-8, ratio=1.6, ncells=6):
    """Strip-exclusion conditions on the centre region y = e^x in (0, e^{-4.5}]:
    per-piece TM domains, geometric subranges down to y0 plus the closing piece
    [0, y0]."""
    t0 = time.time()
    ce = _centre_expansion()
    pieces, ok, worst = [], True, None
    for lo, hi in y_subranges(y0, ratio):
        try:
            rep = centre_piece(at, ce, lo, hi, b1, eps_cells_default(ncells), prec)
        except (ArithmeticError, RuntimeError) as e:
            rep = dict(ok=False, err=str(e), cells=[])
        ok = ok and rep["ok"]
        mg = max((c.get("margin") for c in rep.get("cells", [])
                  if c.get("margin") is not None), default=None)
        pieces.append(dict(lo=lo, hi=hi, ok=rep["ok"], margin=mg,
                           err=rep.get("err"),
                           errs=[c.get("err") for c in rep.get("cells", [])
                                 if c.get("err")][:2] or None))
        if mg is not None:
            worst = mg if worst is None else max(worst, mg)
    return dict(ok=ok, pieces=len(pieces), worst_margin=worst,
                time=time.time() - t0, detail=pieces)


def run_tube(at, tube, b1, x_min=-4.5, x_max=-0.05, prec=320, step_slice=None,
             degrade=False):
    """Strip-exclusion conditions on tube steps covering [x_min, x_max]:
    per-step TM domain sigma = x_k - x in [0, h].  degrade=True: if the
    tau > tstar quadratic phase fails on a step, certify the direct phase
    (tau in [b1, tstar]) and mark the tail open (near-sonic band)."""
    from .tc_tm import step_components
    t0 = time.time()
    steps = [sd for sd in tube.steps
             if float(sd.x) - sd.h >= x_min - 1e-12 and float(sd.x) <= x_max + 1e-12]
    if step_slice is not None:
        steps = steps[step_slice[0]:step_slice[1]]
    pieces, ok, worst = [], True, None
    with precision(prec):
        for sd in steps:
            try:
                comps = step_components(sd)
                ak = AtomKP(at, comps, sd_shift=0)
                pr = build_pairs(ak)
                cm = commute(pr, pr["dT"][-1][0].trim().inv())
                try:
                    rep = analyze_piece(pr, cm, 0.0, float(sd.h), b1,
                                        eps_cells_default(), 0, (0, 0, 0))
                except (ArithmeticError, RuntimeError):
                    if not degrade:
                        raise
                    rep = dict(ok=False)
                if degrade and not rep["ok"]:
                    rep = analyze_piece(pr, cm, 0.0, float(sd.h), b1,
                                        eps_cells_default(), 0, (0, 0, 0),
                                        tq=False)
            except (ArithmeticError, RuntimeError) as e:
                rep = dict(ok=False, err=str(e), cells=[])
            ok = ok and rep["ok"]
            mg = max((c.get("margin") for c in rep.get("cells", [])
                      if c.get("margin") is not None), default=None)
            if os.environ.get("P4_PROGRESS"):
                print("progress tube step x=%.5f ok=%s margin=%s t=%.0fs"
                      % (float(sd.x), rep["ok"], mg, time.time() - t0),
                      file=sys.stderr, flush=True)
            pieces.append(dict(x=float(sd.x), h=float(sd.h), ok=rep["ok"],
                               margin=mg, err=rep.get("err"),
                               open=rep.get("tau_quad"),
                               errs=[c.get("err") for c in rep.get("cells", [])
                                     if c.get("err")][:2] or None))
            if mg is not None:
                worst = mg if worst is None else max(worst, mg)
    return dict(ok=ok, pieces=len(pieces), worst_margin=worst,
                time=time.time() - t0, detail=pieces)


def xi_subranges(xi0=1e-6, ratio=1.2):
    """[0, xi0] Laurent closing piece + geometric recentred pieces to 0.05
    (capped at absolute width 2e-3 near the outer end)."""
    out, lo = [(0.0, xi0)], xi0
    while lo < XI_MAX - 1e-12:
        hi = min(lo * ratio, lo + 2e-3, XI_MAX)
        out.append((lo, hi))
        lo = hi
    return out


def run_sonic(at, b1, prec=320, xi0=1e-6, piece_slice=None):
    """Strip-exclusion conditions on the sonic region xi = -x in [0, 0.05]:
    pieces with lo >= xi0 are recentred to their own domain (no vanishing
    factors there); the closing piece [0, xi0] keeps the Laurent machinery
    with kdiv leads offset-raised by (1, 2, 3) Dp-powers (identities E1/E3)
    and all conditions xi^2-scaled."""
    from .tc_energy import sonic_expansion_cached
    from .tc_tm import sonic_components, sonic_components_recentred
    t0 = time.time()
    ex = sonic_expansion_cached(prec)
    pieces, ok, worst = [], True, None
    subs = [t for t in xi_subranges(xi0) if t[0] < XI_SONIC_TOP]
    if piece_slice is not None:
        subs = subs[piece_slice[0]:piece_slice[1]]
    with precision(prec):
        for lo, hi in subs:
            try:
                if lo >= xi0 * 0.99 and hi > XI_LAUR:
                    # wide-xi pieces: recentred quasi-step series
                    comps = sonic_components_recentred(ex, lo, hi)
                    ak = AtomKP(at, comps, sd_shift=0)
                    pr = build_pairs(ak)
                    cm = commute(pr, pr["dT"][-1][0].trim().inv())
                    rep = analyze_piece(pr, cm, 0.0, hi - lo, b1,
                                        eps_cells_default(), 0, (0, 0, 0))
                elif lo >= xi0 * 0.99:
                    # small-xi pieces: Laurent machinery on the piece's own
                    # domain [0, hi], evaluated on [lo, hi] (lo > 0: all
                    # Laurent offset checks are finite).  If the tau > tstar
                    # quadratic bounds fail, degrade to the direct phase
                    # (tau in [b1, tstar] certified; the tail marked open).
                    comps = sonic_components(ex, hi)
                    ak = AtomKP(at, comps, sd_shift=1)
                    pr = sonic_raise_inputs(build_pairs(ak))
                    cm = commute(pr, pr["dT"][-1][0].trim().inv())
                    try:
                        rep = analyze_piece(pr, cm, lo, hi, b1,
                                            eps_cells_default(), 2, (0, 0, 0),
                                            dp_raise=1)
                    except (ArithmeticError, RuntimeError):
                        rep = dict(ok=False)
                    if not rep["ok"]:
                        rep = analyze_piece(pr, cm, lo, hi, b1,
                                            eps_cells_default(), 2, (0, 0, 0),
                                            dp_raise=1, tq=False)
                else:
                    # the closing piece [0, xi0]: the direct phase (tau in
                    # [b1, tstar]) is certified with ell = 0-safe offset-aware
                    # chains; the tau > tstar tail on this piece is an OPEN
                    # item (the S ell^2 / tau^2 = F(ell tau) corner needs a
                    # joint (w, ell) engine; float truth: F <= -0.9)
                    comps = sonic_components(ex, hi)
                    ak = AtomKP(at, comps, sd_shift=1)
                    pr = sonic_raise_inputs(build_pairs(ak))
                    cm = commute(pr, pr["dT"][-1][0].trim().inv())
                    rep = analyze_piece(pr, cm, lo, hi, b1,
                                        eps_cells_default(), 2, (0, 0, 0),
                                        dp_raise=1, tq=False)
            except (ArithmeticError, RuntimeError) as e:
                rep = dict(ok=False, err=str(e), cells=[])
            ok = ok and rep["ok"]
            mg = max((c.get("margin") for c in rep.get("cells", [])
                      if c.get("margin") is not None), default=None)
            pieces.append(dict(lo=lo, hi=hi, ok=rep["ok"], margin=mg,
                               err=rep.get("err"),
                               open=rep.get("tau_quad"),
                               errs=[c.get("err") for c in rep.get("cells", [])
                                     if c.get("err")][:2] or None))
            if mg is not None:
                worst = mg if worst is None else max(worst, mg)
            if os.environ.get("P4_PROGRESS"):
                print("progress sonic piece %g ok=%s margin=%s t=%.0fs"
                      % (lo, rep["ok"], mg, time.time() - t0),
                      file=sys.stderr, flush=True)
    return dict(ok=ok, pieces=len(pieces), worst_margin=worst,
                time=time.time() - t0, detail=pieces)


TSTAR_FAC = 16.0                 # direct boxes cover tau in [b1, TSTAR_FAC*b1]
XI_LAUR = 2.1e-2                 # sonic pieces below this use Laurent form
XI_SONIC_TOP = 8.5e-3            # sonic region top: the gap tube covers
                                 # xi in [0.0085, 0.055] (overlapping piece 50)


def direct_S_boxes(pr, cm, uqj, remUj, lo, hi, b1, tstar, eps_cells, scale_off,
                   nstrip=5, ntau=20):
    """Direct ball evaluation of S = K + nu^2/4 on (ell-strip) x (eps, tau)
    boxes covering tau in [b1, tstar], mirroring the commutation recursion on
    the small value polynomials (nT/dT depth 2, nB/dB depth 2, tt jet):
        A0 = nT/dT, B0 = nB/dB, lB0 = B0'/B0, A1 = A0 - lB0,
        A1' = A0' - (B0''/B0 - lB0^2), B1 = B0 + A0' - A0 lB0,
        nu = Im A1 + tau tt, nu' = Im A1' + tau tt', rho = Re A1 - c tt,
        f = (c + tau^2/c) tt - (tau/c) nu, K = Re B1 - (tau/c) Im B1
            - (f' + rho f)/2,  S = K + nu^2/4.
    Adaptive in both the ell-strips and the (eps, tau) boxes."""
    from flint import acb
    c = cval()
    worst, ok = None, True
    fails = []
    kps = {}
    for nm, kp in (("nT", pr["nT"]), ("dT", pr["dT"]), ("nB", pr["nB"]),
                   ("dB", pr["dB"]), ("uq", uqj), ("rU", remUj)):
        for sl in (0, 1, 2):
            kps[nm + str(sl)] = [x[sl] for x in kp]
    kps["tt"] = [cm["tt"][0]]
    kps["ttp"] = [cm["tt"][1]]
    # tau edges: fine near b1 (the w-Horner hulls need small relative tau
    # widths where the margins are thin), then geometric to tstar
    edges = [b1 * r for r in (1.0, 1.005, 1.01, 1.02, 1.04, 1.08, 1.15)]
    while edges[-1] < tstar:
        edges.append(min(edges[-1] * 1.15, tstar))

    def do_strip(slo, shi, sdepth):
        nonlocal worst, ok
        # offset-aware hulls: each chain value is (complex hull, ell-offset)
        # so the ell = 0 strip works (S ell^k < 0 iff its hat hull < 0)
        hulls, mo = {}, {}
        ehi, elo_l = arb(max(shi, 0.0)), arb(max(slo, 0.0))

        def elp(d):
            if d == 0:
                return arb(1)
            return (elo_l ** d).union(ehi ** d)
        for n, kp in kps.items():
            hh = []
            for x in kp:
                x = x.trim()
                hh.append((x.tm.eval(slo, shi), x.off))
            hulls[n] = hh
            mo[n] = min(o for _, o in hh)
        degs = dict(nT=6, dT=5, nB=9, dB=7)

        def omul(a, b):
            return (a[0] * b[0], a[1] + b[1])

        def odiv(a, b):
            return (a[0] / b[0], a[1] - b[1])

        def oadd(*ts):
            m = min(t[1] for t in ts)
            v = None
            for t in ts:
                x = t[0] * elp(t[1] - m)
                v = x if v is None else v + x
            return (v, m)

        def oneg(a):
            return (-a[0], a[1])

        def hat(n, wv, deg):
            acc, first = acb(0), True
            m = mo[n]
            for j in range(deg + 1):
                hh = hulls[n]
                x = hh[j][0] * elp(hh[j][1] - m) if j < len(hh) else arb(0)
                acc = x + 0 * wv if first else acc * wv + x
                first = False
            return (acc, m)

        def sbox(elo_, ehi_, t0, t1):
            eb = arb((elo_ + ehi_) / 2, (ehi_ - elo_) / 2)
            tb = arb((t0 + t1) / 2, (t1 - t0) / 2)
            z = acb(eb, tb)
            wv = acb(1) / z
            v = {n + str(sl): hat(n + str(sl), wv, degs[n])
                 for n in ("nT", "dT", "dB") for sl in (0, 1, 2)}
            for dn in ("dT0", "dB0"):
                if not (abs(v[dn][0]).lower() > 0):
                    return None
            zp = (z, 0)
            A0 = omul(zp, odiv(v["nT0"], v["dT0"]))
            A0p = omul(zp, odiv(oadd(omul(v["nT1"], v["dT0"]),
                                     oneg(omul(v["nT0"], v["dT1"]))),
                                omul(v["dT0"], v["dT0"])))
            quadj = []
            for sl in (0, 1, 2):
                hh = hulls["uq" + str(sl)]
                ts = [omul((hh[j][0], hh[j][1]), ((z ** j, 0) if j else
                                                  (acb(1), 0)))
                      for j in range(len(hh))]
                quadj.append(oadd(*ts))
            rUj = [hat("rU" + str(sl), wv, 6) for sl in (0, 1, 2)]
            dBv, dBp, dBpp = v["dB0"], v["dB1"], v["dB2"]
            wp = (wv, 0)
            Uv = oadd(quadj[0], omul(odiv(rUj[0], dBv), wp))
            Up = oadd(quadj[1],
                      omul(odiv(oadd(omul(rUj[1], dBv),
                                     oneg(omul(rUj[0], dBp))),
                                omul(dBv, dBv)), wp))
            Upp = oadd(quadj[2],
                       omul(odiv(oadd(omul(omul(rUj[2], dBv), dBv),
                                      oneg((2 * (rUj[1][0] * dBv[0] * dBp[0]),
                                            rUj[1][1] + dBv[1] + dBp[1])),
                                      (rUj[0][0] * (2 * dBp[0] * dBp[0]
                                                    - dBv[0] * dBpp[0]),
                                       rUj[0][1] + 2 * dBp[1])),
                                 omul(omul(dBv, dBv), dBv)), wp))
            if not (abs(Uv[0]).lower() > 0):
                return None
            lB0 = odiv(Up, Uv)
            lB0p = oadd(odiv(Upp, Uv), oneg(omul(lB0, lB0)))
            A1 = oadd(A0, oneg(lB0))
            A1p = oadd(A0p, oneg(lB0p))
            B1 = oadd(Uv, A0p, oneg(omul(A0, lB0)))
            tt = hat("tt", wv, 0)
            ttp = hat("ttp", wv, 0)
            ttr, ttpr = (tt[0].real, tt[1]), (ttp[0].real, ttp[1])
            tbp = (tb, 0)
            nu = oadd((A1[0].imag, A1[1]), omul(tbp, ttr))
            nup = oadd((A1p[0].imag, A1p[1]), omul(tbp, ttpr))
            rho = oadd((A1[0].real, A1[1]), (-c * ttr[0], ttr[1]))
            cf = ((c + tb * tb / c), 0)
            f = oadd(omul(cf, ttr), ((-tb / c) * nu[0], nu[1]))
            fp = oadd(omul(cf, ttpr), ((-tb / c) * nup[0], nup[1]))
            K = oadd((B1[0].real, B1[1]), ((-tb / c) * B1[0].imag, B1[1]),
                     (arb(-1) / 2 * fp[0], fp[1]),
                     (arb(-1) / 2 * (rho[0] * f[0]), rho[1] + f[1]))
            S = oadd(K, (nu[0] * nu[0] / 4, 2 * nu[1]))
            # S = hat * ell^{S[1]}: for ell in the OPEN (0, shi] the sign of
            # ell^k S-hat is the sign of the hat; test/report the hat hull
            return S[0]

        bad = []

        def check_box(elo_, ehi_, t0, t1, depth=0):
            S = sbox(elo_, ehi_, t0, t1)
            if S is not None and (S < 0):
                return S.upper()
            if depth >= 16:
                bad.append((elo_, ehi_, t0, t1,
                            float(S.upper()) if S is not None else None))
                return S.upper() if S is not None else None
            # split the direction that dominates the hull widths
            w_ = None
            if (ehi_ - elo_) > 0.008 and (ehi_ - elo_) >= 0.5 * (t1 - t0) / t0:
                em = (elo_ + ehi_) / 2
                halves = (((elo_, em), (t0, t1)), ((em, ehi_), (t0, t1)))
            elif (t1 - t0) / t0 > 5e-4:
                tm_ = (t0 * t1) ** 0.5
                halves = (((elo_, ehi_), (t0, tm_)), ((elo_, ehi_), (tm_, t1)))
            elif (ehi_ - elo_) > 0.001:
                em = (elo_ + ehi_) / 2
                halves = (((elo_, em), (t0, t1)), ((em, ehi_), (t0, t1)))
            else:
                bad.append((elo_, ehi_, t0, t1,
                            float(S.upper()) if S is not None else None))
                return S.upper() if S is not None else None
            for ee, ttt in halves:
                u = check_box(ee[0], ee[1], ttt[0], ttt[1], depth + 1)
                if u is not None:
                    w_ = u if w_ is None else max(w_, u)
            return w_

        wloc = None
        for e0, e1 in eps_cells:
            for t0, t1 in zip(edges, edges[1:]):
                u = check_box(e0, e1, t0, t1)
                if u is not None:
                    wloc = u if wloc is None else max(wloc, u)
        if bad:
            # the chain-hull slop scales with the ell-strip width (atom
            # variation amplified through the nu-cancellation), so deep strip
            # refinement -- not (eps, tau) splitting -- is what converges
            if sdepth < 9:
                m_ = (slo + shi) / 2 if slo > 0 else shi / 2
                do_strip(slo, m_, sdepth + 1)
                do_strip(m_, shi, sdepth + 1)
                return
            ok = False
            fails.extend(bad[:2])
        if wloc is not None:
            worst = wloc if worst is None else max(worst, wloc)

    for slo, shi in strip_list(lo, hi, nstrip):
        do_strip(slo, shi, 0)
    return ok, (float(worst) if worst is not None else None), fails[:4]




# ---------------------------------------------------------- Task 1: Lemma D
def run_lemma_d(at, c0=C0_LEMMA_D, prec=320, ycut=YCUT, nsub=1100):
    """Lemma D closure: tc_energy's conditions (a)-(d) certified on the whole
    deep tail x <= -4.5 (y = e^x in (0, e^{-4.5}]) in one Taylor-model
    computation.  The junk-free evaluation rests on the exact layer identity E2
    (check_layer_identity): Z := QNr - (9/4) Dp^3 SDr^3 D2r^2 Bcr^2 has no
    t-layers below 16 with A -> 1 and no odd layers at all, while A - 1 =
    O(t^2) exactly; hence val(Z o u) >= 16, val(QD o u) = 14, and the division
    Q - 9/4 = Z/QD carries the honest t^2 relative scale:
        q2 > 0,  2 c0 q2 + q1 - B1 > 0,  36 q2 + 6 q1 + q0 - B0 > 0
    on every y-subrange, poles via Routh-Hurwitz on the Bcr hat cubic and the
    real A(x), -F_N(x) balls (condition (d)).  B0/B1 are tau-boxed line sups
    on Re kappa = c0 (tc_energy.line_sup) with factor-wise denominators."""
    from .tc_energy import RBAR, kdiv, line_sup
    t0 = time.time()
    pmin, okE2 = check_layer_identity(at)
    if not okE2:
        raise RuntimeError("layer identity E2 failed")
    ce = _centre_expansion()
    out = dict(ok=True, c0=c0, ycut=ycut, layer_identity=[pmin, okE2],
               subranges=[], mins={})
    with precision(prec):
        comps, _ = centre_components(ce, ycut)
        ak = AtomKP(at, comps, sd_shift=0)
        # QN via the compiled exact QNr atom (identical object to E2's Z + QQD)
        QN = ak.value_kp("QNr")
        Dp = ak.value_kp("Dp")
        SDr = ak.value_kp("SDr")
        D2r = ak.value_kp("D2r")
        Bcr = ak.value_kp("Bcr")
        p3 = kpj_mul([[x] for x in Dp], [[x] for x in Dp], 0)
        p3 = kpj_mul(p3, [[x] for x in Dp], 0)
        s3 = kpj_mul([[x] for x in SDr], [[x] for x in SDr], 0)
        s3 = kpj_mul(s3, [[x] for x in SDr], 0)
        d2 = kpj_mul([[x] for x in D2r], [[x] for x in D2r], 0)
        b2 = kpj_mul([[x] for x in Bcr], [[x] for x in Bcr], 0)
        QD = [x[0] for x in kpj_mul(kpj_mul(p3, s3, 0), kpj_mul(d2, b2, 0), 0)]
        # Z = QN - 9/4 QD, exact-shift by y^16 (justified by E2 + A-1 = O(y^2))
        QDp = QD + [lc_zero(ycut)] * (len(QN) - len(QD))
        Z = [(a_ - (arb(9) / 4) * b_).trim() for a_, b_ in zip(QN, QDp)]
        Zs = []
        for x in Z:
            need = 16 - x.off
            Zs.append(lc_raise(x, need) if need > 0 else x)
        QDt = [x.trim() for x in QD]
        mins = {}
        # dense geometric subranges over the outer three decades (the hats
        # vary there); below ycut/1000 one closing subrange suffices (the
        # Taylor models are limit-dominated)
        edges = [ycut * 1.008 ** (-k) for k in range(nsub)][::-1]
        for i in range(len(edges)):
            lo = 0.0 if i == 0 else edges[i - 1]
            hi = edges[i]
            # hat division: Q - 9/4 = y^2 * (Zs/y^14-part) / QD ... work with
            # plain hull coefficient lists (both sides carry y^14; the shifted
            # numerator carries y^16, so the quotient q2,q1,q0-9/4 and the
            # remainder are honestly y^2-scaled)
            num = [x.trim() for x in Zs]
            nch = []
            moffQD = min(x.off for x in QDt)
            for x in num:
                d = x.off - (moffQD + 2)
                if d < 0:
                    raise ArithmeticError("Z coefficient below y^2 over QD")
                pw = arb(1) if d == 0 else \
                    (arb(max(lo, 0.0)) ** d).union(arb(hi) ** d)
                nch.append(x.tm.eval(lo, hi) * pw)
            dch = []
            for x in QDt:
                d = x.off - moffQD
                if d < 0:
                    raise ArithmeticError("QD offset below its minimum")
                pw = arb(1) if d == 0 else \
                    (arb(max(lo, 0.0)) ** d).union(arb(hi) ** d)
                dch.append(x.tm.eval(lo, hi) * pw)
            q2_, q1_, q0s, rem = kdiv(nch, dch)
            q0_ = arb(9) / 4 + q0s
            # dens for line_sup: hat factors (Dp^3 SDr^3) kappa-free, D2r^2, Bcr^2
            def hulled(kp, m):
                o = [None] * len(kp)
                for j, x in enumerate(kp):
                    x = x.trim()
                    d = x.off - m
                    pw = arb(1) if d == 0 else \
                        (arb(max(lo, 0.0)) ** d).union(arb(hi) ** d)
                    o[j] = x.tm.eval(lo, hi) * pw
                return o
            mD = min(x.off for x in Dp) * 3 + min(x.off for x in SDr) * 3
            free = [hulled([x[0] for x in kpj_mul(p3, s3, 0)], mD)[0]]
            d2h = hulled(D2r, min(x.trim().off for x in D2r))
            bch = hulled(Bcr, min(x.trim().off for x in Bcr))
            B0b, B1b = line_sup(rem, dch, c0, dens=[free, d2h, d2h, bch, bch])
            A_, W_ = ak.A.eval(max(lo, 1e-300), hi), ak.W.eval(lo, hi)
            FN_ = A_ - 2 - 2 * W_ / 3
            from .tc_energy import hurwitz_shift as hs
            m = dict(q2=q2_, mB=2 * c0 * q2_ + q1_ - B1b,
                     mC=q2_ * c0 ** 2 + q1_ * c0 + q0_ - B0b,
                     poleA=arb(RBAR) - A_, poleFN=arb(RBAR) + FN_)
            if not hs(bch, RBAR):
                raise ArithmeticError("Hurwitz on the Bcr hat cubic failed")
            okp = all(bool(v > 0) for v in m.values())
            out["ok"] = out["ok"] and okp
            for k, x in m.items():
                mins[k] = x if k not in mins else mins[k].min(x)
            out["subranges"].append(dict(lo=lo, hi=hi, ok=okp))
        out["mins"] = {k: float(v.lower()) for k, v in mins.items()}
    out["time"] = time.time() - t0
    return out


# ------------------------------------------------------------------------ main
def _out_path(name="tc_commutator.json"):
    return os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "results", "theorem_c", name)


def main(argv):
    """Modes:
      checks                       exact identities E1/E2/E3 (~30 s)
      lemma_d                      Task 1: deep-tail closure (writes into
                                   tc_energy.json via tc_energy centre mode)
      centre | sonic [a b] | tube <tube.json> [a b]
                                   region runs (optional piece/step slices,
                                   for worker parallelism)
      merge <part.json> ...        combine region/worker reports and the checks
                                   into results/theorem_c/tc_commutator.json
    The certified statement and its proof are in notes/s4-commutator.md."""
    mode = argv[0] if argv else "checks"
    t0 = time.time()
    at = Atoms(full=True)
    if mode == "checks":
        res = dict(E1_U6top=check_u6_top_identity(at),
                   E2_layers=check_layer_identity(at),
                   E3_sonic_zeros=check_sonic_zeros(at),
                   E5_sonic_vanishing=check_sonic_vanishing(at))
        print(res)
        return res
    if mode == "lemma_d":
        res = run_lemma_d(at)
        print("lemma_d ok=%s mins=%s" % (res["ok"], res["mins"]))
        return res
    if mode == "centre":
        res = run_centre(at, B1_DEFAULT)
    elif mode == "sonic":
        sl = (int(argv[1]), int(argv[2])) if len(argv) > 2 else None
        res = run_sonic(at, B1_DEFAULT, piece_slice=sl)
    elif mode == "tube":
        from .lintube import Tube
        tube = Tube.load(argv[1])
        sl = (int(argv[2]), int(argv[3])) if len(argv) > 3 else None
        res = run_tube(at, tube, B1_DEFAULT, step_slice=sl)
    elif mode == "merge":
        parts = [json.load(open(p)) for p in argv[1:]]
        det = []
        ok = True
        worst = None
        for p in parts:
            ok = ok and p["ok"]
            det.extend(p.get("detail", []))
            if p.get("worst_margin") is not None:
                w = float(p["worst_margin"])
                worst = w if worst is None else max(worst, w)
        opens = [d for d in det if d.get("open")]
        # xi-extent of the open (tau > tstar) pieces: sonic pieces carry
        # (lo, hi) in xi; tube steps carry (x, h) with xi = -x + sigma
        xi_open = max([d["hi"] for d in opens if "hi" in d]
                      + [-d["x"] + d["h"] for d in opens if "x" in d],
                      default=0.0)
        scope = dict(
            full_strip_tau=[B1_DEFAULT, TSTAR_FAC * B1_DEFAULT],
            xi_open_max=xi_open,
            tau_above_tstar=("certified everywhere" if not opens else
                             "certified for xi > %g; OPEN on the near-sonic "
                             "band xi <= %g (see open_note)"
                             % (xi_open, xi_open)))
        res = dict(ok=ok, scope=scope, b1=B1_DEFAULT, m=MC, c=[CNUM, CDEN],
                   eps_max=EPS_MAX, tstar_fac=TSTAR_FAC,
                   open_pieces=len(opens),
                   open_note=(("on %d near-sonic pieces (xi <= %g) only tau "
                               "in [b1, %g] is certified (direct phase); the "
                               "tau > tstar tail there is an open item -- "
                               "see notes/s4-commutator.md") %
                              (len(opens), xi_open, TSTAR_FAC * B1_DEFAULT))
                   if opens else None,
                   checks=dict(E1_U6top=check_u6_top_identity(at),
                               E2_layers=check_layer_identity(at),
                               E3_sonic_zeros=check_sonic_zeros(at),
                               E5_sonic_vanishing=check_sonic_vanishing(at)),
                   pieces=len(det), worst_margin=worst, detail=det,
                   time=time.time() - t0)
        json.dump(res, open(_out_path(), "w"), indent=1, default=str)
        print("tc_commutator.json written: ok=%s pieces=%d worst=%s"
              % (ok, len(det), worst))
        return res
    else:
        raise SystemExit("unknown mode " + mode)
    print("%s: ok=%s pieces=%s worst=%s time=%.0fs"
          % (mode, res["ok"], res.get("pieces"), res.get("worst_margin"),
             res["time"]))
    return res


if __name__ == "__main__":
    main(sys.argv[1:])
