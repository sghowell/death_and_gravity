"""S4 regime L (0 <= Re omega <= Omega_L = 8 kappa_-): validated scaled Frobenius series at both horizons,
uniformly in h = 1/(l + 1/2) in (0, h0] (h^2 in [0, h0^2], h = 0 the Bessel limit), then Riccati
invariance through the forbidden bulk.  Pre-registration section 3.1, Amendments 2 and 5.

Scaling.  At r_h (sign s = +1 at r_+, -1 at r_c) put r = r_h + s h^2 tau, tau >= 0.  With
Delta(r_h + t) = sum_k d_k t^k (d_0 = 0) write Delta = h^2 D(tau), D = sum_k d_k s^k h^{2(k-1)} tau^k,
and Delta'(r) = E(tau) = sum_k k d_k (s h^2 tau)^(k-1).  The radial equation
Delta^2 R'' + Delta Delta' R' + (w^2 r^4 - L Delta) R = 0  becomes, exactly (L h^2 = 1 - h^2/4),

    D^2 R_tt + s D E R_t + (w^2 r^4 - (1 - h^2/4) D) R = 0,      r = r_h + s h^2 tau,

a polynomial ODE in tau with coefficients polynomial in h^2.  Indicial roots +-i w/(2 kappa_h); the
Jost solution is R = tau^rho sum a_n tau^n with rho = -i w/(2 kappa_h).  D^2 has leading coefficient
d_1^2, the recurrence denominator is d_1^2 n (n + 2 rho) and arbw.frobenius' majorant tail bound
applies verbatim to these coefficient polynomials.

Arithmetic (Amendment 5).  omega ranges over a box (a disc about its centre) and h^2 over one of
K_H equal sub-intervals of [0, h0^2]; the recurrence runs in first-order affine forms in the two
deviations (affine.AF), which keeps the linear dependence exactly (plain acb rectangles wrap).

Exit.  y = Phi_x/Phi (Phi = r R) = h^2 D/r^3 + s (D/r^2) R_tau/R.  For an omega box, tau_a is the
smallest multiple of 1/64 with min h^-2 U(r_a) = D/r^4 >= M := max(4 |w|^2_max, C_W) (C_W from
bulk_PB.json); between the exits Re(V - w^2) >= M - M/4 - M/4 > 0, and Riccati invariance carries
Re y_+ > 0 (resp. Re y_c < 0) to x_P.  Certified: Re y_+(tau_a) > 0 and Re y_c(tau_a^c) < 0.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx, fmpq

from .affine import AF, up
from .arbutil import Timer, git_commit, precision
from .arbw import _geom_tails, _shift, _solve_q
from .m3 import fa
from .params import Point

PREC = 128
H0 = Fraction(2, 41)
H0SQ = H0 * H0
K_H = 16
ROOT = Path(__file__).resolve().parents[2]
CERT_DIR = ROOT / "certificates" / "part1"
RES_DIR = ROOT / "results" / "part1"
N_MAX = 3000
TAU_STEP = Fraction(1, 64)
BOX = Fraction(1, 20)
MIN_BOX = Fraction(1, 80)
DESIGN_L = {"Omega_L_over_km": "8", "Im_range_over_km": ["-11/20", "0"], "box_side_over_km": str(BOX),
            "min_box_side_over_km": str(MIN_BOX), "h2_range": ["0", str(H0SQ)], "h2_subintervals": K_H,
            "tau_grid": str(TAU_STEP), "exit_rule": "D/r^4 >= max(4|w|^2_max, C_W)"}


def _hull(lo_: arb, hi_: arb) -> arb:
    return (lo_ + hi_) / 2 + arb(0, 1) * ((hi_ - lo_) / 2)


# ---------------------------------------------------------------------------------------------
# the scaled equation
# ---------------------------------------------------------------------------------------------
def params_af(pt: Point, wbox, h2_lo: Fraction, h2_hi: Fraction):
    """(w, h2) as affine forms: wbox = (re_lo, re_hi, im_lo, im_hi) in kappa_- units (Fractions)."""
    km = pt.km
    re_lo, re_hi, im_lo, im_hi = (fa(x) for x in wbox)
    cr, ci = ((re_lo + re_hi) / 2 * km), ((im_lo + im_hi) / 2 * km)
    wc = acb(cr.mid(), ci.mid())
    hr, hi_ = (re_hi - re_lo) / 2 * km, (im_hi - im_lo) / 2 * km
    Rw = up((hr * hr + hi_ * hi_).sqrt()) + up(abs(acb(cr, ci) - wc))
    Rw = up(Rw + arb(2) ** -200)
    hc = fa((h2_lo + h2_hi) / 2)
    Rh = up(fa((h2_hi - h2_lo) / 2) + arb(2) ** -200)
    w = AF(wc, 1, 0, 0, Rw, Rh)
    h2 = AF(hc, 0, 1, 0, Rw, Rh)
    return w, h2


def tau_polys(pt: Point, side: int, w: AF, h2: AF):
    """Ascending tau-coefficients (affine forms) of A = D^2, B = s D E, C = w^2 r^4 - (1 - h2/4) D;
    also D's coefficients and r = r0 + r1 tau."""
    rh = pt.rp if side > 0 else pt.rc
    d = _shift(pt.dco, rh)                      # Delta(r_h + t), d[0] = 0 exactly (not used)
    s = side
    one = w.const(1)
    hp = [one, h2, h2 * h2, h2 * h2 * h2]
    Dc = [one * 0] + [hp[k - 1] * (d[k] * s ** k) for k in range(1, 5)]
    Ec = [hp[k] * ((k + 1) * d[k + 1] * s ** k) for k in range(4)]
    A = [one * 0 for _ in range(9)]
    for i in range(1, 5):
        for j in range(1, 5):
            A[i + j] = A[i + j] + Dc[i] * Dc[j]
    B = [one * 0 for _ in range(8)]
    for i in range(1, 5):
        for j in range(4):
            B[i + j] = B[i + j] + (Dc[i] * Ec[j]) * s
    r0 = one * acb(rh)
    r1 = h2 * s
    r2 = [r0 * r0, (r0 * r1) * 2, r1 * r1]
    r4 = [one * 0 for _ in range(5)]
    for i in range(3):
        for j in range(3):
            r4[i + j] = r4[i + j] + r2[i] * r2[j]
    ww = w * w
    C = [ww * r4[k] - (one - h2 * acb(fmpq(1, 4))) * Dc[k] for k in range(5)]
    return A, B, C, Dc, (r0, r1)


def frob_tau(A, B, C, rho: AF, tau: fmpq, min_terms=30):
    """(S, Sd) with R = tau^rho S, R_tau = tau^rho Sd at tau > 0, as affine forms (tail included)."""
    A2 = A[2]
    ta = arb(tau)
    a = [rho.const(1)]
    S = rho.const(1)
    Sd = rho * (1 / acb(ta))
    tk = arb(1)
    small = 0
    eps = arb(2) ** (-(ctx.prec - 12))
    for n in range(1, N_MAX):
        acc = rho.const(0)
        for k in range(max(0, n - 6), n):
            j2, j1, j0 = n - k + 2, n - k + 1, n - k
            term = rho.const(0)
            if j2 < len(A):
                term = term + A[j2] * ((rho + k) * (rho + (k - 1)))
            if j1 < len(B):
                term = term + B[j1] * (rho + k)
            if j0 < len(C):
                term = term + C[j0]
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
        raise ArithmeticError("frob_tau: N too small relative to |rho|")
    A2l = abs(A2.c0).lower() - A2.dev()
    if not (A2l > 0):
        raise ArithmeticError("frob_tau: |A_2| not bounded below")
    g = []
    for j in range(1, 7):
        Aj = A[j + 2].absu() if j + 2 < len(A) else arb(0)
        Bj = B[j + 1].absu() if j + 1 < len(B) else arb(0)
        Cj = C[j].absu() if j < len(C) else arb(0)
        gj = (Aj * (1 + Rr / Nn) ** 2 / (1 - 2 * Rr / Nn) + Bj * (Nn + Rr) / (Nn * (Nn - 2 * Rr))
              + Cj / (Nn * (Nn - 2 * Rr))) / A2l
        g.append(up(gj))
    q = arb(_solve_q(g))
    if not (sum((gj * q ** (-(j + 1)) for j, gj in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("frob_tau majorant root not verified")
    K = max((a[m].absu() / q ** m for m in range(M - 6, M)), key=lambda x: x.upper())
    s0, s1 = _geom_tails(up(K), q, abs(ta), M, Rr)
    S = S + AF(0, 0, 0, s0, S.Rw, S.Rh)
    Sd = Sd + AF(0, 0, 0, s1, S.Rw, S.Rh)
    return S, Sd, {"terms": M, "x": float((q * abs(ta)).upper())}


def scaled_U(pt: Point, side: int, tau, h2):
    """h^-2 U(r_h + s h^2 tau) = D(tau)/r^4 (arb; h2 an arb interval)."""
    rh = pt.rp if side > 0 else pt.rc
    d = _shift(pt.dco, rh)
    ta = arb(tau)
    D = sum((d[k] * side ** k * h2 ** (k - 1) * ta ** k for k in range(1, 5)), arb(0))
    r = rh + side * h2 * ta
    return D / r ** 4


def h2_pieces():
    return [(H0SQ * j / K_H, H0SQ * (j + 1) / K_H) for j in range(K_H)]


def scaled_U_min(pt: Point, side: int, tau) -> arb:
    """Certified lower bound of D/r^4 over h^2 in [0, h0^2] (min over the K_H sub-intervals)."""
    lo_ = None
    for a, b in h2_pieces():
        v = scaled_U(pt, side, tau, _hull(fa(a), fa(b))).lower()
        lo_ = v if lo_ is None or v < lo_ else lo_
    return arb(lo_)


def exit_tau(pt: Point, side: int, Mbox: arb, tau_max=Fraction(4096)):
    """Smallest multiple of TAU_STEP (found by doubling and bisection on the grid) with
    min D/r^4 >= Mbox over h^2 in [0, h0^2]; the returned point is re-checked."""
    def ok(tau):
        return bool(scaled_U_min(pt, side, fmpq(tau.numerator, tau.denominator)) >= Mbox)
    hi_ = TAU_STEP
    while not ok(hi_):
        hi_ *= 2
        if hi_ > tau_max:
            raise ArithmeticError("exit not found")
    lo_ = hi_ / 2 if hi_ > TAU_STEP else Fraction(0)
    while hi_ - lo_ > TAU_STEP:
        mid = Fraction(int(((lo_ + hi_) / 2) / TAU_STEP)) * TAU_STEP
        if mid <= lo_:
            mid = lo_ + TAU_STEP
        if mid >= hi_:
            break
        if ok(mid):
            hi_ = mid
        else:
            lo_ = mid
    assert ok(hi_)
    return hi_


def exit_tau_point(pt: Point, side: int, Mbox: arb, h2v: Fraction, tau_max=Fraction(2 ** 20)):
    """The exit rule at one value of h^2 (N3, K2): smallest multiple of TAU_STEP with D/r^4 >= Mbox."""
    def ok(tau):
        return bool(scaled_U(pt, side, fmpq(tau.numerator, tau.denominator), fa(h2v)) >= Mbox)
    hi_ = TAU_STEP
    while not ok(hi_):
        hi_ *= 2
        if hi_ > tau_max:
            raise ArithmeticError("exit not found")
    lo_ = hi_ / 2 if hi_ > TAU_STEP else Fraction(0)
    while hi_ - lo_ > TAU_STEP:
        mid = Fraction(int(((lo_ + hi_) / 2) / TAU_STEP)) * TAU_STEP
        if mid <= lo_:
            mid = lo_ + TAU_STEP
        if mid >= hi_:
            break
        if ok(mid):
            hi_ = mid
        else:
            lo_ = mid
    return hi_


def y_exit(pt: Point, side: int, wbox, h2_lo: Fraction, h2_hi: Fraction, tau: Fraction):
    """Affine enclosure of y = Phi_x/Phi at r = r_h + s h^2 tau over the omega box and the h^2 piece."""
    w, h2 = params_af(pt, wbox, h2_lo, h2_hi)
    kap = pt.kp if side > 0 else pt.kc
    rho = w * (-acb(0, 1) / (2 * kap))
    A, B, C, Dc, (r0, r1) = tau_polys(pt, side, w, h2)
    tq = fmpq(tau.numerator, tau.denominator)
    S, Sd, info = frob_tau(A, B, C, rho, tq)
    ta = arb(tq)
    Dt = Dc[1] * ta
    for k in range(2, 5):
        Dt = Dt + Dc[k] * ta ** k
    r = r0 + r1 * ta
    y = h2 * Dt / (r * r * r) + (Dt / (r * r)) * (Sd / S) * side
    return y, info


class LChecker:
    def __init__(self, pt: Point, W_bound: Fraction):
        self.pt = pt
        self.CW = 4 * arb(fa(W_bound)).max(pt.kp * pt.kp)

    def side_check(self, side: int, wbox, tau: Fraction):
        """min (side +) / max (side c) over the h^2 pieces of the Re y bounds."""
        worst = None
        terms = 0
        for a, b in h2_pieces():
            y, info = y_exit(self.pt, side, wbox, a, b, tau)
            v = y.re_lower() if side > 0 else y.re_upper()
            terms = max(terms, info["terms"])
            if worst is None or (v < worst if side > 0 else v > worst):
                worst = v
        return worst, terms

    def box(self, wbox, exit_scale=None):
        """omega/kappa_- in wbox = (re_lo, re_hi, im_lo, im_hi).  exit_scale (N5 only) places the exit
        at D/r^4 >= exit_scale (Re w)_min^2 instead of the exit rule."""
        pt = self.pt
        km = pt.km
        re_lo, re_hi, im_lo, im_hi = wbox
        re_ = _hull(fa(re_lo) * km, fa(re_hi) * km)
        im_ = _hull(fa(im_lo) * km, fa(im_hi) * km)
        wabs2 = up(abs(acb(re_, im_)) ** 2)
        M = up((4 * wabs2).max(self.CW))
        rec = {"box": [str(x) for x in wbox], "M": M.str(8)}
        target = M
        if exit_scale is not None:
            target = arb((fa(exit_scale) * (fa(re_lo) * km) ** 2).lower())
            rec["exit_target"] = target.str(6)
        ok = True
        for side, name in ((1, "plus"), (-1, "c")):
            try:
                tau = exit_tau(pt, side, target)
                y, terms = self.side_check(side, wbox, tau)
                good = bool(y > 0) if side > 0 else bool(y < 0)
                Umin = scaled_U_min(pt, side, fmpq(tau.numerator, tau.denominator))
                exit_rule = bool(Umin >= M)
                rec[name] = {"tau_a": str(tau), "Re_y_bound": y.str(6), "y_sign_ok": good,
                             "exit_rule_ok": exit_rule, "terms": terms}
                if exit_scale is not None:
                    # the exit is in the allowed region: D/r^4 < (Re w)_min^2 at tau_a for all h^2
                    Umax = max((scaled_U(pt, side, fmpq(tau.numerator, tau.denominator), _hull(fa(a), fa(b))).upper()
                                for a, b in h2_pieces()), key=float)
                    rec[name]["exit_in_allowed_region"] = bool(arb(Umax) < (fa(re_lo) * km) ** 2)
                ok = ok and good and exit_rule
            except (ArithmeticError, ZeroDivisionError, ValueError) as ex:
                rec[name] = {"error": str(ex)}
                ok = False
        rec["poles_avoided"] = bool(im_.lower() > -pt.kp)      # kappa_+ < kappa_c
        rec["pass"] = bool(ok and rec["poles_avoided"])
        return rec


def boxes_L():
    """Production omega-boxes: [0, 8] x [-11/20, 0] (kappa_- units), side 1/20."""
    out = []
    n_re = int(Fraction(8) / BOX)
    n_im = int(Fraction(11, 20) / BOX)
    for i in range(n_re):
        for j in range(n_im):
            out.append((i * BOX, (i + 1) * BOX, -(j + 1) * BOX, -j * BOX))
    return out


def certify_box_refined(ch: LChecker, b, min_side=MIN_BOX):
    """Certify a box; on failure halve it in both directions, down to min_side."""
    rec = ch.box(b)
    if rec["pass"] or (b[1] - b[0]) / 2 < min_side:
        return [rec]
    out = []
    mr, mi = (b[0] + b[1]) / 2, (b[2] + b[3]) / 2
    for sb in ((b[0], mr, b[2], mi), (mr, b[1], b[2], mi), (b[0], mr, mi, b[3]), (mr, b[1], mi, b[3])):
        out.extend(certify_box_refined(ch, sb, min_side))
    return out


# ---------------------------------------------------------------------------------------------
# continuation path (N3) and Wronskian values (K2)
# ---------------------------------------------------------------------------------------------
def jost_at(pt: Point, side: int, l: int, w: acb, tau: Fraction):
    """(R, dR/dr) of the FORMULATION-normalised Jost solution at r = r_h + s h^2 tau, h = 1/(l + 1/2),
    w an acb ball/box: R = (h^2 tau)^rho S(tau), dR/dr = (h^2 tau)^rho s h^-2 Sd(tau)."""
    h2v = Fraction(1) / (Fraction(l) + Fraction(1, 2)) ** 2
    Rw = up(abs(acb(w.real.rad(), w.imag.rad())) + arb(2) ** -200)
    tiny = up(arb(2) ** -200)
    wa = AF(acb(w.real.mid(), w.imag.mid()), 1, 0, 0, Rw, tiny)
    h2a = AF(fa(h2v), 0, 1, 0, Rw, tiny)
    kap = pt.kp if side > 0 else pt.kc
    rho = wa * (-acb(0, 1) / (2 * kap))
    A, B, C, _Dc, _r = tau_polys(pt, side, wa, h2a)
    tq = fmpq(tau.numerator, tau.denominator)
    S, Sd, info = frob_tau(A, B, C, rho, tq)
    rho_b = -acb(0, 1) * w / (2 * kap)
    pref = (rho_b * acb(fa(h2v) * arb(tq)).log()).exp()
    return pref * S.to_acb(), pref * side * Sd.to_acb() / fa(h2v), info


def _shift_af(P, c: arb):
    """Coefficients of P(c + s) (P ascending, affine-form coefficients; c exact)."""
    from math import comb
    n = len(P)
    out = []
    for k in range(n):
        acc = P[0].const(0)
        cp = arb(1)
        for j in range(k, n):
            acc = acc + P[j] * (comb(j, k) * cp)
            cp = cp * c
        out.append(acc)
    return out


def taylor_tau(A, B, C, tau_c: Fraction, R, dR, step: Fraction, min_terms=12):
    """Continue (R, R_tau) of A R'' + B R' + C R = 0 (tau-polynomials, affine forms) from tau_c to
    tau_c + step at an ordinary point; recurrence and majorant tail as arbw.taylor_step."""
    c = fa(tau_c)
    As, Bs, Cs = _shift_af(A, c), _shift_af(B, c), _shift_af(C, c)
    dd = fa(step)
    bco = [R, dR]
    S = R + dR * dd
    Sd = dR * 1
    dk = dd
    small = 0
    eps = arb(2) ** (-(ctx.prec - 12))
    A0 = As[0]
    for n in range(N_MAX):
        acc = R.const(0)
        for j in range(1, min(len(As), n + 3)):
            m = n - j + 2
            acc = acc + As[j] * (m * (m - 1)) * bco[m]
        for j in range(min(len(Bs), n + 2)):
            m = n - j + 1
            acc = acc + Bs[j] * m * bco[m]
        for j in range(min(len(Cs), n + 1)):
            acc = acc + Cs[j] * bco[n - j]
        bn2 = -(acc / (A0 * ((n + 2) * (n + 1))))
        bco.append(bn2)
        dk = dk * dd
        S = S + bn2 * dk
        Sd = Sd + bn2 * ((n + 2) * (dk / dd))
        if n >= min_terms and abs(bn2.c0).upper() * abs(dk) < eps * abs(S.c0).lower():
            small += 1
            if small >= 6:
                break
        else:
            small = 0
    M = len(bco)
    N0 = M - 2
    A0l = abs(A0.c0).lower() - A0.dev()
    if not (A0l > 0):
        raise ArithmeticError("taylor_tau: A_0 not bounded below")
    g = []
    for j in range(1, 10):
        Aj = As[j].absu() if j < len(As) else arb(0)
        Bj = Bs[j - 1].absu() if j - 1 < len(Bs) else arb(0)
        Cj = Cs[j - 2].absu() if 0 <= j - 2 < len(Cs) else arb(0)
        g.append(up((Aj + Bj / (N0 + 1) + Cj / ((N0 + 1) * (N0 + 2))) / A0l))
    q = arb(_solve_q(g))
    if not (sum((gj * q ** (-(j + 1)) for j, gj in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("taylor_tau majorant root not verified")
    K = max((bco[m].absu() / q ** m for m in range(M - 9, M)), key=lambda x: x.upper())
    s0, s1 = _geom_tails(up(K), q, abs(dd), M, arb(0))
    S = S + AF(0, 0, 0, s0, S.Rw, S.Rh)
    Sd = Sd + AF(0, 0, 0, s1, S.Rw, S.Rh)
    return S, Sd


def y_continued_tau(pt: Point, side: int, wbox_acb: acb, h2v: Fraction, tau0: Fraction, tau_a: Fraction):
    """y at r_h + s h^2 tau_a for one value of h^2 and an omega ball/box: scaled Frobenius to tau0, then
    affine Taylor steps in tau (step <= tau/2 and <= a quarter of the distance to the other roots)."""
    w = wbox_acb
    Rw = up(abs(acb(w.real.rad(), w.imag.rad())) + arb(2) ** -200)
    tiny = up(arb(2) ** -200)
    wa = AF(acb(w.real.mid(), w.imag.mid()), 1, 0, 0, Rw, tiny)
    h2a = AF(fa(h2v), 0, 1, 0, Rw, tiny)
    kap = pt.kp if side > 0 else pt.kc
    rho = wa * (-acb(0, 1) / (2 * kap))
    A, B, C, Dc, (r0, r1) = tau_polys(pt, side, wa, h2a)
    S, Sd, _ = frob_tau(A, B, C, rho, fmpq(tau0.numerator, tau0.denominator))
    # (R, R_tau) up to the common factor tau0^rho
    rh = pt.rp if side > 0 else pt.rc
    dist = min(float(abs(rh - ri).lower()) for ri in pt.roots if not ri.overlaps(rh)) / float(h2v)
    tau = tau0
    R, dR = S, Sd
    nsteps = 0
    while tau < tau_a:
        # the Delta^2-type majorant at tau_c has radius about (sqrt 2 - 1) tau_c: steps of tau_c/4
        st = min(tau / 4, Fraction(int(dist / 4 * 1024), 1024), tau_a - tau)
        st = Fraction(int(st * 2 ** 20), 2 ** 20) if st != tau_a - tau else st
        while True:
            try:
                R2, dR2 = taylor_tau(A, B, C, tau, R, dR, st)
                break
            except ArithmeticError:
                st = st / 2
                if st < Fraction(1, 2 ** 24):
                    raise
        R, dR = R2, dR2
        tau = tau + st
        nsteps += 1
    ta = fa(tau_a)
    Dt = Dc[1] * ta
    for k in range(2, 5):
        Dt = Dt + Dc[k] * ta ** k
    r = r0 + r1 * ta
    y = h2a * Dt / (r * r * r) + (Dt / (r * r)) * (dR / R) * side
    return y, nsteps


def y_continued(pt: Point, side: int, l: int, w: acb, tau0: Fraction, r_exit):
    """y at r_exit for integer l: scaled Frobenius to tau0, then arbw Taylor steps in r (N3 path)."""
    from .arbw import path, taylor_step
    L = l * (l + 1)
    h2v = Fraction(1) / (Fraction(l) + Fraction(1, 2)) ** 2
    R, dR, _ = jost_at(pt, side, l, w, tau0)
    rh = pt.rp if side > 0 else pt.rc
    c = rh + side * fa(h2v * tau0)
    for nxt in path(pt, c, r_exit, 0.5, None):
        R, dR, _ = taylor_step(pt, c, R, dR, w, L, arb(nxt) - c)
        c = arb(nxt)
    r = arb(r_exit)
    F = pt.F(r)
    return F * (1 / r + dR / R)


def wronskian_L(pt: Point, l: int, w: acb, rmatch: fmpq | None = None):
    """W_l(w) (FORMULATION section 1.2 normalisation) from the scaled Frobenius values, continued by
    arbw Taylor steps to the matching point (K2)."""
    from .arbw import _rat, path, taylor_step
    L = l * (l + 1)
    h2v = Fraction(1) / (Fraction(l) + Fraction(1, 2)) ** 2
    rm = _rat(pt.rP.mid(), 20) if rmatch is None else rmatch
    out = {}
    for side in (1, -1):
        # start where the scaled series is comfortably inside its disc: t = 0.02 (r_+) / 0.2 (r_c)
        tstart = Fraction(1, 50) if side > 0 else Fraction(1, 5)
        tau0 = tstart / h2v
        R, dR, _ = jost_at(pt, side, l, w, tau0)
        rh = pt.rp if side > 0 else pt.rc
        c = rh + side * fa(tstart)
        for nxt in path(pt, c, rm, 0.5, None):
            R, dR, _ = taylor_step(pt, c, R, dR, w, L, arb(nxt) - c)
            c = arb(nxt)
        out[side] = (R, dR)
    (Rp, dRp), (Rc, dRc) = out[1], out[-1]
    rma = arb(rm)
    W = acb(pt.Delta(rma)) * (Rp * dRc - dRp * Rc)
    return W * (-acb(0, 1) * w * acb(pt.c_plus) + acb(0, 1) * w * acb(pt.c_cos)).exp() / acb(pt.rp * pt.rc)


# ---------------------------------------------------------------------------------------------
# resumable runs
# ---------------------------------------------------------------------------------------------
def run_L(worker: int = 0, nworkers: int = 1):
    from .part1 import load_bulk
    W_bound, _ = load_bulk()
    RES_DIR.mkdir(parents=True, exist_ok=True)
    out = RES_DIR / f"L_PB.w{worker}.jsonl"
    done = set()
    if out.exists():
        for line in out.read_text().splitlines():
            if line.strip():
                done.add(tuple(json.loads(line)["box"]))
    pt = Point("PB", 256)
    with precision(PREC):
        ch = LChecker(pt, W_bound)
        for k, b in enumerate(boxes_L()):
            key = tuple(str(x) for x in b)
            if k % nworkers != worker or key in done:
                continue
            with Timer() as tm:
                leaves = certify_box_refined(ch, b)
            row = {"box": list(key), "pass": all(x["pass"] for x in leaves), "leaves": leaves,
                   "commit": git_commit(), "timing_s": {"wall": round(tm.wall, 2), "cpu": round(tm.cpu, 2)},
                   "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}
            with out.open("a") as f:
                f.write(json.dumps(row) + "\n")
            if not row["pass"]:
                print("L FAIL", key, flush=True)


def assemble():
    rows = {}
    for f in sorted(RES_DIR.glob("L_PB.w*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[tuple(r["box"])] = r
    want = [tuple(str(x) for x in b) for b in boxes_L()]
    missing = [b for b in want if b not in rows]
    leaves = [x for b in want if b in rows for x in rows[b]["leaves"]]
    ok = (not missing) and all(rows[b]["pass"] for b in want if b in rows)
    commits = sorted({rows[b]["commit"] for b in want if b in rows})
    cpu = round(sum(rows[b]["timing_s"]["cpu"] for b in want if b in rows), 1)
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    certs = {}
    for name, sgn in (("plus", 1), ("c", -1)):
        vals = [x[name] for x in leaves if "Re_y_bound" in x.get(name, {})]
        side_ok = (not missing) and len(vals) == len(leaves) and all(v["y_sign_ok"] and v["exit_rule_ok"] for v in vals)
        ext = [float(arb(v["Re_y_bound"].split("+/-")[0].strip("[ ")).mid()) for v in vals]
        cert = {"regime": "L", "side": name, "point": "PB", "h0": str(H0), "precision_bits": PREC,
                "design": DESIGN_L, "boxes": len(want), "leaves": len(leaves), "missing": missing,
                ("min_Re_y_plus_lower" if sgn > 0 else "max_Re_y_c_upper"): (min(ext) if sgn > 0 else max(ext)) if ext else None,
                "max_tau_a": str(max((Fraction(v["tau_a"]) for v in vals), default=0)),
                "poles_avoided_all": all(x["poles_avoided"] for x in leaves),
                "cpu_s_total_both_sides": cpu, "commits": commits,
                "pass": bool(side_ok and all(x["poles_avoided"] for x in leaves)),
                "leaves_detail": [{"box": x["box"], "M": x["M"], **x[name]} for x in leaves]}
        (CERT_DIR / f"L_PB_{name}.json").write_text(json.dumps(cert, indent=1) + "\n")
        certs[name] = cert
    return {"pass": ok and all(c["pass"] for c in certs.values()), "plus": certs["plus"]["pass"], "c": certs["c"]["pass"]}
