"""Validated Wronskian W_l(omega) in Arb (ball arithmetic, outward rounding).

For an omega given as an acb ball (a point or a whole box), ``wronskian`` returns an acb ball
containing W_l(omega') for EVERY omega' in the ball (FORMULATION section 1.2 normalisation,
including the Jost factor exp(-i w c_+ + i w c_c)/(r_+ r_c)).

Ingredients (all enclosures are uniform over the omega ball):
  * Frobenius series at r_+ and r_c, R = |t|^rho sum a_k t^k, rho = -i w/(2 kappa), with the
    exact relations Delta(r_h) = 0 (so A_0 = A_1 = B_0 = 0 exactly) and the factored indicial
    polynomial ind(n) = Delta'(r_h)^2 n (n + 2 rho).  Tail: majorant recursion.  For n >= N,
      |a_n| <= sum_{j=1}^{6} g_j |a_{n-j}|,
      g_j = [ |A_{j+2}| (1+R/N)^2/(1-2R/N) + |B_{j+1}| (N+R)/(N(N-2R)) + |C_j|/(N(N-2R)) ] / |A_2|,
    R >= |rho| (each factor is decreasing in n, so the bound holds for all n >= N).  With q > 0
    such that sum g_j q^{-j} <= 1 (checked in Arb), induction gives |a_n| <= K q^n for n >= N - 6,
    K = max_{N-6 <= m < N} |a_m| q^{-m}; the tails of sum a_n t^n and sum (n+rho) a_n t^{n-1} are
    bounded by geometric series in x = q|t| < 1.
  * Taylor continuation at ordinary points c (Delta(c) != 0): b_{n+2} from the recurrence of
    A y'' + B y' + C y = 0; for n >= N0, |b_{n+2}| <= sum_{j=1}^{8} g_j |b_{n+2-j}|,
      g_j = ( |A_j| + |B_{j-1}|/(N0+1) + |C_{j-2}|/((N0+1)(N0+2)) ) / |A_0|,
    and the same geometric tail argument (window of 8 coefficients).
  A = Delta^2, B = Delta Delta', C = w^2 r^4 - L Delta  (the radial ODE times Delta).
"""
from __future__ import annotations

from fractions import Fraction
from math import comb

import mpmath as mp
from flint import acb, arb, ctx, fmpq

from .arbutil import precision
from .params import Point

N_FROB_MAX = 4000
N_TAYLOR_MAX = 4000


def _shift(dco, c):
    """Ascending coefficients of Delta(c + t) for Delta = sum dco_k r^k (arb), c arb."""
    co = list(dco)
    n = len(co)
    out = []
    for k in range(n):
        # value of the k-th derivative / k!
        acc = arb(0)
        for i in reversed(range(k, n)):
            acc = acc * c + co[i] * _binom(i, k)
        out.append(acc)
    return out


def _binom(i, k):
    return comb(i, k)


def _pmul(x, y):
    out = [acb(0)] * (len(x) + len(y) - 1)
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            out[i + j] += xi * yj
    return out


def ode_polys(pt: Point, c, w, L, singular=False):
    d = [acb(v) for v in _shift(pt.dco, c)]
    if singular:
        d[0] = acb(0)                      # Delta(r_h) = 0 exactly
    dp = [d[k + 1] * (k + 1) for k in range(4)]
    A = _pmul(d, d)
    B = _pmul(d, dp)
    c4 = [c ** 4, 4 * c ** 3, 6 * c ** 2, 4 * c, arb(1)]
    C = [w * w * acb(c4[k]) - L * d[k] for k in range(5)]
    if singular:
        A[0] = A[1] = B[0] = acb(0)
    return A, B, C, d


def _absu(z):
    return arb(abs(z).upper())


def _solve_q(g):
    """Small float q with sum g_j q^-j <= 1 (then verified in Arb by the caller)."""
    gf = [float(x.upper()) for x in g]

    def f(q):
        return sum(gj * q ** -(j + 1) for j, gj in enumerate(gf) if gj > 0)

    hi_ = 1.0
    while f(hi_) > 1:
        hi_ *= 2
    while hi_ > 1e-200 and f(hi_ / 2) <= 1:
        hi_ /= 2
    lo_ = hi_ / 2
    for _ in range(100):
        mid = (lo_ + hi_) / 2
        if f(mid) > 1:
            lo_ = mid
        else:
            hi_ = mid
        if hi_ - lo_ < 1e-13 * hi_:
            break
    return hi_ * (1 + 1e-9)


def _geom_tails(K, q, t_abs, M, rho_abs):
    """Bounds for sum_{n>=M} K q^n |t|^n and sum_{n>=M} (n+|rho|) K q^n |t|^(n-1)."""
    x = q * t_abs
    if not (x < 1):
        raise ArithmeticError("majorant ratio q|t| >= 1")
    xM = x ** M
    s0 = K * xM / (1 - x)
    s1 = K / t_abs * xM * ((M + rho_abs) / (1 - x) + x / (1 - x) ** 2)
    return s0, s1


def frobenius(pt: Point, which: str, w, L, t, min_terms=20):
    """(R, dR/dr) at r = r_h + t for the Frobenius solution at r_h in {'+', 'c'}; t exact (fmpq)."""
    rh = pt.rp if which == "+" else pt.rc
    A, B, C, _d = ode_polys(pt, rh, w, L, singular=True)
    kap = pt.kp if which == "+" else pt.kc
    rho = -acb(0, 1) * w / (2 * kap)
    A2 = A[2]
    ta = arb(t)
    tabs = abs(ta)
    a = [acb(1)]
    S = acb(1)
    Sd = rho / ta
    tk = arb(1)
    small = 0
    eps = arb(2) ** (-(ctx.prec - 8))
    for n in range(1, N_FROB_MAX):
        acc = acb(0)
        for k in range(max(0, n - 6), n):
            j2, j1, j0 = n - k + 2, n - k + 1, n - k
            term = acb(0)
            if j2 < len(A):
                term += A[j2] * (k + rho) * (k + rho - 1)
            if j1 < len(B):
                term += B[j1] * (k + rho)
            if j0 < len(C):
                term += C[j0]
            acc += a[k] * term
        an = -acc / (A2 * n * (n + 2 * rho))
        a.append(an)
        tk = tk * ta
        tv = an * tk
        S += tv
        Sd += an * (n + rho) * tk / ta
        if n >= min_terms and _absu(tv) < eps * abs(S).lower():
            small += 1
            if small >= 8:
                break
        else:
            small = 0
    M = len(a)                                   # computed a_0 .. a_{M-1}
    # majorant constants for n >= M
    Rr = _absu(rho)
    Nn = arb(M)
    if not (Nn > 2 * Rr + 1):
        raise ArithmeticError("Frobenius: N too small relative to |rho|")
    A2l = abs(A2).lower()
    g = []
    for j in range(1, 7):
        Aj = _absu(A[j + 2]) if j + 2 < len(A) else arb(0)
        Bj = _absu(B[j + 1]) if j + 1 < len(B) else arb(0)
        Cj = _absu(C[j]) if j < len(C) else arb(0)
        gj = (Aj * (1 + Rr / Nn) ** 2 / (1 - 2 * Rr / Nn) + Bj * (Nn + Rr) / (Nn * (Nn - 2 * Rr))
              + Cj / (Nn * (Nn - 2 * Rr))) / A2l
        g.append(arb(gj.upper()))
    qf = _solve_q(g)
    q = arb(qf)
    if not (sum((gj * q ** (-(j + 1)) for j, gj in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("Frobenius majorant root not verified")
    K = max((_absu(a[m]) / q ** m for m in range(M - 6, M)), key=lambda x: x.upper())
    K = arb(K.upper())
    s0, s1 = _geom_tails(K, q, tabs, M, Rr)
    S += acb(arb(0, s0.upper()), arb(0, s0.upper()))
    Sd += acb(arb(0, s1.upper()), arb(0, s1.upper()))
    pref = (rho * acb(tabs).log()).exp()       # |t|^rho
    return pref * S, pref * Sd, {"terms": M, "q": qf, "tail": float(s0.upper())}


def taylor_step(pt: Point, c, R, dR, w, L, dstep, min_terms=10):
    """Continue (R, R') from c to c + dstep (c, dstep arb; ordinary points)."""
    A, B, C, _ = ode_polys(pt, c, w, L)
    eps = arb(2) ** (-(ctx.prec - 8))
    bco = [R, dR]
    dd = acb(dstep)
    S = R + dR * dd
    Sd = acb(dR)
    dk = dd
    small = 0
    A0 = A[0]
    for n in range(N_TAYLOR_MAX):
        acc = acb(0)
        for j in range(1, min(len(A), n + 3)):
            m = n - j + 2
            if m >= 0:
                acc += A[j] * (m * (m - 1)) * bco[m]
        for j in range(min(len(B), n + 2)):
            m = n - j + 1
            if m >= 0:
                acc += B[j] * m * bco[m]
        for j in range(min(len(C), n + 1)):
            acc += C[j] * bco[n - j]
        bn2 = -acc / (A0 * ((n + 2) * (n + 1)))
        bco.append(bn2)
        dk = dk * dd
        tv = bn2 * dk
        S += tv
        Sd += (n + 2) * bn2 * dk / dd
        if n >= min_terms and _absu(tv) < eps * (abs(S) + abs(Sd * dd)).lower():
            small += 1
            if small >= 8:
                break
        else:
            small = 0
    M = len(bco)                                 # b_0 .. b_{M-1}
    N0 = M - 2
    A0l = abs(A0).lower()
    g = []
    for j in range(1, 9):
        Aj = _absu(A[j]) if j < len(A) else arb(0)
        Bj = _absu(B[j - 1]) if j - 1 < len(B) else arb(0)
        Cj = _absu(C[j - 2]) if 0 <= j - 2 < len(C) else arb(0)
        g.append(arb(((Aj + Bj / (N0 + 1) + Cj / ((N0 + 1) * (N0 + 2))) / A0l).upper()))
    qf = _solve_q(g)
    q = arb(qf)
    if not (sum((gj * q ** (-(j + 1)) for j, gj in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("Taylor majorant root not verified")
    K = max((_absu(bco[m]) / q ** m for m in range(M - 8, M)), key=lambda x: x.upper())
    K = arb(K.upper())
    s0, s1 = _geom_tails(K, q, abs(arb(dstep)), M, arb(0))
    S += acb(arb(0, s0.upper()), arb(0, s0.upper()))
    Sd += acb(arb(0, s1.upper()), arb(0, s1.upper()))
    return S, Sd, {"terms": M, "q": qf, "x": float((q * abs(arb(dstep))).upper())}


def _rat(x, bits=24):
    """Exact dyadic rational near the float x."""
    return fmpq(round(float(x) * 2 ** bits), 2 ** bits)


def majorant_radius(c, roots):
    """x* with prod_i (1 + x/|c - r_i|)^2 = 2: the radius up to which the Delta^2 majorant of a
    Taylor step at c can certify convergence (float; only used to choose step sizes)."""
    ds = [abs(c - x) for x in roots]
    f = lambda x: 1.0
    lo_, hi_ = 0.0, min(ds)
    for _ in range(100):
        mid = (lo_ + hi_) / 2
        f = 1.0
        for d in ds:
            f *= (1 + mid / d) ** 2
        if f > 2:
            hi_ = mid
        else:
            lo_ = mid
    return lo_


def path(pt: Point, r_start, r_end, frac, kfac):
    """Exact-rational step endpoints from r_start (arb) to r_end (fmpq): step <= frac * x*(r)
    (majorant radius) and <= kfac (float design choices; rigor is checked in each step)."""
    pts = []
    rs = [float(x.mid()) for x in pt.roots]
    r = float(r_start.mid())
    target = float(arb(r_end).mid())
    while True:
        h = frac * majorant_radius(r, rs)
        if kfac is not None:
            h = min(h, kfac)
        if abs(target - r) <= h * 1.0000001:
            pts.append(fmpq(r_end))
            return pts
        r = r + h if target > r else r - h
        pts.append(_rat(r))
        r = float(arb(pts[-1]).mid())


def wronskian(pt: Point, l: int, w, rmatch=None, t_frac=Fraction(1, 5), frac=0.5, jost=True, info=None):
    """Enclosure of W_l over the acb ball w (absolute units, not kappa_- units)."""
    L = l * (l + 1)
    w = acb(w)
    rm_ = pt.rP if rmatch is None else arb(rmatch)
    rmatch_q = _rat(rm_.mid(), 20) if rmatch is None else fmpq(rmatch)
    rs = [float(x.mid()) for x in pt.roots]
    # step-size control by the majorant radius and the local frequency (design only)
    wabs = float(abs(w).upper())
    xp = majorant_radius(rs[2], [rs[0], rs[1], rs[3]])
    xc = majorant_radius(rs[3], [rs[0], rs[1], rs[2]])
    rho_p = wabs / (2 * float(pt.kp.mid()))
    rho_c = wabs / (2 * float(pt.kc.mid()))
    tp = _rat(min(float(t_frac) * 2.5 * xp, 2.5 * xp / max(1.0, rho_p)), 30)
    tc = _rat(min(float(t_frac) * 2.5 * xc, 2.5 * xc / max(1.0, rho_c)), 30)
    Rp, dRp, i1 = frobenius(pt, "+", w, L, tp)
    c = pt.rp + arb(tp)
    stats = {"frob_p": i1, "steps": []}
    for nxt in path(pt, c, rmatch_q, frac, None):
        h = arb(nxt) - c
        Rp, dRp, inf = taylor_step(pt, c, Rp, dRp, w, L, h)
        stats["steps"].append(inf)
        c = arb(nxt)
    Rc, dRc, i2 = frobenius(pt, "c", w, L, -tc)
    stats["frob_c"] = i2
    c = pt.rc - arb(tc)
    for nxt in path(pt, c, rmatch_q, frac, None):
        h = arb(nxt) - c
        Rc, dRc, inf = taylor_step(pt, c, Rc, dRc, w, L, h)
        stats["steps"].append(inf)
        c = arb(nxt)
    rmatch_a = arb(rmatch_q)
    W = acb(pt.Delta(rmatch_a)) * (Rp * dRc - dRp * Rc)
    if jost:
        W = W * (-acb(0, 1) * w * acb(pt.c_plus) + acb(0, 1) * w * acb(pt.c_cos)).exp() / acb(pt.rp * pt.rc)
    if info is not None:
        info.update(stats)
    return W


def mp_to_acb(z):
    return acb(arb(mp.nstr(mp.re(z), 60)), arb(mp.nstr(mp.im(z), 60)))


def w_at(name: str, l: int, w_over_km, prec=256):
    """Convenience: enclosure of W_l at omega = w_over_km * kappa_- (point ball)."""
    pt = Point(name, prec)
    with precision(prec):
        w = acb(w_over_km) * acb(pt.km)
        return wronskian(pt, l, w)
