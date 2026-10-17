"""Validated W_l(omega) from the phase-factored Jost solutions (S3, part iii).

Factored solutions.  With r* = sum_i A_i log|r - r_i| (A_+ = 1/(2 kappa_+), A_c = -1/(2 kappa_c)) put
    Z_+ := exp(+i w (r* - c_+)) R_+ ,    Z_c := exp(-i w (r* - c_c)) R_c ,
where R_+ = (r - r_+)^{-i w/(2 kappa_+)} S_+(r), R_c = (r_c - r)^{-i w/(2 kappa_c)} S_c(r) are the
FORMULATION section 1.2 Jost solutions (S_h analytic, S_h(r_h) = 1) and c_+, c_c are
`Point.c_plus`, `Point.c_cos`.  Then Z_+ = exp(i w Rt_+(r)) S_+(r) and Z_c = exp(-i w Rt_c(r)) S_c(r)
(Rt_h(r) = r* - A_h log|r - r_h| - c_h, Rt_h(r_h) = 0) are analytic at r_+ resp. r_c with value 1,
and they solve
    Delta Z'' + (Delta' - 2 i s w r^2) Z' - (2 i s w r + L) Z = 0 ,   s = +1 (Z_+), s = -1 (Z_c),
(substitute R = exp(-i s w r*) Z into (Delta R')' + (w^2 r^4/Delta - L) R = 0; r*' = r^2/Delta).
The exponential factors cancel exactly in the Jost-normalised Wronskian (arbw.wronskian):
    W_l(w) = Delta (R_+ R_c' - R_+' R_c) exp(-i w c_+ + i w c_c)/(r_+ r_c)
           = [Delta (Z_+ Z_c' - Z_+' Z_c) + 2 i w r^2 Z_+ Z_c] / (r_+ r_c)        (any r in (r_+, r_c)).
So no exp(+-i w r*) is ever evaluated, and over an omega-disc the propagated quantities vary only
through the potential (reflection, tunnelling), not through the free phase.

Series.  At r_h: Z = sum a_n t^n (t = r - r_h), a_0 = 1, indicial factor D_m = m (Delta'(r_h)(m - 1)
+ B_0) = Delta'(r_h) m (m - i s w/kappa_h sign), which vanishes only on E_+ resp. E_c (the poles of W).
Tail: for m >= N, |a_m| <= sum_{k=1}^{3} g_k |a_{m-k}|, g_k = (|A_{k+1}| m^2 + |B_k| m + |C_{k-1}|)/|D_m|,
each term decreasing in m (|D_m| >= |A_1| m (m - |z|), z the second exponent); then |a_n| <= K q^n.
Taylor steps at ordinary points: as arbw.taylor_step with the polynomials A = Delta, B, C above.

Enclosures over an omega-disc (Lohner frames with disc radii).  The centre trajectory y_c (point
omega_c) is propagated with the point step matrices T_c; the deviation y(w) - y_c over the disc is
kept as {B z : |z_j| <= e_j} (B a point frame, the parallelepiped of the centre propagation), with
    e' = |C| e + |B'^-1 d| + |B'^-1 r1| + |B'^-1 (T_box - T_c) B| e ,   C = B'^-1 T_c B ,
d = T(box) y_c - T_c y_c (vector enclosure over the box), r1 the rounding of the centre step.  All
quantities are balls; the frame choice only affects sharpness.
"""
from __future__ import annotations

from fractions import Fraction

from flint import acb, arb, ctx, fmpq

from .arbw import _geom_tails, _rat, _shift, _solve_q, majorant_radius
from .params import Point

N_MAX = 6000
I = acb(0, 1)


def _absu(z):
    return arb(abs(z).upper())


def _mid(z):
    return acb(z.real.mid(), z.imag.mid())


def zpolys(pt: Point, c, w, L, sgn, singular=False):
    """Ascending coefficients in t = r - c of A = Delta, B = Delta' - 2 i s w r^2, C = -(2 i s w r + L)."""
    d = [acb(v) for v in _shift(pt.dco, c)]
    if singular:
        d[0] = acb(0)                                   # Delta(r_h) = 0 exactly
    dp = [d[k + 1] * (k + 1) for k in range(4)]
    cc = acb(c)
    r2 = [cc * cc, 2 * cc, acb(1)]                      # (c + t)^2
    k = 2 * I * sgn * w
    B = [dp[0] - k * r2[0], dp[1] - k * r2[1], dp[2] - k * r2[2], dp[3]]
    C = [-(k * cc + L), -k]
    return d, B, C


def frob_z(pt: Point, which: str, w, L, t, min_terms=20):
    """(Z, Z') at r = r_h + t for the analytic solution at r_h (Z(r_h) = 1); t exact (fmpq)."""
    rh = pt.rp if which == "+" else pt.rc
    sgn = 1 if which == "+" else -1
    A, B, C = zpolys(pt, rh, acb(w), L, sgn, singular=True)
    A1 = A[1]
    zest = float((2 * abs(acb(w)) * abs(rh * rh) / abs(A1)).upper())
    min_terms = max(min_terms, int(2 * zest) + 8)        # the tail bound needs N > 2|z| + 2
    ta = arb(t)
    tabs = abs(ta)
    a = [acb(1)]
    S = acb(1)
    Sd = acb(0)
    tk = arb(1)
    small = 0
    eps = arb(2) ** (-(ctx.prec - 8))
    for m in range(1, N_MAX):
        acc = acb(0)
        for n in range(max(0, m - 3), m):
            j2, j1, j0 = m + 1 - n, m - n, m - 1 - n
            term = acb(0)
            if j2 < len(A):
                term += A[j2] * (n * (n - 1))
            if j1 < len(B):
                term += B[j1] * n
            if 0 <= j0 < len(C):
                term += C[j0]
            acc += a[n] * term
        Dm = m * (A1 * (m - 1) + B[0])
        am = -acc / Dm
        a.append(am)
        Sd += am * m * tk                                # uses t^(m-1)
        tk = tk * ta
        tv = am * tk
        S += tv
        if m >= min_terms and _absu(tv) < eps * abs(S).lower():
            small += 1
            if small >= 6:
                break
        else:
            small = 0
    M = len(a)
    # second exponent z = (2 i s w r_h^2)/A1 ; |D_m| >= |A1| m (m - |z|)
    z = 2 * I * sgn * acb(w) * acb(rh * rh) / A1
    za = _absu(z)
    Nn = arb(M)
    if not (Nn > 2 * za + 2):
        raise ArithmeticError("frob_z: N too small relative to the second exponent")
    A1l = abs(A1).lower()
    g = []
    for k in range(1, 4):
        Ak = _absu(A[k + 1]) if k + 1 < len(A) else arb(0)
        Bk = _absu(B[k]) if k < len(B) else arb(0)
        Ck = _absu(C[k - 1]) if k - 1 < len(C) else arb(0)
        gk = (Ak * Nn * Nn + Bk * Nn + Ck) / (A1l * Nn * (Nn - za))
        g.append(arb(gk.upper()))
    qf = _solve_q(g)
    q = arb(qf)
    if not (sum((gk * q ** (-(k + 1)) for k, gk in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("frob_z majorant root not verified")
    K = max((_absu(a[m]) / q ** m for m in range(M - 3, M)), key=lambda x: x.upper())
    K = arb(K.upper())
    s0, s1 = _geom_tails(K, q, tabs, M, arb(0))
    S += acb(arb(0, s0.upper()), arb(0, s0.upper()))
    Sd += acb(arb(0, s1.upper()), arb(0, s1.upper()))
    return S, Sd


def taylor_z(pt: Point, c, z0, dz0, w, L, h, sgn, min_terms=8):
    """Continue (Z, Z') of the s-equation from c to c + h (ordinary points; c, h arb)."""
    A, B, C = zpolys(pt, c, acb(w), L, sgn)
    eps = arb(2) ** (-(ctx.prec - 8))
    bco = [acb(z0), acb(dz0)]
    dd = acb(h)
    S = bco[0] + bco[1] * dd
    Sd = acb(bco[1])
    dk = dd
    small = 0
    A0 = A[0]
    for n in range(N_MAX):
        acc = acb(0)
        for j in range(1, min(len(A), n + 3)):
            m = n - j + 2
            acc += A[j] * (m * (m - 1)) * bco[m]
        for j in range(min(len(B), n + 2)):
            m = n - j + 1
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
            if small >= 6:
                break
        else:
            small = 0
    M = len(bco)
    N0 = M - 2
    A0l = abs(A0).lower()
    g = []
    for j in range(1, 5):
        Aj = _absu(A[j]) if j < len(A) else arb(0)
        Bj = _absu(B[j - 1]) if j - 1 < len(B) else arb(0)
        Cj = _absu(C[j - 2]) if 0 <= j - 2 < len(C) else arb(0)
        g.append(arb(((Aj + Bj / (N0 + 1) + Cj / ((N0 + 1) * (N0 + 2))) / A0l).upper()))
    qf = _solve_q(g)
    q = arb(qf)
    if not (sum((gj * q ** (-(j + 1)) for j, gj in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("taylor_z majorant root not verified")
    K = max((_absu(bco[m]) / q ** m for m in range(M - 4, M)), key=lambda x: x.upper())
    K = arb(K.upper())
    s0, s1 = _geom_tails(K, q, abs(arb(h)), M, arb(0))
    S += acb(arb(0, s0.upper()), arb(0, s0.upper()))
    Sd += acb(arb(0, s1.upper()), arb(0, s1.upper()))
    return S, Sd


# ---------------------------------------------------------------------------------------------
# path, scaling, frames

def _delta_f(rs, lam3, r):
    D = -lam3
    for x in rs:
        D *= (r - x)
    return D


def scale_at(pt: Point, r, wcf: complex, L) -> fmpq:
    """Exact positive rational s(r) for the coordinates (Z, Z'/s) (design only; any s > 0 is valid)."""
    rs = [float(x.mid()) for x in pt.roots]
    lam3 = float((pt.Lam / 3).mid())
    rf = float(arb(r).mid())
    D = abs(_delta_f(rs, lam3, rf))
    sv = ((2 * abs(wcf) * rf * rf) ** 2 + (L + 1) * D) ** 0.5 / D
    return _rat(sv, 30) if sv > 2 ** -20 else fmpq(1, 2 ** 20)


def path(pt: Point, r0, r1, wabs, frac=0.5, theta=1.0):
    """Exact-rational step endpoints: step <= frac * majorant radius and <= theta / (2|w| r^2/|Delta|)."""
    rs = [float(x.mid()) for x in pt.roots]
    lam3 = float((pt.Lam / 3).mid())
    pts = []
    r = float(arb(r0).mid())
    target = float(arb(r1).mid())
    while True:
        rate = 2 * wabs * r * r / abs(_delta_f(rs, lam3, r))
        h = frac * majorant_radius(r, rs)
        if rate > 0:
            h = min(h, theta / rate)
        if abs(target - r) <= h * 1.0000001:
            pts.append(fmpq(r1))
            return pts
        r = r + h if target > r else r - h
        pts.append(_rat(r, 40))
        r = float(arb(pts[-1]).mid())


def _mv(X, v):
    return [X[0][0] * v[0] + X[0][1] * v[1], X[1][0] * v[0] + X[1][1] * v[1]]


def _mm(X, Y):
    return [[X[i][0] * Y[0][j] + X[i][1] * Y[1][j] for j in range(2)] for i in range(2)]


def _inv2(X):
    det = X[0][0] * X[1][1] - X[0][1] * X[1][0]
    return [[X[1][1] / det, -X[0][1] / det], [-X[1][0] / det, X[0][0] / det]]


def _frame(X):
    """Point matrix: mid(X) with unit columns (the parallelepiped frame)."""
    out = [[None, None], [None, None]]
    for j in range(2):
        v = [_mid(X[0][j]), _mid(X[1][j])]
        n = arb(((abs(v[0]) ** 2 + abs(v[1]) ** 2).sqrt()).mid())
        out[0][j] = _mid(v[0] / n)
        out[1][j] = _mid(v[1] / n)
    return out


def _qr_frame(X):
    """Point matrix with orthonormal columns (Gram-Schmidt on mid(X), the longer column first)."""
    v = [[_mid(X[0][j]), _mid(X[1][j])] for j in range(2)]
    n = [arb(((abs(x[0]) ** 2 + abs(x[1]) ** 2).sqrt()).mid()) for x in v]
    if float(n[1]) > float(n[0]):
        v = [v[1], v[0]]
        n = [n[1], n[0]]
    q1 = [_mid(v[0][0] / n[0]), _mid(v[0][1] / n[0])]
    dot = q1[0].conjugate() * v[1][0] + q1[1].conjugate() * v[1][1]
    u = [v[1][0] - dot * q1[0], v[1][1] - dot * q1[1]]
    nu = arb(((abs(u[0]) ** 2 + abs(u[1]) ** 2).sqrt()).mid())
    q2 = [_mid(u[0] / nu), _mid(u[1] / nu)]
    return [[q1[0], q2[0]], [q1[1], q2[1]]]


def _choose_frame(TB, min_sin=0.1):
    """The parallelepiped frame while it is well conditioned (|det| of the unit-column frame, i.e.
    the sine of the angle between the columns, >= min_sin); an orthonormal frame otherwise."""
    Bn = _frame(TB)
    det = Bn[0][0] * Bn[1][1] - Bn[0][1] * Bn[1][0]
    if float(abs(det).mid()) >= min_sin:
        return Bn
    return _qr_frame(TB)


def _step_matrix(pt, c, w, L, h, sgn, a0, a1):
    R1, dR1 = taylor_z(pt, c, acb(1), acb(0), w, L, h, sgn)
    R2, dR2 = taylor_z(pt, c, acb(0), acb(1), w, L, h, sgn)
    return [[R1, R2 * a0], [dR1 / a1, dR2 * a0 / a1]]


def propagate(pt, which, wc, wbox, L, theta=1.0, frac=0.5, t_frac=Fraction(1, 5), rmatch=None, stats=None):
    """Centre (Z, Z') at the matching point for omega_c, and a radius dev with
    |Z(w) - Z_c| <= dev, |Z'(w) - Z'_c| <= dev * s_m for every w in the box wbox (or 0 if wbox is None).
    Returns (Zc, dZc, dev, s_m)."""
    sgn = 1 if which == "+" else -1
    wc = acb(wc)
    wcf = complex(float(wc.real.mid()), float(wc.imag.mid()))
    rm = _rat(pt.rP.mid(), 20) if rmatch is None else fmpq(rmatch)
    rs = [float(x.mid()) for x in pt.roots]
    wabs = float(abs(wbox if wbox is not None else wc).upper())
    if which == "+":
        xr = majorant_radius(rs[2], [rs[0], rs[1], rs[3]])
        kap = float(pt.kp.mid())
    else:
        xr = majorant_radius(rs[3], [rs[0], rs[1], rs[2]])
        kap = float(pt.kc.mid())
    zexp = wabs / kap                                   # |second exponent|
    t = _rat(min(float(t_frac) * 2.5 * xr, 2.5 * xr / max(1.0, zexp)), 30)
    if which == "c":
        t = -t
    rh = pt.rp if which == "+" else pt.rc
    start = rh + arb(t)
    Z, dZ = frob_z(pt, which, wc, L, t)
    s0 = scale_at(pt, start, wcf, L)
    y = [Z, dZ / arb(s0)]
    yc = [_mid(y[0]), _mid(y[1])]
    rnd0 = [y[0] - yc[0], y[1] - yc[1]]
    e = [_absu(rnd0[0]), _absu(rnd0[1])]
    if wbox is not None:
        Zb, dZb = frob_z(pt, which, wbox, L, t)
        e[0] = arb((e[0] + _absu(Zb - yc[0])).upper())
        e[1] = arb((e[1] + _absu(dZb / arb(s0) - yc[1])).upper())
    B = [[acb(1), acb(0)], [acb(0), acb(1)]]
    c = start
    pts = path(pt, start, rm, wabs, frac, theta)
    for nxt in pts:
        h = arb(nxt) - c
        s1 = scale_at(pt, nxt, wcf, L)
        a0, a1 = arb(s0), arb(s1)
        Tc = _step_matrix(pt, c, wc, L, h, sgn, a0, a1)
        yn = _mv(Tc, yc)
        yc_new = [_mid(yn[0]), _mid(yn[1])]
        r1 = [yn[0] - yc_new[0], yn[1] - yc_new[1]]
        TB = _mm(Tc, B)
        Bn = _choose_frame(TB)
        Bi = _inv2(Bn)
        Cm = _mm(Bi, TB)
        add = _mv(Bi, r1)
        en = [_absu(Cm[0][0]) * e[0] + _absu(Cm[0][1]) * e[1] + _absu(add[0]),
              _absu(Cm[1][0]) * e[0] + _absu(Cm[1][1]) * e[1] + _absu(add[1])]
        if wbox is not None:
            Tb = _step_matrix(pt, c, wbox, L, h, sgn, a0, a1)
            vb = taylor_z(pt, c, yc[0], yc[1] * a0, wbox, L, h, sgn)
            d = [vb[0] - yn[0], vb[1] / a1 - yn[1]]
            dB = _mv(Bi, d)
            DT = [[Tb[i][j] - Tc[i][j] for j in range(2)] for i in range(2)]
            G = _mm(Bi, _mm(DT, B))
            en[0] += _absu(dB[0]) + _absu(G[0][0]) * e[0] + _absu(G[0][1]) * e[1]
            en[1] += _absu(dB[1]) + _absu(G[1][0]) * e[0] + _absu(G[1][1]) * e[1]
        e = [arb(en[0].upper()), arb(en[1].upper())]
        yc, B = yc_new, Bn
        c = arb(nxt)
        s0 = s1
    # deviation set {B z : |z_j| <= e_j}; columns of B have unit 2-norm
    dev = arb((e[0] + e[1]).upper()) * arb("1.0000001")
    if stats is not None:
        stats[which] = {"steps": len(pts), "t0": float(t), "e": [float(e[0].upper()), float(e[1].upper())],
                        "absy": float((abs(yc[0]) + abs(yc[1])).upper())}
    return yc[0], yc[1] * arb(s0), dev, arb(s0)


def wronskian_z(pt: Point, l: int, wc, wbox=None, theta=1.0, frac=0.5, stats=None):
    """Without wbox: (an acb ball containing W_l(w) for every w in the ball wc, None).
    With wbox: (a ball for the centre, widened by the box deviation [not sharp], and an acb ball
    containing W_l(w) for every w in wbox).  Use a separate call for a sharp centre value."""
    L = l * (l + 1)
    wc = acb(wc)
    rm = _rat(pt.rP.mid(), 20)
    Zp, dZp, ep, sm = propagate(pt, "+", wc, wbox, L, theta, frac, rmatch=rm, stats=stats)
    Zc, dZc, ec, _ = propagate(pt, "c", wc, wbox, L, theta, frac, rmatch=rm, stats=stats)
    rma = arb(rm)
    D = acb(pt.Delta(rma))
    den = acb(pt.rp * pt.rc)

    def wform(zp, dzp, zc, dzc, w):
        return (D * (zp * dzc - dzp * zc) + 2 * I * w * acb(rma * rma) * zp * zc) / den

    def widen(z, r):
        return z + acb(arb(0, r.upper()), arb(0, r.upper()))

    # the deviation radii also contain the accumulated rounding of the centre trajectory
    args = (widen(Zp, ep), widen(dZp, ep * sm), widen(Zc, ec), widen(dZc, ec * sm))
    Wc = wform(*args, wc)
    if wbox is None:
        return Wc, None
    return Wc, wform(*args, acb(wbox))
