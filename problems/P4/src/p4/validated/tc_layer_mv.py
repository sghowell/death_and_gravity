"""Theorem C, S4-4: centred (mean-value) enclosures of the layer engine's
Sigma on tube steps, for tc_layer.certify_piece(sigma_fn=...).

Why.  Away from the sonic point at moderate tau (14 <= tau <~ 40) the natural
interval extension of Sigma (tc_layer.sigma_box) loses the correlation between
blocks that are large and nearly cancel: on deep/deep-064 at eps ~ 4.4,
tau ~ 17.5 the block H = V' - V (log B0)' is ~ -85 with Im H ~ 0.04, and Im H
enters Sigma multiplied by 32/tau; a box of width 0.15 in eps, 0.1 in tau and
2e-4 in sigma then has a hull of width ~10 around a true value of -4.6, and
depth 20 of the adaptive cover still fails.  The quantities vary little across
the box, so a centred form, whose excess width is quadratic in the box size,
removes most of that loss.

sigma_box_mv (variables eps, iota = 1/tau).  On a tube step (no Laurent
offsets) Sigma is, for each fixed sigma = x_k - x in the strip and each fixed
choice of the Taylor-model coefficient values, a real C^1 function of
(eps, iota) on the closed box E x I, I = [1/tau_hi, 1/tau_lo] (iota = 0
included: w = 1/kappa = iota/(eps iota + i) is smooth there).  By the mean
value theorem along the segment from the centre c = (eps_c, iota_c),

    Sigma in Sigma(c) + [dSigma/deps](box) (E - eps_c) + [dSigma/diota](box) (I - iota_c),

where Sigma(c) is the natural enclosure at the thin point c with the whole
sigma-strip, and the gradient hulls come from forward-mode differentiation
in ball arithmetic over the whole box (class DV), with
w = iota q, q = 1/(eps iota + i), dw/deps = -w^2, dw/diota = i q^2.

sigma_box_mv3 (variables sigma, eps, iota).  A Taylor model f(sigma) in
sum_k [c_k] sigma^k on [0, R] is not differentiable as an enclosure, but
f = P + delta with P the midpoint polynomial (exact coefficients) and
|delta(sigma)| <= r := sum_k rad(c_k) R^k (sigma >= 0).  Write Sigma =
Phi(f(sigma), eps, iota) with f the vector of all coefficient functions
(those of the kappa-polynomials, tt, tt', u2, and lam = xi = xi_k + sigma).
For fixed delta* = delta(sigma) the map (s, e, i) -> Phi(P(s) + delta*, e, i)
is C^1 on the box, so

    Phi(f(sigma), eps, iota) in Phi(P(sigma_c) + [-r, r], eps_c, iota_c)
        + [dPhi/dsigma](box) (S - sigma_c) + [dPhi/deps](box) (E - eps_c)
        + [dPhi/diota](box) (I - iota_c),

with the coefficient values P(S) + [-r, r] and derivatives P'(S) in the
gradient pass (P' exact from the midpoint coefficients).

Both functions return the intersection with the natural enclosure, so they
are never wider; a denominator hull containing 0 makes the centred form
non-finite and the natural enclosure is returned (non-finite there as well,
so the box is rejected as before).  Laurent (sonic) pieces are passed to
tc_layer.sigma_box unchanged.  This module does not change tc_layer; the
certified statement of a part computed with it is the same: dT, dB, nB != 0
and Sigma <= certified_bound < 0 on the piece.
"""
from __future__ import annotations

from flint import acb, arb

from . import tc_layer as tl
from .tc_tm import TM

INF = float("inf")


def _acb(x):
    return x if isinstance(x, acb) else acb(x)


class DV:
    """A complex quantity and its partial derivatives in n directions:
    v (acb) and d = tuple of n acb.  Forward-mode differentiation with ball
    arithmetic: every operation encloses the true value and derivatives for
    every point of the box the inputs enclose."""

    __slots__ = ("v", "d")

    def __init__(self, v, d=()):
        self.v = _acb(v)
        self.d = tuple(_acb(x) for x in d)

    @staticmethod
    def of(x, n):
        if isinstance(x, DV):
            return x
        return DV(x, (acb(0),) * n)

    def _n(self, o):
        return len(o.d) if isinstance(o, DV) else len(self.d)

    def __add__(self, o):
        o = DV.of(o, len(self.d))
        n = max(len(self.d), len(o.d))
        a, b = self._pad(n), o._pad(n)
        return DV(self.v + o.v, tuple(x + y for x, y in zip(a, b)))

    __radd__ = __add__

    def _pad(self, n):
        return self.d + (acb(0),) * (n - len(self.d))

    def __neg__(self):
        return DV(-self.v, tuple(-x for x in self.d))

    def __sub__(self, o):
        return self + (-DV.of(o, len(self.d)))

    def __rsub__(self, o):
        return DV.of(o, len(self.d)) - self

    def __mul__(self, o):
        if isinstance(o, DV):
            n = max(len(self.d), len(o.d))
            a, b = self._pad(n), o._pad(n)
            return DV(self.v * o.v, tuple(x * o.v + self.v * y for x, y in zip(a, b)))
        o = _acb(o)
        return DV(self.v * o, tuple(x * o for x in self.d))

    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, DV):
            n = max(len(self.d), len(o.d))
            a, b = self._pad(n), o._pad(n)
            q = self.v / o.v
            return DV(q, tuple((x - q * y) / o.v for x, y in zip(a, b)))
        o = _acb(o)
        return DV(self.v / o, tuple(x / o for x in self.d))

    def __rtruediv__(self, o):
        return DV.of(o, len(self.d)) / self

    def __pow__(self, k):
        if not (isinstance(k, int) and k >= 0):
            raise ValueError("DV power: non-negative integers only")
        out = DV(acb(1), (acb(0),) * len(self.d))
        for _ in range(k):
            out = out * self
        return out

    @property
    def real(self):
        return DV(acb(self.v.real), tuple(acb(x.real) for x in self.d))

    @property
    def imag(self):
        return DV(acb(self.v.imag), tuple(acb(x.imag) for x in self.d))


# ------------------------------------------------------------------ algebra
def assemble(V, Vp, Y, Lh, Lp, t, tp, u2h, lam, E, I):
    """Sigma from the normalized blocks: the formula of tc_layer.sigma_box,
    written once for DV (or plain ball) arguments."""
    c = tl.cval()
    W, Wp = V - Lh, Vp - Lp
    G = Y + t * Lh - tp
    H = Vp - V * Lh
    nu, nup = W.imag, Wp.imag
    rho = W.real - (E + c) * t
    cf = c * I * I + 1 / c
    reB1 = u2h * lam * (E * E * I * I - 1) + E * I * I * G.real - I * G.imag + I * I * H.real
    imB1 = -(2 * E * u2h * lam + G.real + E * I * G.imag + I * H.imag) / c
    fp = -cf * tp / 2 + (I / c) * nup / 2
    rhof = -rho * (cf * t - (I / c) * nu) / 2
    nu2 = nu * nu * I * I / 4
    return reB1 + imB1 + fp + rhof + nu2


def _blocks(ld, R):
    """(V, Vp, Y, Lh, Lp) from a ratio evaluator R(num, dens, k)."""
    V = R("nV", ("dT",), 1)
    Y = R("nY", ("kdB",), 2)
    if getattr(ld, "combine", False):
        Vp = R("NVp", ("dT", "dT"), 2)
        Lh = R("NL", ("nB", "dB"), 1)
        Lp = R("Mn", ("nB", "nB"), 2) - R("Md", ("dB", "dB"), 2)
    else:
        Vp = R("nVp", ("dT",), 2) - V * R("dT1", ("dT",), 1)
        dBp, dBpp = R("dB1", ("dB",), 1), R("dB2", ("dB",), 2)
        nBp, nBpp = R("nB1", ("nB",), 1), R("nB2", ("nB",), 2)
        Lh = nBp - dBp
        Lp = (nBpp - nBp * nBp) - (dBpp - dBp * dBp)
    return V, Vp, Y, Lh, Lp


def _qform(F, coef, w):
    """Q-form of a pure w-form, Q = a_d + sum_{j<d} a_j w^{d-j}, from the
    coefficient values coef[j] (DV or balls) at a DV w."""
    if F.p != F.d:
        raise ArithmeticError("centred form: %s is not a pure w-form" % F.name)
    acc = DV.of(acb(0), len(w.d))
    for j in range(0, F.p):
        cj = coef.get(j)
        acc = acc * w + (cj if cj is not None else acb(0))
    cd = coef.get(F.d)
    return acc * w + (cd if cd is not None else acb(0))


def _ratio(N, Ds, k, cN, cDs, lam, w):
    e = sum(D.d for D in Ds) - N.d
    if e < 0:
        raise ArithmeticError("centred form: numerator degree above the denominator's")
    q = _qform(N, cN, w)
    for D, cD in zip(Ds, cDs):
        q = q / _qform(D, cD, w)
    if e:
        q = q * (w ** e)
    return q * (lam ** k) if k else q


def _finite(x):
    return x.mid().is_finite() and x.rad().is_finite()


def _intersect(S_nat, S_mv):
    if not _finite(S_mv):
        return S_nat
    if not _finite(S_nat):
        return S_mv
    lo = max(S_nat.lower(), S_mv.lower())
    hi = min(S_nat.upper(), S_mv.upper())
    if lo > hi:
        raise ArithmeticError("centred and natural enclosures of Sigma are disjoint")
    for S in (S_mv, S_nat):                # one contains the intersection exactly
        if S.lower() >= lo and S.upper() <= hi:
            return S
    return arb(lo).union(arb(hi))          # outward-rounded hull of [lo, hi]


def _centre(g, box):
    """(eps_c, iota_c, thin geometry dict at the centre with the same strip)."""
    tlo, thi = g["tau"]
    ec = (box[2] + box[3]) / 2
    ic = (1.0 / tlo + (0.0 if thi == INF else 1.0 / thi)) / 2
    Tc = arb(1) / arb(ic)
    return ec, ic, Tc


def _w_dv(g, n, ie, ii):
    """w over the box as a DV with d/deps at slot ie and d/diota at slot ii."""
    Eb, itb, wv = g["E"], g["itau"], g["w"]
    q = acb(1) / acb(Eb * itb, arb(1))
    d = [acb(0)] * n
    d[ie] = -(wv * wv)
    d[ii] = acb(0, 1) * q * q
    return DV(wv, d)


# --------------------------------------------------- centred in (eps, iota)
def sigma_box_mv(ld, box, geo=None):
    """Centred form in (eps, iota) intersected with tc_layer.sigma_box (tube
    steps); tc_layer.sigma_box on sonic pieces."""
    g = tl.box_geometry(ld, box) if geo is None else geo
    if g is None:
        return None
    S_nat = tl.sigma_box(ld, box, g)
    if ld.laurent:
        return S_nat
    L, lam, Eb, itb = g["L"], g["lam"], g["E"], g["itau"]
    ec, ic, Tc = _centre(g, box)
    gc = dict(L=L, E=arb(ec), w=acb(1) / acb(arb(ec), Tc), itau=arb(ic), lam=lam,
              z=None, zeta=None, use_large=False)
    try:
        Sc = tl.sigma_box(ld, box, gc)
        F = ld.forms
        w = _w_dv(g, 2, 0, 1)
        strip = lambda nm: F[nm].strip(L)
        coef = lambda nm: dict(list(strip(nm)[0].items()) + list(strip(nm)[1].items()))
        R = lambda n, ds, k: _ratio(F[n], [F[x] for x in ds], k, coef(n),
                                    [coef(x) for x in ds], lam, w)
        V, Vp, Y, Lh, Lp = _blocks(ld, R)
        t, tp, u2h = ld.real_hulls(L, lam)
        Sd = assemble(V, Vp, Y, Lh, Lp, t, tp, u2h, lam,
                      DV(acb(Eb), (acb(1), acb(0))), DV(acb(itb), (acb(0), acb(1))))
        S_mv = Sc + Sd.d[0].real * (Eb - arb(ec)) + Sd.d[1].real * (itb - arb(ic))
    except (ArithmeticError, ZeroDivisionError, ValueError):
        return S_nat
    return _intersect(S_nat, S_mv)


# ---------------------------------------------- centred in (sigma, eps, iota)
class _MidPoly:
    """f = P + delta on [0, R]: P the midpoint polynomial of a TM (exact
    coefficients), P' its derivative, r = sum_k rad(c_k) R^k >= |delta|."""

    __slots__ = ("P", "dP", "r")

    def __init__(self, tm):
        mids = [arb(x.mid()) for x in tm.c]
        self.P = TM(mids, tm.R)
        self.dP = TM([(k + 1) * mids[k + 1] for k in range(len(mids) - 1)] or [arb(0)], tm.R)
        Rb = arb(tm.R)
        r = arb(0)
        for k, x in enumerate(tm.c):
            r += arb(x.rad()) * Rb ** k
        self.r = arb(0, r.upper())

    def val(self, lo, hi):
        return self.P.eval(lo, hi) + self.r

    def der(self, lo, hi):
        return self.dP.eval(lo, hi)


USED_COMBINE = ("nV", "dT", "nY", "kdB", "NVp", "NL", "nB", "dB", "Mn", "Md")
USED_PLAIN = ("nV", "dT", "nY", "kdB", "nVp", "dT1", "dB1", "dB", "dB2", "nB1", "nB", "nB2")


def _mid_data(ld):
    """Midpoint polynomials of every coefficient TM the tube-step Sigma uses
    (cached on the LayerData, with per-strip and per-centre value caches)."""
    md = getattr(ld, "_mv3", None)
    if md is not None:
        return md
    md = dict(forms={}, real={}, strip={}, centre={})
    for nm in (USED_COMBINE if getattr(ld, "combine", False) else USED_PLAIN):
        F = ld.forms[nm]
        if any(F.o[j] != 0 for j in F.idx):
            raise ArithmeticError("centred form: Laurent offsets in %s" % nm)
        md["forms"][nm] = {j: _MidPoly(F.tm[j]) for j in F.idx}
    for nm, lc in (("tt", ld.tt), ("ttp", ld.ttp), ("u2", ld.u2)):
        if lc.off != 0:
            raise ArithmeticError("centred form: Laurent offset in %s" % nm)
        md["real"][nm] = _MidPoly(lc.tm)
    ld._mv3 = md
    return md


def _sigma_dv(ld, md, coef_of, real_of, lam, w, E, I):
    F = ld.forms
    R = lambda n, ds, k: _ratio(F[n], [F[x] for x in ds], k, coef_of(n),
                                [coef_of(x) for x in ds], lam, w)
    V, Vp, Y, Lh, Lp = _blocks(ld, R)
    tt, ttp, u2 = real_of("tt"), real_of("ttp"), real_of("u2")
    t, tp, u2h = tt * lam, ttp * lam * lam, u2 * lam
    return assemble(V, Vp, Y, Lh, Lp, t, tp, u2h, lam, E, I)


def sigma_box_mv3(ld, box, geo=None):
    """Centred form in (sigma, eps, iota) intersected with tc_layer.sigma_box
    (tube steps); tc_layer.sigma_box on sonic pieces."""
    g = tl.box_geometry(ld, box) if geo is None else geo
    if g is None:
        return None
    S_nat = tl.sigma_box(ld, box, g)
    if ld.laurent:
        return S_nat
    try:
        md = _mid_data(ld)
        L, Eb, itb = g["L"], g["E"], g["itau"]
        lo, hi = L
        sc = (lo + hi) / 2
        ec, ic, Tc = _centre(g, box)
        # centre: coefficient balls P(sigma_c) + [-r, r], thin (eps_c, iota_c)
        z0 = ()
        if sc not in md["centre"]:
            md["centre"][sc] = (
                {nm: {j: DV(mp.val(sc, sc), z0) for j, mp in fm.items()}
                 for nm, fm in md["forms"].items()},
                {nm: DV(mp.val(sc, sc), z0) for nm, mp in md["real"].items()})
        cc, rc = md["centre"][sc]
        lam_c = DV(acb(ld.xik + arb(sc)), z0)
        wc = DV(acb(1) / acb(arb(ec), Tc), z0)
        Sc = _sigma_dv(ld, md, lambda n: cc[n], lambda n: rc[n], lam_c, wc,
                       DV(acb(arb(ec)), z0), DV(acb(arb(ic)), z0)).v.real
        # gradient over the box: slots (sigma, eps, iota)
        if (lo, hi) not in md["strip"]:
            md["strip"][(lo, hi)] = (
                {nm: {j: DV(mp.val(lo, hi), (acb(mp.der(lo, hi)), acb(0), acb(0)))
                      for j, mp in fm.items()} for nm, fm in md["forms"].items()},
                {nm: DV(mp.val(lo, hi), (acb(mp.der(lo, hi)), acb(0), acb(0)))
                 for nm, mp in md["real"].items()})
        cs, rs = md["strip"][(lo, hi)]
        lam = DV(acb(g["lam"]), (acb(1), acb(0), acb(0)))
        w = _w_dv(g, 3, 1, 2)
        Sd = _sigma_dv(ld, md, lambda n: cs[n], lambda n: rs[n], lam, w,
                       DV(acb(Eb), (acb(0), acb(1), acb(0))),
                       DV(acb(itb), (acb(0), acb(0), acb(1))))
        Sb = arb(lo).union(arb(hi))
        S_mv = (Sc + Sd.d[0].real * (Sb - arb(sc)) + Sd.d[1].real * (Eb - arb(ec))
                + Sd.d[2].real * (itb - arb(ic)))
    except (ArithmeticError, ZeroDivisionError, ValueError):
        return S_nat
    return _intersect(S_nat, S_mv)
