"""Theorem C, S4-6: the mean-value (centred) form of the centre sigma.

tc_layer_centre.sigma_box_centre encloses Sigma_c = S/(1 + (y tau)^2) by
graded natural interval evaluation.  On wide boxes its hulls carry the box's
ranges independently through many terms (the dependency problem).  This
module adds the centred form, as tc_layer_mv's mv3 does for tube steps:

    Sigma_c(v) in Sigma_c(c) + sum_i [dSigma_c/dv_i](box) (V_i - c_i),

over three variables v = (g, h, eps) per regime of tc_layer_centre:

    regime L:  g = y,  h = sigma = 1/s,  eps;
    regime S:  g = r = 1/tau,  h = s,  eps.

Taylor-model coefficients.  A TM with ball coefficients [m_k +/- r_k] on
[0, R] encloses f pointwise; so f(y) = P(y) + delta(y) with P the exact
midpoint polynomial and |delta(y)| <= rho = sum r_k y_hi^k for
0 <= y <= y_hi.  For each fixed delta* = delta(y(v)) (frozen) the map
v -> Phi(v, P(y(v)) + delta*) is C^1 on the box, so the mean value theorem
gives the form above with the centre value computed from P(y_c) + [-rho, rho]
and the gradient hulls from P(Y) + [-rho, rho] (values) and P'(Y)
(derivatives), Y covering y(v) on the box (tc_layer_mv, mv3).  In regime S,
y = g h, so Y = [g_lo h_lo, g_hi h_hi]; segments may leave the y-strip, which
is harmless (P is a polynomial and delta* is frozen; the enclosure is used
only at true points, where y lies in the strip).

Gradients.  Forward-mode differentiation in ball arithmetic (class DV: value
and the three partials), carried through the graded arithmetic of
tc_layer_centre (class GVD: integer exponents of g with DV coefficients).
The total g-derivative of sum_a g^a c_a is sum_a (a g^{a-1} c_a + g^a dc_a/dg);
the coefficients' own g-dependence (the TM argument, D, Dr) is in their DV
partials, and folding g^a c into the bucket g^{K+1} multiplies by the DV of
g^{a-K-1}, so the bucket's partials keep that dependence.  The Neumann
remainder (-u)^{K+2}/(1 + u) of a division uses the DV of 1/(1 + u) built
from u's value and gradient hulls.  Negative exponents (lam = g^{-1}) must
end with exactly zero coefficients, value and partials (as in
tc_layer_centre; checked).

The centred form is intersected with the natural enclosure of the same
regime (never wider).

Status: correct (point values lie in its enclosures on random boxes, tests)
but not used by the production driver.  On c29 its gradient hulls carry the
same dependency as the values, so it is only sometimes tighter than the
natural enclosure, at about six times the cost; the natural enclosure on the
identity-raised layer data certified c29 at tau >= 14 first
(notes/s4-6-centre.md).  Kept as an optional evaluator (sigma_fn).
"""
from __future__ import annotations

from flint import acb, arb

from . import tc_layer as tl
from . import tc_layer_centre as tcc
from .tc_tm import TM

INF = float("inf")
_Z = acb(0)


# =====================================================================
# forward-mode values
# =====================================================================
class DV:
    """value v and partials d = (d/dg, d/dh, d/deps), complex balls."""
    __slots__ = ("v", "d")

    def __init__(self, v, d=None):
        self.v = v if isinstance(v, acb) else acb(v)
        self.d = d if d is not None else (_Z, _Z, _Z)

    def __add__(self, o):
        if not isinstance(o, DV):
            return DV(self.v + acb(o), self.d)
        return DV(self.v + o.v, (self.d[0] + o.d[0], self.d[1] + o.d[1], self.d[2] + o.d[2]))

    __radd__ = __add__

    def __neg__(self):
        return DV(-self.v, (-self.d[0], -self.d[1], -self.d[2]))

    def __sub__(self, o):
        return self + (-o if isinstance(o, DV) else -acb(o))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, DV):
            o = acb(o)
            return DV(self.v * o, (self.d[0] * o, self.d[1] * o, self.d[2] * o))
        v = self.v * o.v
        return DV(v, (self.d[0] * o.v + self.v * o.d[0], self.d[1] * o.v + self.v * o.d[1],
                      self.d[2] * o.v + self.v * o.d[2]))

    __rmul__ = __mul__

    def inv(self):
        if self.v.contains(0):
            raise ArithmeticError("DV inverse: value contains 0")
        iv = acb(1) / self.v
        q = -(iv * iv)
        return DV(iv, (self.d[0] * q, self.d[1] * q, self.d[2] * q))

    @property
    def real(self):
        return DV(acb(self.v.real), tuple(acb(x.real) for x in self.d))

    @property
    def imag(self):
        return DV(acb(self.v.imag), tuple(acb(x.imag) for x in self.d))

    def is_zero(self):
        return self.v.is_zero() and all(x.is_zero() for x in self.d)


def _dv(v, d=None):
    return v if isinstance(v, DV) else DV(v, d)


# =====================================================================
# graded values with DV coefficients
# =====================================================================
class GVD:
    """sum_a g^a c[a] with DV coefficients; g in the hull gh >= 0."""
    __slots__ = ("c", "gh")

    def __init__(self, c, gh):
        self.c, self.gh = c, gh

    @classmethod
    def const(cls, v, gh):
        return cls({0: _dv(v)}, gh)

    @classmethod
    def mono(cls, e, v, gh):
        return cls({e: _dv(v)}, gh)

    def _gpow(self, k):
        """DV of g^k (k >= 1) on the hull: value gh^k, d/dg k gh^(k-1)."""
        return DV(acb(self.gh ** k), (acb(k * self.gh ** (k - 1)), _Z, _Z))

    def _put(self, e, v, out):
        if e > tcc.KMAX + 1:
            v = v * self._gpow(e - tcc.KMAX - 1)
            e = tcc.KMAX + 1
        out[e] = out[e] + v if e in out else v

    def __add__(self, o):
        if not isinstance(o, GVD):
            o = GVD.const(o, self.gh)
        out = dict(self.c)
        for e, v in o.c.items():
            out[e] = out[e] + v if e in out else v
        return GVD(out, self.gh)

    __radd__ = __add__

    def __neg__(self):
        return GVD({e: -v for e, v in self.c.items()}, self.gh)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, GVD):
            return GVD({e: c * o for e, c in self.c.items()}, self.gh)
        out = {}
        for ea, ca in self.c.items():
            for eb, cb in o.c.items():
                self._put(ea + eb, ca * cb, out)
        return GVD(out, self.gh)

    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, GVD):
            return self * o.inv()
        return self * (acb(1) / acb(o))

    def __pow__(self, n):
        r = GVD.const(1, self.gh)
        for _ in range(n):
            r = r * self
        return r

    def shift(self, k):
        out = {}
        for e, v in self.c.items():
            self._put(e + k, v, out)
        return GVD(out, self.gh)

    @property
    def real(self):
        return GVD({e: v.real for e, v in self.c.items()}, self.gh)

    @property
    def imag(self):
        return GVD({e: v.imag for e, v in self.c.items()}, self.gh)

    def total(self):
        """(value hull, gradient hulls) of sum_a g^a c[a] over the box.
        Negative exponents need exactly-zero coefficients (value and
        partials)."""
        gh = self.gh
        val = acb(0)
        dg, dh, de = acb(0), acb(0), acb(0)
        for e, c in self.c.items():
            if e < 0:
                if c.is_zero():
                    continue
                raise ArithmeticError("graded DV: nonzero coefficient at g^%d" % e)
            if e == 0:
                val += c.v
                dg += c.d[0]
                dh += c.d[1]
                de += c.d[2]
            else:
                ge = gh ** e
                val += ge * c.v
                dg += ge * c.d[0] + (e * gh ** (e - 1)) * c.v
                dh += ge * c.d[1]
                de += ge * c.d[2]
        return val, (dg, dh, de)

    def inv(self):
        if any(e < 0 and not v.is_zero() for e, v in self.c.items()):
            raise ArithmeticError("graded DV inverse: negative exponent")
        q0 = self.c.get(0)
        if q0 is None or q0.v.contains(0):
            raise ArithmeticError("graded DV inverse: leading coefficient contains 0")
        iq = q0.inv()
        u = GVD({e: v * iq for e, v in self.c.items() if e > 0}, self.gh)
        if not u.c:
            return GVD.const(iq, self.gh)
        uv, ug = u.total()
        U = DV(uv, ug)
        one_u = U + 1
        if one_u.v.contains(0):
            raise ArithmeticError("graded DV inverse: 1 + u contains 0")
        acc = GVD.const(1, self.gh)
        term = GVD.const(1, self.gh)
        mu = -u
        for _ in range(tcc.KMAX + 1):
            term = term * mu
            acc = acc + term
        rem = (term * mu) * one_u.inv()
        return (acc + rem) * iq


# =====================================================================
# the graded variables of one evaluation (box or centre point)
# =====================================================================
class TMSplit:
    """The midpoint polynomial P of a TM and the radius bound rho of its
    ball part on [0, y_hi]: f = P + delta, |delta| <= rho."""
    __slots__ = ("P", "dP", "rad")

    def __init__(self, tm):
        mids = [arb(c.mid()) for c in tm.c]
        self.P = TM(mids, tm.R, tm.ve)
        self.dP = TM([mids[k] * k for k in range(1, len(mids))] or [arb(0)], tm.R, tm.ve)
        self.rad = [arb(c.rad()) for c in tm.c]

    def rho(self, y_hi):
        yb = arb(y_hi)
        acc = arb(0)
        for k, r in enumerate(self.rad):
            if not r.is_zero():
                acc += r * yb ** k
        return arb(0, acc.upper())


def _split(ld, key, tm, cache):
    k = ("split",) + key
    if k not in cache:
        cache[k] = TMSplit(tm)
    return cache[k]


class GradeD:
    """Graded DV variables of one evaluation.  hull: the box (gradient pass)
    or the centre point (centre pass); y_hi bounds delta over the true box."""

    def __init__(self, regime, gh, hh, eh, yh, y_hi, centre):
        self.regime, self.gh, self.centre = regime, gh, centre
        self.y_hi = y_hi
        self.yh = yh                          # hull of the TM argument y
        one = GVD.const(1, gh)
        self.one = one
        E = DV(acb(eh), (_Z, _Z, acb(1)))
        self.E = E
        if regime == "L":
            sg = DV(acb(hh), (_Z, acb(1), _Z))                  # sigma = 1/s
            self.sigma = sg
            # D = 1 + (g eps sigma)^2  (g = y)
            ges = DV(acb(gh * eh * hh), (acb(eh * hh), acb(gh * eh), acb(gh * hh)))
            D = ges * ges + 1
            iD = D.inv()
            self.ZETA = GVD({0: DV(acb(0, 1)) * (-1) * sg * iD, 1: E * sg * sg * iD}, gh)
            self.W = self.ZETA.shift(1)
            self.IT = GVD.mono(1, sg, gh)
            self.Y = GVD.mono(1, 1, gh)
            s2 = sg * sg
            self.norm = (s2 + 1).inv()                           # s^2/(1 + s^2)
        else:
            s = DV(acb(hh), (_Z, acb(1), _Z))
            self.s = s
            ge = DV(acb(gh * eh), (acb(eh), _Z, acb(gh)))        # r eps
            Dr = ge * ge + 1
            iDr = Dr.inv()
            self.W = GVD({1: DV(acb(0, -1)) * iDr, 2: E * iDr}, gh)
            self.ZETA = None
            self.IT = GVD.mono(1, 1, gh)
            self.Y = GVD.mono(1, s, gh)                          # y = r s
            self.norm = (s * s + 1).inv()
        self.LAM = GVD.mono(-1, 1, gh)

    def ypow(self, e):
        if e == 0:
            return self.one
        if e > 0:
            return self.Y ** e
        if self.regime == "L":
            return GVD.mono(e, 1, self.gh)
        sp = self.s
        q = GVD.const(1, self.gh)
        for _ in range(-e):
            q = q * sp.inv()
        return q.shift(e)

    def tmval(self, sp):
        """DV of a TM coefficient f(y) = P(y) + delta*: value P(Y) + [-rho, rho],
        derivatives through P'(Y) and y(v)."""
        lo, hi = self.yh
        Pv = sp.P.eval(lo, hi) + sp.rho(self.y_hi)
        dP = sp.dP.eval(lo, hi)
        if self.regime == "L":
            return DV(acb(Pv), (acb(dP), _Z, _Z))
        return DV(acb(Pv), (acb(dP * self.s.v.real), acb(dP * self.gh), _Z))


def _form_small(F, gd, cache, key):
    k = ("form", key, F.name)
    if k in cache:
        return cache[k]
    acc = None
    for j in F.idx:
        h = gd.tmval(_split(None, (F.name, j), F.tm[j], cache))
        if j >= F.p:
            t = gd.ypow(F.e[j]) * h
            if j > F.p:
                raise ArithmeticError("centre MV form: z-power in %s" % F.name)
        else:
            t = gd.ypow(F.g[j]) * h * gd.W ** (F.p - j)
        acc = t if acc is None else acc + t
    cache[k] = acc
    return acc


def _ratio(N, D, k, gd, cache, key):
    """lam^k N/D for pure w-forms (pivot = top, the raised centre data)."""
    q = _form_small(N, gd, cache, key) / _form_small(D, gd, cache, key)
    zp = N.p - D.p
    ye = N.m - D.m
    if zp < 0:
        q = q * gd.W ** (-zp)
        ye += zp
    elif zp > 0:
        raise ArithmeticError("centre MV form: positive z-power %s/%s" % (N.name, D.name))
    return q * gd.ypow(ye) * (gd.LAM ** k if k else gd.one)


def _real(lc, name, k, gd, cache):
    h = gd.tmval(_split(None, ("real", name), lc.tm, cache))
    return gd.ypow(lc.off) * h * (gd.LAM ** k if k else gd.one)


def sigma_graded_dv(ld, gd, cache, key):
    """Sigma = S lam^2 it^2 as a GVD (the formula of tc_layer.sigma_box)."""
    F = ld.forms
    cent = getattr(ld, "centring", {})

    def R(n, d, k):
        if (n, d) in cent:
            c = cent[(n, d)]
            r = _ratio(F["C_%s_%s" % (n, d)], F[d], k, gd, cache, key)
            if c:
                r = r + (gd.LAM ** k if k else gd.one) * arb(c)
            return r
        return _ratio(F[n], F[d], k, gd, cache, key)
    c = tl.cval()
    ic = arb(1) / c
    t = _real(ld.tt, "tt", 1, gd, cache)
    tp = _real(ld.ttp, "ttp", 2, gd, cache)
    u2h = _real(ld.u2, "u2", 1, gd, cache)
    lam, it = gd.LAM, gd.IT
    Eg = GVD.const(gd.E, gd.gh)
    V = R("nV", "dT", 1)
    Y = R("nY", "kdB", 2)
    Vp = R("nVp", "dT", 2) - V * R("dT1", "dT", 1)
    dBp, dBpp = R("dB1", "dB", 1), R("dB2", "dB", 2)
    nBp, nBpp = R("nB1", "nB", 1), R("nB2", "nB", 2)
    Lh = nBp - dBp
    Lp = (nBpp - nBp * nBp) - (dBpp - dBp * dBp)
    W, Wp = V - Lh, Vp - Lp
    G = Y + t * Lh - tp
    H = Vp - V * Lh
    nu, nup = W.imag, Wp.imag
    rho = W.real - (Eg + c) * t
    cf = c * it * it + ic
    return (u2h * lam * (Eg * Eg * it * it - 1) + Eg * it * it * G.real
            - it * G.imag + it * it * H.real
            - (2 * Eg * u2h * lam + G.real + Eg * it * G.imag + it * H.imag) * ic
            - cf * tp * arb(0.5) + (it * ic) * nup * arb(0.5)
            - rho * (cf * t - (it * ic) * nu) * arb(0.5)
            + nu * nu * it * it * arb(0.25))


def _variables(regime, g):
    """(g-hull, h-hull, eps-hull, y-hull for the TM argument, y_hi) of the
    box from tl.box_geometry's effective ranges."""
    L, E = g["L"], g["E"]
    sl, sh = g["s"]
    tlo, thi = g["tau"]
    if regime == "L":
        if not sl > 0:
            raise ArithmeticError("regime L needs s > 0")
        gh = tl.iv(L[0], L[1])
        hh = (arb(0) if sh == INF else arb(1) / arb(sh)).union(arb(1) / arb(sl))
        return gh, hh, E, (L[0], L[1]), L[1]
    if sh == INF:
        raise ArithmeticError("regime S needs bounded s")
    gh = (arb(0) if thi == INF else arb(1) / arb(thi)).union(arb(1) / arb(tlo))
    hh = tl.iv(sl, sh)
    yl = float((gh * arb(sl)).lower()) * (1 - 1e-14)
    yu = float((gh * arb(sh)).upper()) * (1 + 1e-14)
    return gh, hh, E, (max(yl, 0.0), yu), L[1]


def _mid_point(x):
    return arb(x.mid())


def sigma_mv(ld, box, geo, regime):
    """Centred form of Sigma_c on the box in one regime (a ball), or raises."""
    gh, hh, eh, yh, y_hi = _variables(regime, geo)
    cache = {}
    # gradient pass over the box
    gdb = GradeD(regime, gh, hh, eh, yh, y_hi, False)
    tot = sigma_graded_dv(ld, gdb, cache, "box") * gdb.norm
    val, grad = tot.total()
    # centre pass at the exact midpoint (delta widened to the box's rho)
    gc, hc, ec = _mid_point(gh), _mid_point(hh), _mid_point(eh)
    yc = gc if regime == "L" else gc * hc
    ycf = float(yc.mid())                 # bracket the exact centre y for TM.eval
    gdc = GradeD(regime, gc, hc, ec, (ycf * (1 - 1e-14), ycf * (1 + 1e-14)), y_hi, True)
    totc = sigma_graded_dv(ld, gdc, cache, "centre") * gdc.norm
    cv, _ = totc.total()
    out = cv.real
    for dvar, rng, cen in ((grad[0], gh, gc), (grad[1], hh, hc), (grad[2], eh, ec)):
        out = out + dvar.real * (rng - cen)
    nat = val.real
    if not out.overlaps(nat):
        raise ArithmeticError("centre MV: centred and natural enclosures disjoint")
    return tcc._intersect(out, nat)


def sigma_box_centre_mv(ld, box, geo=None):
    """Ball enclosing Sigma_c = S/(1 + (y tau)^2) over the box: the centred
    form in every applicable regime, intersected with each other and with
    the regime's natural enclosure (computed in the same pass).  Plugs into
    tc_layer.certify_piece / run_layer_piece as sigma_fn."""
    g = tl.box_geometry(ld, box) if geo is None else geo
    if g is None:
        return None
    out, err = None, None
    for rg in tcc.regimes_of(g["s"]):
        try:
            v = sigma_mv(ld, box, g, rg)
        except (ArithmeticError, ZeroDivisionError, ValueError) as ex:
            err = ex
            continue
        if out is None:
            out = v
        else:
            if not out.overlaps(v):
                raise ArithmeticError("disjoint regime enclosures (a bug)")
            out = tcc._intersect(out, v)
    if out is None:
        raise err
    return out
