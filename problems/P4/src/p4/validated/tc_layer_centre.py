"""Theorem C, S4-6: the centre variant of the layer engine (probe stage).

The centre region is y = e^x in (0, e^{-4.5}], 30 pieces in y (the plan's
centre/c00..c29; c00 = [0, y0] closes at y = 0).  The commuted energy
condition is the one of tc_layer (notes/s4-commutator.md section 3):

    S = Re B1 - (tau/c) Im B1 - (f' + rho f)/2 + nu^2/4 < 0,

with the same pairs (tc_commutator.build_pairs over tc_tm.centre_components,
built from the certified (a*, mu*) balls of ec_centre).  Near the centre S
scales in the layer variable s = y tau, not in xi = -x:

    S = G(s, eps) + O(y),   G(0+, eps) = -6,   G(s, eps) ~ k2(eps) s^2 (s large),

k2(0) = -19.0, k2(1) = -144.5, k2(6) = -772 (float data, notes/s4-6-centre.md).
So the quantity certified here is

    Sigma_c := S / (1 + s^2),          s = y tau,

which stays O(1) and bounded away from 0 as y -> 0 and as tau -> infinity.

Pivot forms.  Every kappa-polynomial P = sum_j y^{o_j} TM_j(y) kappa^j of the
layer data is the tc_layer.PForm in the Laurent variable y, with
z = y kappa and w = 1/kappa:

    P = y^m z^p [ sum_{j>=p} TM_j y^{e_j} z^{j-p} + sum_{j<p} TM_j y^{g_j} w^{p-j} ].

tc_layer.LayerData builds the blocks (nV, dT, nVp, dT1, nY, kdB, dB, dB1, dB2,
nB, nB1, nB2, tt, tt', u2) and checks that every denominator has its pivot
and its top on the layer line; it is used unchanged (laurent=True with the
TM variable ell = y).

Graded evaluation.  The sonic engine (tc_layer.sigma_box) normalizes with
lam = ell, and every block's ell-power cancels exactly.  At the centre no
single power of y does that: the blocks V, Lh, ... are O(1), while tt, tt'
and u2 are O(y^2) and the t-linear terms of S are O(1) only through tau y = s.
So each quantity is carried as a polynomial in one grading variable g with
exactly tracked integer exponents and ball coefficients (class GV), and the
y-powers cancel by exponent bookkeeping, not by interval arithmetic:

  regime L (s >= S_SPLIT, any y including y -> 0):  g = y,
      it = 1/tau = g / s,  z = i s + g eps,  w = 1/kappa = g zeta,
      zeta = 1/z = (g eps - i s)/(s^2 + y^2 eps^2),  lam = g^{-1};
      the tc_layer identity S lam^2 it^2 = T(lam) gives Sigma = S/s^2.
  regime S (s <= S_SPLIT):  g = r = 1/tau,  y = g s,
      it = g,  z = i s + g s eps,  w = g (g eps - i)/(1 + r^2 eps^2),  lam = g^{-1};
      then Sigma = S.

All are exact identities (z w = y, tau y = s); g -> 0 is the corner y -> 0
(L) or tau -> infinity at bounded s (S).  The exponent-0 coefficients of the
pivot forms are real (only even powers of z appear without a y factor, and
z = i s at g = 0), so the negative exponents created by lam cancel exactly:
GV.value() refuses a box if any negative-exponent coefficient of the final
Sigma is not exactly zero.  Exponents above KMAX are folded into KMAX + 1
with the g-hull (sum_{a > K} g^a c_a = g^{K+1} sum g^{a-K-1} c_a, contained
in g^{K+1} * sum gh^{a-K-1} c_a).  Division by a pivot-form denominator
Q = q0 (1 + u), u = O(g): the Neumann series to order KMAX + 1 plus the exact
remainder (-u)^{KMAX+2}/(1 + u), bounded with the hull of u; a box whose q0
or 1 + u hull contains 0 is rejected (so dT, dB, nB != 0 on every certified
box with y > 0).  Finally Sigma_c = Sigma * s^2/(1 + s^2) (L) or
Sigma/(1 + s^2) (S).

The cover is tc_layer.certify_piece / run_layer_piece, unchanged, with
sigma_fn = sigma_box_centre and boxes (y-strip) x (eps-cell) x (s = y tau) x
(tau): tc_layer.box_geometry with laurent=True reads ell as y and xi(ell) =
ell, so its s is y tau.

Not done at this stage (see notes/s4-6-centre.md): a centred (mean-value)
form for Laurent coefficients, and an exact identity for the near-zero
leading Taylor coefficients of nB, nB1, nB2 and nY (balls around 0 of radius
~1e3 against coefficients ~1e14-1e15, from the width of the certified
(a*, mu*) balls); they dominate the pivot forms of the closing piece when
s = y tau is below ~1e-5.
"""
from __future__ import annotations

from flint import acb, arb

from . import tc_layer as tl
from .arbseries import precision

KMAX = 4              # exponents tracked exactly up to KMAX; higher folded into KMAX + 1
S_SPLIT = 0.5         # regime S below, regime L above (in s = y tau)
Z_LARGE = 8.0         # regime L: small-z forms for s <= Z_LARGE, large-z (zeta) forms above
INF = float("inf")


# =====================================================================
# graded values
# =====================================================================
class GV:
    """sum_a g^a c[a]: integer exponents a (possibly negative), complex ball
    coefficients, g a real variable in the hull gh (gh >= 0).  Exponents
    above KMAX are folded into KMAX + 1 (exactly, by containment)."""
    __slots__ = ("c", "gh")

    def __init__(self, c, gh):
        self.c, self.gh = c, gh

    @classmethod
    def const(cls, v, gh):
        return cls({0: acb(v)}, gh)

    @classmethod
    def mono(cls, e, v, gh):
        return cls({e: acb(v)}, gh)

    def _fold(self, e, v, out):
        if e > KMAX + 1:
            v = v * self.gh ** (e - KMAX - 1)
            e = KMAX + 1
        out[e] = out[e] + v if e in out else v

    def __add__(self, o):
        if not isinstance(o, GV):
            o = GV.const(o, self.gh)
        out = dict(self.c)
        for e, v in o.c.items():
            out[e] = out[e] + v if e in out else v
        return GV(out, self.gh)

    __radd__ = __add__

    def __neg__(self):
        return GV({e: -v for e, v in self.c.items()}, self.gh)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, GV):
            v = acb(o)
            return GV({e: c * v for e, c in self.c.items()}, self.gh)
        out = {}
        for ea, ca in self.c.items():
            for eb, cb in o.c.items():
                self._fold(ea + eb, ca * cb, out)
        return GV(out, self.gh)

    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, GV):
            return self * o.inv()
        v = acb(1) / acb(o)
        return self * v

    def __pow__(self, n):
        r = GV.const(1, self.gh)
        for _ in range(n):
            r = r * self
        return r

    def shift(self, k):
        """g^k * self (exact)."""
        out = {}
        for e, v in self.c.items():
            self._fold(e + k, v, out)
        return GV(out, self.gh)

    @property
    def real(self):
        return GV({e: acb(v.real) for e, v in self.c.items()}, self.gh)

    @property
    def imag(self):
        return GV({e: acb(v.imag) for e, v in self.c.items()}, self.gh)

    def min_exp(self):
        nz = [e for e, v in self.c.items() if not v.is_zero()]
        return min(nz) if nz else None

    def value(self, allow_negative=False):
        """Ball containing sum_a g^a c[a] for every g in gh.  Negative
        exponents must carry exactly-zero coefficients (unless
        allow_negative, which needs gh bounded away from 0)."""
        acc = acb(0)
        for e, v in self.c.items():
            if e < 0:
                if v.is_zero():
                    continue
                if not allow_negative:
                    raise ArithmeticError("graded value: nonzero coefficient at g^%d" % e)
                if not (self.gh > 0):
                    raise ArithmeticError("graded value: g^%d with g-hull touching 0" % e)
            acc = acc + (v * self.gh ** e if e else v)
        return acc

    def inv(self):
        """1/self for self = q0 + O(g), q0 = c[0] != 0, no negative exponents:
        (1/q0) [sum_{k <= K+1} (-u)^k + (-u)^{K+2}/(1 + u)], u = self/q0 - 1."""
        if any(e < 0 and not v.is_zero() for e, v in self.c.items()):
            raise ArithmeticError("graded inverse: negative exponent")
        q0 = self.c.get(0)
        if q0 is None or q0.contains(0):
            raise ArithmeticError("graded inverse: leading coefficient contains 0")
        iq = acb(1) / q0
        u = GV({e: v * iq for e, v in self.c.items() if e > 0}, self.gh)
        if not u.c:
            return GV.const(iq, self.gh)
        uv = u.value()
        one_u = acb(1) + uv
        if one_u.contains(0):
            raise ArithmeticError("graded inverse: 1 + u contains 0")
        acc = GV.const(1, self.gh)
        term = GV.const(1, self.gh)
        mu = -u
        for _ in range(KMAX + 1):
            term = term * mu
            acc = acc + term
        rem = (term * mu) * (acb(1) / one_u)
        return (acc + rem) * iq


# =====================================================================
# layer data at the centre
# =====================================================================
class CentreLayerData(tl.LayerData):
    """tc_layer.LayerData of the centre pairs (Laurent variable y = e^x, TM
    domain [0, hi]).  The construction and its checks are tc_layer's
    (laurent=True); the derivative raise of tc_layer for (tt dT_k)' does not
    fire on the centre pairs (checked: the computed offsets are already at
    the target), and it would be sound here too, since d/dx = y d/dy keeps
    Laurent offsets."""

    centre = True
    # (numerator, denominator) pairs evaluated as c + (N - c D)/D: at the
    # centre N ~ c D to leading order (log-derivatives of y-homogeneous
    # leading parts: nB ~ y^10, dB ~ y^8, dT ~ y^7, and V = nV/dT -> 3), so
    # the difference N - c D, formed coefficientwise on the Taylor models
    # before any hull is taken, is small and carries the dependence on
    # (y, s) that the separate hulls of N and D would not cancel.
    CENTRED = (("nV", "dT"), ("nVp", "dT"), ("dT1", "dT"), ("dB1", "dB"),
               ("dB2", "dB"), ("nB1", "nB"), ("nB2", "nB"))

    def __init__(self, pr, label=""):
        super().__init__(pr, True, label=label)
        from .tc_tm import LC
        self.centring = {}
        for nn, dn in self.CENTRED:
            N, D = self.forms[nn], self.forms[dn]
            c = centring_constant(N, D)
            kp = []
            for j in range(max(N.d, D.d) + 1):
                a = LC(N.tm[j], N.o[j]) if j in N.tm else None
                b = LC(D.tm[j], D.o[j]) * arb(c) if j in D.tm else None
                if a is None and b is None:
                    kp.append(tl._lc_zero_like(LC(D.tm[D.d], 0)))
                elif b is None:
                    kp.append(a)
                elif a is None:
                    kp.append(-b)
                else:
                    kp.append(a - b)
            self.forms["C_%s_%s" % (nn, dn)] = tl.PForm(kp, "C_%s_%s" % (nn, dn))
            self.centring[(nn, dn)] = c


def centring_constant(N, D):
    """A fixed exact constant c ~ N/D at the centre: the ratio of the top
    coefficients' constant terms, rounded to an integer when within 1e-6 of
    one, else to 2^-20.  Any constant gives an exact identity; this one makes
    N - c D small."""
    a = float(N.tm[N.d].c[0].mid())
    b = float(D.tm[D.d].c[0].mid())
    if N.d != D.d or b == 0.0:
        return 0
    r = a / b
    if abs(r - round(r)) < 1e-6 * max(1.0, abs(r)):
        return int(round(r))
    return round(r * 2 ** 20) / 2 ** 20


def centre_layer_data(at, ce, hi, prec=320, label=None):
    """CentreLayerData of a centre piece [lo, hi]: centre_components(ce, hi)
    on y in [0, hi], built exactly as tc_commutator.centre_piece builds its
    pairs (AtomKP with sd_shift = 0)."""
    from . import tc_commutator as tc
    from .tc_tm import AtomKP, centre_components
    with precision(prec):
        comps, _ = centre_components(ce, hi)
        ak = AtomKP(at, comps, sd_shift=0)
        pr = tc.build_pairs(ak)
        return CentreLayerData(pr, label=label or "centre[0,%g]" % hi)


# =====================================================================
# the graded box variables
# =====================================================================
class Grade:
    """The graded variables of one box: g-hull, Y (y as a GV), Z (z = y kappa),
    W (w = 1/kappa), ZETA (1/z, regime L only), IT (1/tau), LAM (lambda) and
    the final normalization factor (a real ball) taking Sigma to
    Sigma_c = S/(1 + s^2)."""

    def __init__(self, regime, L, E, sb, tb):
        self.regime = regime
        e = E
        if regime == "L":
            ylo, yhi = L
            gh = tl.iv(ylo, yhi)
            sl, sh = sb
            if not sl > 0:
                raise ArithmeticError("regime L needs s > 0")
            # sigma = 1/s on [1/sh, 1/sl] (sh may be inf)
            isg = (arb(0) if sh == INF else arb(1) / arb(sh)).union(arb(1) / arb(sl))
            self.gh = gh
            one = GV.const(1, gh)
            self.Y = GV.mono(1, 1, gh)
            self.IT = GV.mono(1, isg, gh)
            D = (arb(1) + (gh * e * isg) ** 2)          # (s^2 + y^2 eps^2)/s^2, real
            # zeta = 1/z = (y eps - i s)/(s^2 + y^2 eps^2) = sigma (-i + y eps sigma)/D
            self.ZETA = GV({0: acb(0, -isg / D), 1: acb(e * isg * isg / D)}, gh)
            self.W = self.ZETA.shift(1)                 # w = y zeta
            if sh <= Z_LARGE:
                sbox = tl.iv(sl, sh)
                self.Z = GV({0: acb(0, sbox), 1: acb(e)}, gh)
            else:
                self.Z = None
            self.LAM = GV.mono(-1, 1, gh)
            self.norm = arb(1) / (arb(1) + isg * isg)   # s^2/(1 + s^2)
            self.one = one
            self.sbox = None if sh == INF else tl.iv(sl, sh)
        elif regime == "S":
            tlo, thi = tb
            rg = (arb(0) if thi == INF else arb(1) / arb(thi)).union(arb(1) / arb(tlo))
            gh = rg
            sl, sh = sb
            if sh == INF:
                raise ArithmeticError("regime S needs bounded s")
            sbox = tl.iv(sl, sh)
            self.gh = gh
            self.one = GV.const(1, gh)
            self.Y = GV.mono(1, sbox, gh)                # y = r s
            self.IT = GV.mono(1, 1, gh)                  # 1/tau = r
            Dr = arb(1) + (rg * e) ** 2
            # w = 1/kappa = r/(i + r eps) = r (r eps - i)/(1 + r^2 eps^2)
            self.W = GV({1: acb(0, -arb(1) / Dr), 2: acb(e / Dr)}, gh)
            self.Z = GV({0: acb(0, sbox), 1: acb(sbox * e)}, gh)
            self.ZETA = None
            self.LAM = GV.mono(-1, 1, gh)
            self.norm = arb(1) / (arb(1) + sbox * sbox)  # 1/(1 + s^2)
            self.sbox = sbox
        else:
            raise ValueError(regime)
        self.E = E

    def ypow(self, e):
        if e == 0:
            return self.one
        if e > 0:
            return self.Y ** e
        if self.regime == "L":
            return GV.mono(e, 1, self.gh)
        # regime S: y^e = g^e s^e, e < 0
        if not (self.sbox > 0):
            raise ArithmeticError("negative y-power at s -> 0 (regime S)")
        return GV.mono(e, self.sbox ** e, self.gh)


def _form_small(F, L, gd, cache):
    """Q = P/(y^m z^p) as a GV (small-z form)."""
    key = ("s", F.name)
    if key in cache:
        return cache[key]
    acc = None
    for j in F.idx:
        h = acb(F.tm[j].eval(L[0], L[1]))
        if j >= F.p:
            t = gd.ypow(F.e[j]) * h
            if j > F.p:
                if gd.Z is None:
                    raise ArithmeticError("small form needs bounded z")
                t = t * gd.Z ** (j - F.p)
        else:
            t = gd.ypow(F.g[j]) * h * gd.W ** (F.p - j)
        acc = t if acc is None else acc + t
    cache[key] = acc
    return acc


def _form_large(F, L, gd, cache):
    """Q / z^{d-p} = P/(y^m z^d) as a GV (large-z form, zeta = 1/z)."""
    key = ("l", F.name)
    if key in cache:
        return cache[key]
    if gd.ZETA is None:
        raise ArithmeticError("large form needs zeta")
    acc = None
    for j in F.idx:
        h = acb(F.tm[j].eval(L[0], L[1]))
        if j >= F.p:
            t = gd.ypow(F.e[j]) * h
            if j < F.d:
                t = t * gd.ZETA ** (F.d - j)
        else:
            t = gd.ypow(F.g[j]) * h * gd.W ** (F.p - j)
            if F.d > F.p:
                t = t * gd.ZETA ** (F.d - F.p)
        acc = t if acc is None else acc + t
    cache[key] = acc
    return acc


def ratio_gv(N, D, k, L, gd, use_large, cache):
    """lam^k N/D as a GV:  N/D = y^{mN-mD} z^{pN-pD} QN/QD (small form) or
    y^{mN-mD} z^{dN-dD} QN^L/QD^L (large form); a negative z-power is
    rewritten exactly as y^{-1} w per power (z w = y)."""
    if use_large:
        q = _form_large(N, L, gd, cache) / _form_large(D, L, gd, cache)
        zp = N.d - D.d
        if zp > 0:
            raise ArithmeticError("positive z-power in a large-form ratio %s/%s" % (N.name, D.name))
        if zp < 0:
            q = q * gd.ZETA ** (-zp)
        ye = N.m - D.m
    else:
        q = _form_small(N, L, gd, cache) / _form_small(D, L, gd, cache)
        zp = N.p - D.p
        ye = N.m - D.m
        if zp < 0:
            q = q * gd.W ** (-zp)
            ye += zp
        elif zp > 0:
            if gd.Z is None:
                raise ArithmeticError("z-power needs bounded z")
            q = q * gd.Z ** zp
    return q * gd.ypow(ye) * (gd.LAM ** k if k else gd.one)


def _real_gv(lc, k, L, gd):
    """lam^k * LC (value = y^off tm(y)) as a GV."""
    h = acb(lc.tm.eval(L[0], L[1]))
    return gd.ypow(lc.off) * h * (gd.LAM ** k if k else gd.one)


def regime_of(sb):
    sl, sh = sb
    if sl >= S_SPLIT:
        return "L"
    if sh <= S_SPLIT or sl <= 0:
        if sh == INF:
            raise ArithmeticError("s-range [0, inf]: split s")
        return "S"
    return "L"


def sigma_graded(ld, box, geo=None, regime=None, terms=False):
    """The graded Sigma (GV) of the box and its Grade; Sigma = S lam^2 it^2
    (S/s^2 in regime L, S in regime S)."""
    g = tl.box_geometry(ld, box) if geo is None else geo
    if g is None:
        return None, None
    L, E = g["L"], g["E"]
    sb, tb = g["s"], g["tau"]
    rg = regime or regime_of(sb)
    gd = Grade(rg, L, E, sb, tb)
    use_large = rg == "L"
    F = ld.forms
    cache = {}
    cent = getattr(ld, "centring", {})

    def R(n, d, k):
        if (n, d) in cent:
            c = cent[(n, d)]
            r = ratio_gv(F["C_%s_%s" % (n, d)], F[d], k, L, gd, use_large, cache)
            if c:                         # N/D = c + (N - c D)/D, exactly
                r = r + (gd.LAM ** k if k else gd.one) * arb(c)
            return r
        return ratio_gv(F[n], F[d], k, L, gd, use_large, cache)
    c = tl.cval()
    t = _real_gv(ld.tt, 1, L, gd)
    tp = _real_gv(ld.ttp, 2, L, gd)
    u2h = _real_gv(ld.u2, 1, L, gd)
    lam = gd.LAM
    it = gd.IT
    Eg = GV.const(E, gd.gh)
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
    cf = c * it * it + (arb(1) / c)
    T = dict(reB1=u2h * lam * (Eg * Eg * it * it - 1) + Eg * it * it * G.real
             - it * G.imag + it * it * H.real,
             imB1=-(2 * Eg * u2h * lam + G.real + Eg * it * G.imag + it * H.imag) * (arb(1) / c),
             fp=-cf * tp * arb(0.5) + (it * (arb(1) / c)) * nup * arb(0.5),
             rhof=-rho * (cf * t - (it * (arb(1) / c)) * nu) * arb(0.5),
             nu2=nu * nu * it * it * arb(0.25))
    if terms:
        return T, gd
    tot = T["reB1"] + T["imB1"] + T["fp"] + T["rhof"] + T["nu2"]
    return tot, gd


def sigma_box_centre(ld, box, geo=None):
    """Ball enclosing Sigma_c = S/(1 + (y tau)^2) over the box (or None if
    empty).  Raises ArithmeticError if a denominator hull contains 0 or the
    graded cancellation is not exact.  Plugs into tc_layer.certify_piece /
    run_layer_piece as sigma_fn."""
    tot, gd = sigma_graded(ld, box, geo)
    if tot is None:
        return None
    v = tot.value()
    if not v.imag.contains(0):
        raise ArithmeticError("Sigma not real")
    return v.real * gd.norm


# =====================================================================
# point evaluation (the float design check; not a certificate)
# =====================================================================
def point_pairs_centre(at, ce, y, prec=320):
    """build_pairs at the single point y (constant Taylor models), from the
    certified centre series; the component LCs carry their Laurent offsets."""
    from . import tc_commutator as tc
    from .tc_tm import LC, TM, AtomKP, centre_components
    with precision(prec):
        comps, _ = centre_components(ce, y)
        vals = [cp.tm.eval(y, y) * arb(y) ** cp.off for cp in comps]
        cc = tuple(LC(TM.const(v, y), 0) for v in vals)
        ak = AtomKP(at, cc, sd_shift=0)
        return tc.build_pairs(ak)


def point_sigma_direct(pr, y, eps, tau, prec=320):
    """Sigma_c = S/(1 + (y tau)^2) at a point from the direct chain
    (tc_layer.point_S_direct; independent of the pivot-form algebra)."""
    S, _ = tl.point_S_direct(pr, eps, tau, prec)
    with precision(prec):
        s = arb(y) * arb(tau)
        return S / (1 + s * s), S


def point_box(y, eps, tau):
    """A degenerate box at (y, eps, tau) for sigma_box_centre."""
    s = y * tau
    return (y, y, eps, eps, s * (1 - 1e-15), s * (1 + 1e-15), tau, tau)
