"""Taylor-model (TM) core for the S4-3 certificates (Lemma D closure + commutator).

A TM is a polynomial in the local variable ell >= 0 with arb-ball coefficients,
valid on ell in [0, R]: it encloses a function f pointwise,
    f(ell) in sum_k c_k ell^k  for every ell in [0, R]
(coefficients are balls; truncation spill is absorbed into the coefficient at
CAP-1 using ell^k <= ell^{CAP-1} R^{k-CAP+1}, which is pointwise valid since
ell >= 0).  TMs are never differentiated: x-derivatives are carried as separate
TMs built at the component level (exact series ops + Cauchy tails), and all
composite derivatives are formed by jet (Leibniz) arithmetic.

Layers:
  TM              ball-coefficient polynomial on [0, R]
  LC = (TM, off)  Laurent coefficient: value = ell^off * TM(ell), off in Z
  kp              kappa-polynomial: python list of LC (ascending kappa powers)
  jet             list [f, f', f'', ...] of kp (x-derivative slots)
  pair            (num_jets, den_jets): rational function num/den, never divided

The three regions build component TMs (and their x-derivative jets) from the
certified data: centre series (ell = t = e^x), sonic series (ell = xi = -x),
tube steps (ell = sigma = x_k - x).  Atom kappa-polys are assembled from the
compiled fmpq term lists of ``tc_energy.Atoms``.
"""
from __future__ import annotations

from flint import arb

CAP = 80                      # TM truncation order (coefficients 0..CAP-1)
SPLIT_REL = [1e-12]           # near-zero threshold for LC.split_junk (mutable knob)


class TM:
    """ve: the enclosure is f = P + e with P the ball-coefficient polynomial
    below order ve and val(e) >= ve (exact-coefficient semantics up to ve);
    drop_lead(k) is sound only for k <= ve (plus an exact identity)."""

    __slots__ = ("c", "R", "ve")

    def __init__(self, c, R, ve=CAP - 1):
        self.c, self.R, self.ve = list(c), R, min(ve, CAP - 1)
        if len(self.c) > CAP:
            spill = arb(0)
            Rb = arb(R)
            for k in range(CAP, len(self.c)):
                spill += abs(self.c[k]) * Rb ** (k - CAP + 1)
            self.c = self.c[:CAP]
            while len(self.c) < CAP:
                self.c.append(arb(0))
            self.c[CAP - 1] += arb(0, spill.upper())
        if not all(x.is_finite() for x in self.c):
            raise ArithmeticError("non-finite TM coefficient")

    @classmethod
    def const(cls, v, R):
        return cls([v if isinstance(v, arb) else arb(v)], R)

    def _ve2(self, o):
        return min(self.ve, o.ve if isinstance(o, TM) else CAP - 1)

    def __add__(self, o):
        if not isinstance(o, TM):
            o = TM.const(o, self.R)
        n = max(len(self.c), len(o.c))
        c = [arb(0)] * n
        for i, x in enumerate(self.c):
            c[i] += x
        for i, x in enumerate(o.c):
            c[i] += x
        return TM(c, self.R, self._ve2(o))

    def __neg__(self):
        return TM([-x for x in self.c], self.R, self.ve)

    def __sub__(self, o):
        return self + (-o if isinstance(o, TM) else TM.const(-arb(o), self.R))

    def __mul__(self, o):
        if not isinstance(o, TM):
            return TM([x * o for x in self.c], self.R, self.ve)
        c = [arb(0)] * (len(self.c) + len(o.c) - 1)
        for i, x in enumerate(self.c):
            if x.is_zero():
                continue
            for j, y in enumerate(o.c):
                c[i + j] += x * y
        return TM(c, self.R, self._ve2(o))

    __rmul__ = __mul__

    def shift(self, k):
        """ell^k * self (k >= 0): exact coefficient shift."""
        return TM([arb(0)] * k + self.c, self.R, min(self.ve + k, CAP - 1))

    def drop_lead(self, k, tol=None):
        """Divide by ell^k: requires the k lowest coefficients to enclose 0
        (their exact vanishing must be justified by the caller -- an exact
        identity).  Residual widths are absorbed pointwise: |c_j ell^j| <=
        |c_j| ell^k R^{j-k}... only valid j<k contributions are the ball widths,
        added to the new constant term as width/R^{k-j} <= width * (1/R)^{k-j}
        -- we bound |c_j ell^{j-k}| on [0,R] only for j>=k; for j<k the term
        c_j ell^{j-k} is unbounded at 0, so the enclosure-of-zero widths must be
        certified exactly zero by the caller; we assert containment of 0 and
        drop them."""
        if k > self.ve:
            raise ArithmeticError("drop_lead beyond the exact-coefficient range ve")
        for j in range(k):
            if j < len(self.c) and not self.c[j].contains(arb(0)):
                raise ArithmeticError("drop_lead: low coefficient does not contain 0")
        return TM(self.c[k:] if len(self.c) > k else [arb(0)], self.R, self.ve - k)

    def eval(self, lo, hi, J=8):
        """Ball enclosing self(ell) for ell in [lo, hi] subset [0, R]:
        recentred Taylor form via J exact synthetic divisions by (x - mid),
        f(x) = sum_{j<=J} g_j (x-mid)^j + (x-mid)^{J+1} q(x),
        so the g_j come out point-tight (ball cancellations happen exactly)
        and the interval-evaluated q is suppressed by delta^{J+1}."""
        if hi <= lo or len(self.c) <= 2:
            xb = arb((lo + hi) / 2, (hi - lo) / 2)
            acc = arb(0)
            for x in reversed(self.c):
                acc = acc * xb + x
            return acc
        mid = arb((lo + hi) / 2)
        delta = arb(0, (hi - lo) / 2)
        c = list(self.c)
        g = []
        for _ in range(min(J + 1, len(c))):
            # synthetic division by (x - mid): b_k = c_k + mid b_{k+1}
            b_ = [arb(0)] * len(c)
            b_[-1] = c[-1]
            for k in range(len(c) - 2, -1, -1):
                b_[k] = c[k] + mid * b_[k + 1]
            g.append(b_[0])
            c = b_[1:]
            if not c:
                break
        acc = arb(0)
        if c:
            xb = arb((lo + hi) / 2, (hi - lo) / 2)
            for x in reversed(c):
                acc = acc * xb + x
        out = acc
        for gj in reversed(g):
            out = out * delta + gj
        return out

    def abs_up(self, lo, hi):
        return abs(self.eval(lo, hi)).upper()

    def inv(self):
        """1/self on [0, R]; needs |self - c0| / |c0| hull < 1."""
        c0 = self.c[0]
        if not (abs(c0).lower() > 0):
            raise ArithmeticError("TM inverse: constant term contains 0")
        q = TM([arb(0)] + [-x / c0 for x in self.c[1:]], self.R)   # self = c0(1-q)
        qh = abs(q.eval(0.0, self.R))
        if not (qh.upper() < 1):
            raise ArithmeticError("TM inverse: |q| hull >= 1 (R=%g len=%d "
                                  "c0=%s qh=%s)" % (self.R, len(self.c),
                                                    str(self.c[0])[:20],
                                                    str(qh)[:20]))
        out, p = TM.const(arb(1), self.R), TM.const(arb(1), self.R)
        n = max(2, len(self.c))
        for _ in range(n - 1):
            p = p * q
            out = out + p
        # geometric tail |q|^n/(1-|q|): q has exact valuation vq, so the tail
        # is a function of valuation >= n*vq; place the ball there (pointwise
        # valid: |tail(ell)| <= tailbound and, refined, <= tailbound*(ell/R)^m).
        vq = 0
        for x in q.c:
            if x.is_zero():
                vq += 1
            else:
                break
        vq = min(vq, len(q.c) - 1)
        m = min(n * vq, CAP - 1)
        qhat = abs(TM(q.c[vq:] or [arb(0)], self.R).eval(0.0, self.R)).upper()
        tail = (qhat ** n * arb(self.R) ** (n * vq - m) / (1 - qh).lower()).upper()
        out = out + TM([arb(0)] * m + [arb(0, tail)], self.R, m)
        return out * (1 / c0)


def tm_zero(R):
    return TM([arb(0)], R)


# ------------------------------------------------------------------ Laurent layer
class LC:
    """value = ell^off * tm(ell)."""
    __slots__ = ("tm", "off")

    def __init__(self, tm, off=0):
        self.tm, self.off = tm, off

    def __mul__(self, o):
        if isinstance(o, LC):
            return LC(self.tm * o.tm, self.off + o.off)
        return LC(self.tm * o, self.off)

    __rmul__ = __mul__

    def __neg__(self):
        return LC(-self.tm, self.off)

    def __add__(self, o):
        if not isinstance(o, LC):
            o = LC(TM.const(o, self.tm.R), 0)
        if self.off == o.off:
            return LC(self.tm + o.tm, self.off)
        if self.off < o.off:
            return LC(self.tm + o.tm.shift(o.off - self.off), self.off)
        return LC(self.tm.shift(self.off - o.off) + o.tm, o.off)

    def __sub__(self, o):
        return self + (-o)

    def trim(self):
        """Strip leading exactly-zero TM coefficients into the offset (sound:
        dropping exact zeros needs no identity; capped by ve)."""
        k = 0
        for x in self.tm.c:
            if x.is_zero():
                k += 1
            else:
                break
        k = min(k, len(self.tm.c) - 1)
        if k <= 0:
            return self
        # dropping exactly-zero coefficients is pointwise-free (they contribute
        # nothing to any hull), no ve constraint applies
        tm = TM(self.tm.c[k:], self.tm.R, max(self.tm.ve - k, 0))
        return LC(tm, self.off + k)

    def split_junk(self, rel=None):
        """Exact decomposition  self = junk + ell^k * content  (as LCs), where
        the k leading TM coefficients are near-zero balls (contain 0 and have
        radius <= rel * coefficient scale).  Nothing is dropped: both parts are
        returned, so no identity is needed.  Returns (content_LC, junk_LC)."""
        if rel is None:
            rel = SPLIT_REL[0]
        t = self.trim()
        S = max((abs(x).upper() for x in t.tm.c), default=arb(0))
        k = 0
        for x in t.tm.c[:-1]:
            if x.contains(arb(0)) and float(x.rad()) <= rel * float(S):
                k += 1
            else:
                break
        if k == 0:
            return t, None
        junk = LC(TM(t.tm.c[:k], t.tm.R, ve=0), t.off)
        cont = LC(TM(t.tm.c[k:], t.tm.R, ve=max(t.tm.ve - k, 0)), t.off + k)
        return cont, junk

    def inv(self):
        t = self.trim()
        return LC(t.tm.inv(), -t.off)

    def normalized(self, off):
        """Re-express at the given (lower or equal) offset."""
        if off > self.off:
            raise ValueError("cannot raise offset without drop_lead")
        return LC(self.tm.shift(self.off - off), off) if off < self.off else self

    def eval(self, lo, hi):
        """Ball for ell in [lo, hi]; needs lo > 0 if off < 0."""
        if self.off < 0 and lo <= 0:
            raise ArithmeticError("negative Laurent offset evaluated at 0")
        xb = arb((lo + hi) / 2, (hi - lo) / 2)
        return self.tm.eval(lo, hi) * xb ** self.off

    def abs_up(self, lo, hi):
        return abs(self.eval(lo, hi)).upper()


def lc_zero(R):
    return LC(tm_zero(R), 0)


def lc_raise(x, k):
    """Re-express the LC at offset off+k: sound iff the k lowest TM coefficients
    are exactly zero (caller supplies the identity; drop_lead checks ve and 0-containment)."""
    return LC(x.tm.drop_lead(k), x.off + k)


# ------------------------------------------------------------- jet (JLC) layer
# A JLC is a list [f, f', f'', ...] of LC: the value and x-derivatives of one
# scalar quantity.  Never differentiated; combined by Leibniz.
_BIN = ((1,), (1, 1), (1, 2, 1), (1, 3, 3, 1))


def jlc_mul(x, y, depth=None):
    n = min(len(x), len(y)) if depth is None else depth + 1
    return [sum((_BIN[k][i] * (x[i] * y[k - i]) for i in range(1, k + 1)),
                _BIN[k][0] * (x[0] * y[k])) for k in range(n)]


def jlc_add(*xs):
    n = min(len(x) for x in xs)
    return [sum((x[k] for x in xs[1:]), xs[0][k]) for k in range(n)]


def jlc_neg(x):
    return [-v for v in x]


def jlc_scale(x, c):
    return [v * c for v in x]


def jlc_inv(x):
    """Jet of 1/f to depth len(x)-1 <= 2."""
    g0 = x[0].inv()
    out = [g0]
    if len(x) > 1:
        out.append(-(x[1] * g0 * g0))
    if len(x) > 2:
        out.append(2 * (x[1] * x[1] * g0 * g0 * g0) - x[2] * g0 * g0)
    return out


def jlc_deriv(x):
    return x[1:]


# kappa-polynomials with JLC coefficients (ascending kappa)
def kpj_mul(a, b, depth=None):
    out = [None] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            t = jlc_mul(x, y, depth)
            out[i + j] = t if out[i + j] is None else jlc_add(out[i + j], t)
    return out


def kpj_add(*ps):
    n = max(len(p) for p in ps)
    zR = ps[0][0][0].tm.R
    d = min(len(p[0]) for p in ps)
    out = [[lc_zero(zR) for _ in range(d)] for _ in range(n)]
    for p in ps:
        for i, x in enumerate(p):
            out[i] = jlc_add(out[i], x[:d])
    return out


def kpj_scale(p, c):
    return [jlc_scale(x, c) for x in p]


def kpj_neg(p):
    return [jlc_neg(x) for x in p]


def kpj_deriv(p):
    return [jlc_deriv(x) for x in p]


def kpj_depth(p, d):
    return [x[:d + 1] for x in p]


def kpj_kmul(p, jlc0=None):
    """Multiply by kappa (shift up); jlc0 fills the constant slot."""
    zR = p[0][0].tm.R
    d = len(p[0])
    z = [lc_zero(zR) for _ in range(d)]
    return [jlc0[:d] if jlc0 is not None else z] + list(p)


# ------------------------------------------------------------- atom evaluation
class AtomKP:
    """Assemble kappa-polys of JLC for the compiled atoms of tc_energy.Atoms,
    from component value-LCs (A, N, W, V).  Derivative slots use the exact
    D1-atoms via  X' = dX/sd,  X'' = (ddX*sd - dX*dsd)/sd^3  where sd = SD-value
    and the numerators' extra vanishing at a singular end is handled by the
    caller through sd_shift (Laurent offset raise justified by exact identities:
    Delta~(u0) = 0 via D(u0) = 0, PW(u0) = PV(u0) = 0)."""

    def __init__(self, at, comps, sd_shift=0):
        self.at, self.sd_shift = at, comps and sd_shift
        self.A, self.N, self.W, self.V = comps
        self.R = self.A.tm.R
        self._pw = {}
        self._mon = {}
        self.sd_raw = self.value_kp("SD")[0]
        self.dsd_raw = self.value_kp("dSD")[0]
        sd = lc_raise(self.sd_raw, sd_shift) if sd_shift else self.sd_raw
        self.sd = sd                          # Delta~ ~ xi at the sonic point
        self.sdi = sd.inv()
        self.sdi3 = self.sdi * self.sdi * self.sdi

    def _pow(self, comp, n):
        key = (comp, n)
        if key not in self._pw:
            if n == 0:
                self._pw[key] = LC(TM.const(arb(1), self.R), 0)
            else:
                self._pw[key] = self._pow(comp, n - 1) * getattr(self, comp)
        return self._pw[key]

    def _monom(self, a, b, c, e):
        key = (a, b, c, e)
        if key not in self._mon:
            self._mon[key] = (self._pow("A", a) * self._pow("N", b)) * \
                (self._pow("W", c) * self._pow("V", e))
        return self._mon[key]

    def value_kp(self, name):
        """kappa-poly (list of LC) of the atom's value."""
        terms, dmax = self.at.cmp[name]
        kd = max(t[0] for t in terms)
        out = [None] * (kd + 1)
        for j, coef, a, b, c, e in terms:
            t = self._monom(a, b, c, e) * arb(coef)
            out[j] = t if out[j] is None else out[j] + t
        return [x if x is not None else lc_zero(self.R) for x in out]

    def jet_kp(self, name, depth, num_shift=0):
        """kappa-poly of JLC to the given depth (<= 2).  num_shift: extra exact
        vanishing of the D1-numerators at the singular end (= sd_shift there)."""
        v0 = self.value_kp(name)
        out = [[x] for x in v0]
        if depth >= 1:
            d1 = self.value_kp("d" + name)
            for j in range(len(out)):
                x = d1[j] if j < len(d1) else lc_zero(self.R)
                if num_shift:
                    x = lc_raise(x, num_shift)
                out[j].append(x * self.sdi)
        if depth >= 2:
            d1 = self.value_kp("d" + name)
            d2 = self.value_kp("dd" + name)
            for j in range(len(out)):
                x1 = d1[j] if j < len(d1) else lc_zero(self.R)
                x2 = d2[j] if j < len(d2) else lc_zero(self.R)
                # X'' = (ddX*SD - dX*dSD)/SD^3; the numerator's composite value
                # is SD^3 X'' (X'' analytic), so it vanishes to order 3*shift.
                num = x2 * self.sd_raw - x1 * self.dsd_raw
                if num_shift:
                    num = lc_raise(num, 3 * num_shift)
                out[j].append(num * self.sdi3)
        return out


# ----------------------------------------------------------- region components
def centre_components(ce, ymax, prec_tail=None):
    """Component LCs (A, N, W, V) on y = e^x in [0, ymax] from the certified
    centre expansion (balls + Banach tail at order K+1)."""
    K = ce.K
    tb = ce.cert.tail_bound(arb(ymax).upper()).upper()
    tms = []
    for i in range(3):
        c = [ce.balls[k][i] for k in range(K + 1)]
        c += [arb(0, (tb / arb(ymax) ** (K + 1)).upper())]
        tms.append(TM(c, ymax, ve=K + 1))
    n, w, v = tms
    y2 = TM([arb(0), arb(0), arb(1)], ymax)
    S = TM.const(arb(1), ymax) - (v * v) * y2
    Ahat = (2 * w) * (TM.const(arb(1), ymax) + (v * v) * y2 * (arb(1) / 3)
                      + (n * v) * (arb(4) / 3)) * S.inv()
    A = TM.const(arb(1), ymax) + Ahat.shift(2)
    return (LC(A, 0), LC(n, -1), LC(w, 2), LC(v, 1)), (n, w, v)


def sonic_components(ex, xmax):
    """Component LCs (A, N, W, V) on xi = -x in [0, xmax] from the certified
    sonic expansion; coefficients sign-flipped (u(x) = sum b_n (-xi)^n)."""
    K = ex.K
    rho = arb(ex.cert.nu) * arb(0.75)
    if not (rho > arb(xmax)):
        raise ArithmeticError("sonic Cauchy radius does not cover the region")
    tb = (ex.cert.tail_bound(rho) * (arb(xmax) / rho) ** (K + 1)
          / (1 - arb(xmax) / rho)).upper()
    out = []
    for i in range(4):
        c = [ex.balls[k][i] * (-1 if k % 2 else 1) for k in range(K + 1)]
        c += [arb(0, (arb(tb) / arb(xmax) ** (K + 1)).upper())]
        out.append(LC(TM(c, xmax, ve=K + 1), 0))
    return tuple(out)


def sonic_components_recentred(ex, lo, hi):
    """Component LCs (A, N, W, V) on ell = xi - lo in [0, hi - lo] from the
    certified sonic expansion, Taylor-shifted to x0 = -lo (exact synthetic
    shifts of the ball coefficients; tails via the certified Cauchy radius)."""
    K = ex.K
    rho = arb(ex.cert.nu) * arb(0.75)
    if not (rho > arb(hi)):
        raise ArithmeticError("sonic Cauchy radius does not cover the piece")
    R = hi - lo
    tb = (ex.cert.tail_bound(rho) * (arb(hi) / rho) ** (K + 1)
          / (1 - arb(hi) / rho)).upper()
    out = []
    for i in range(4):
        c = [ex.balls[k][i] for k in range(K + 1)]
        # u(x) = sum b_k x^k, x = -(lo + ell): shift then flip signs
        for _ in range(1):
            pass
        # synthetic Taylor shift to x0 = -lo: g_n = sum_k C(k,n) b_k x0^{k-n}
        x0 = arb(-lo)
        g = []
        cc = list(c)
        for n in range(K + 1):
            b_ = [arb(0)] * len(cc)
            if len(cc) == 1:
                g.append(cc[0])
                break
            b_[-1] = cc[-1]
            for k in range(len(cc) - 2, -1, -1):
                b_[k] = cc[k] + x0 * b_[k + 1]
            g.append(b_[0])
            cc = b_[1:]
        # d/d ell = -d/dx: coefficient sign flip
        gg = [g[n] * ((-1) ** n) for n in range(len(g))]
        gg += [arb(0, (arb(tb) / arb(R) ** (K + 1)).upper())]
        out.append(LC(TM(gg, R, ve=K + 1), 0))
    return tuple(out)


def step_components(sd, hi=None):
    """Component LCs (A, N, W, V) on sigma = x_k - x in [0, h] from a certified
    tube step (point Taylor coefficients + eps_z enclosure)."""
    h = float(sd.h) if hi is None else hi
    tms = []
    for i in range(4):
        c = [sd.co[k][i] * (-1 if k % 2 else 1) for k in range(len(sd.co))]
        tm = TM(c, h, ve=0) + TM([arb(0, sd.eps_z)], h, ve=0)
        tms.append(tm)
    n, w, v, t = tms
    ti = t.inv()
    N, W, V = n * ti, w * (t * t), v * t
    S = TM.const(arb(1), h) - V * V
    A = TM.const(arb(1), h) + (2 * W) * (TM.const(arb(1), h) + V * V * (arb(1) / 3)
                                         + N * V * (arb(4) / 3)) * S.inv()
    return (LC(A, 0), LC(N, 0), LC(W, 0), LC(V, 0))
