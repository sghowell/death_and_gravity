"""Truncated power series over ``flint.arb`` balls.

A ``Series`` wraps an ``arb_poly`` together with an optional truncation length
``cap`` (``None`` = keep every coefficient, used for the untruncated residual of
a polynomial approximation).  All arithmetic is ball arithmetic, so every
coefficient of a computed series is a rigorous enclosure of the exact
coefficient of the same expression evaluated on the exact coefficients.

Conventions: ``s[i]`` is the coefficient of t^i (0 for i >= length);
``s.deriv()`` is d/dt, ``s.euler()`` is t d/dt.  ``l1nu(s, nu)`` is the weighted
norm  sum_n |s_n| nu^n  (an upper-bound ball).
"""
from __future__ import annotations

from contextlib import contextmanager

from flint import arb, arb_poly, ctx, fmpq

DEFAULT_PREC = 256


@contextmanager
def precision(prec=DEFAULT_PREC):
    """Temporarily set the Arb working precision (bits)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        yield
    finally:
        ctx.prec = old


def to_arb(x):
    """Exact-as-possible conversion of ints/fmpq/str/arb to ``arb`` (floats are exact binary)."""
    if isinstance(x, arb):
        return x
    if isinstance(x, (int, fmpq, str)):
        return arb(x)
    if isinstance(x, float):
        return arb(x)          # exact binary value; callers wanting a decimal should pass str
    raise TypeError(f"cannot convert {type(x)} to arb")


def ball(mid, rad=0):
    """arb with the given midpoint (decimal string allowed) and radius."""
    m = to_arb(mid)
    return m if rad == 0 else m + arb(0, rad)


class Series:
    __slots__ = ("cap", "p")

    def __init__(self, coeffs, cap=None):
        if isinstance(coeffs, arb_poly):
            p = coeffs
        else:
            p = arb_poly([to_arb(c) for c in coeffs])
        if cap is not None and p.length() > cap:
            p = arb_poly(p.coeffs()[:cap])
        self.p, self.cap = p, cap

    # -- constructors -------------------------------------------------------
    @classmethod
    def const(cls, c, cap=None):
        return cls([to_arb(c)], cap)

    @classmethod
    def var(cls, cap=None):
        return cls([arb(0), arb(1)], cap)

    # -- basic access -------------------------------------------------------
    def __len__(self):
        return self.p.length()

    def __getitem__(self, i):
        if i < 0 or i >= self.p.length():
            return arb(0)
        return self.p[i]

    def coeffs(self, n=None):
        """List of the first ``n`` coefficients (zero-padded)."""
        c = self.p.coeffs()
        if n is None:
            return c
        return (c + [arb(0)] * n)[:n]

    def with_cap(self, cap):
        return Series(self.p, cap)

    # -- arithmetic ---------------------------------------------------------
    def _cap(self, o):
        if self.cap is None:
            return o.cap
        if o.cap is None:
            return self.cap
        return min(self.cap, o.cap)

    def __add__(self, o):
        if isinstance(o, Series):
            return Series(self.p + o.p, self._cap(o))
        return Series(self.p + arb_poly([to_arb(o)]), self.cap)

    __radd__ = __add__

    def __neg__(self):
        return Series(-self.p, self.cap)

    def __sub__(self, o):
        if isinstance(o, Series):
            return Series(self.p - o.p, self._cap(o))
        return Series(self.p - arb_poly([to_arb(o)]), self.cap)

    def __rsub__(self, o):
        return Series(arb_poly([to_arb(o)]) - self.p, self.cap)

    def __mul__(self, o):
        if isinstance(o, Series):
            return Series(self.p * o.p, self._cap(o))
        if not isinstance(o, (int, float, fmpq, arb)):
            return NotImplemented                     # e.g. BiSeries: use its __rmul__
        return Series(self.p * arb_poly([to_arb(o)]), self.cap)

    __rmul__ = __mul__

    def __pow__(self, n):
        n = int(n)
        assert n >= 0
        out = Series.const(1, self.cap)
        base = self
        while n:
            if n & 1:
                out = out * base
            n >>= 1
            if n:
                base = base * base
        return out

    def deriv(self):
        """d/dt: coefficient i is (i+1) s_{i+1}."""
        c = self.p.coeffs()
        return Series([c[i] * i for i in range(1, len(c))] or [arb(0)],
                      None if self.cap is None else self.cap)

    def euler(self):
        """t d/dt: coefficient i is i s_i."""
        c = self.p.coeffs()
        return Series([c[i] * i for i in range(len(c))] or [arb(0)], self.cap)

    def __call__(self, x):
        """Evaluate the (finite) polynomial at the ball ``x`` by Horner's rule."""
        c = self.p.coeffs()
        acc = arb(0)
        for ci in reversed(c):
            acc = acc * x + ci
        return acc

    # -- ring helpers used generically by the evaluator ----------------------
    def zero(self):
        return Series([arb(0)], self.cap)

    def scalar(self, c):
        return Series([to_arb(c)], self.cap)

    def inv(self):
        """1/s as a truncated series (cap required); s_0 must be provably nonzero."""
        assert self.cap is not None
        c = self.coeffs(self.cap)
        if not c[0] != 0:
            raise ZeroDivisionError("series inverse: constant term not provably nonzero")
        q = [1 / c[0]]
        for n in range(1, self.cap):
            s = arb(0)
            for j in range(1, n + 1):
                s += c[j] * q[n - j]
            q.append(-s / c[0])
        return Series(q, self.cap)

    def __truediv__(self, o):
        if isinstance(o, Series):
            return self * o.inv()
        return Series(self.p * arb_poly([1 / to_arb(o)]), self.cap)

    def sqrt(self):
        """Truncated series square root; s_0 must be provably positive."""
        assert self.cap is not None
        c = self.coeffs(self.cap)
        if not c[0] > 0:
            raise ValueError("series sqrt: constant term not provably positive")
        r = [c[0].sqrt()]
        for n in range(1, self.cap):
            s = arb(0)
            for j in range(1, n):
                s += r[j] * r[n - j]
            r.append((c[n] - s) / (2 * r[0]))
        return Series(r, self.cap)

    def __repr__(self):
        return f"Series(len={len(self)}, cap={self.cap})"


# ----------------------------------------------------------------------------
# norms
# ----------------------------------------------------------------------------
def abs_upper(x):
    """A (nonnegative, exact-endpoint) arb bounding |x| from above."""
    return x.abs_upper()


def l1nu(s, nu, start=0, stop=None):
    """Upper bound (as an arb) of  sum_{n=start}^{stop-1} |s_n| nu^n."""
    nu = to_arb(nu)
    c = s.coeffs()
    if stop is None:
        stop = len(c)
    tot = arb(0)
    w = nu ** start
    for n in range(start, min(stop, len(c))):
        tot += abs_upper(c[n]) * w
        w *= nu
    return tot


def max_norm(vec):
    """|v|_inf as an upper-bound arb for a list of arb."""
    m = arb(0)
    for v in vec:
        m = m.max(abs_upper(v))
    return m


def contains_float(b, x, tol=0.0):
    """True if the ball ``b`` (enlarged by ``tol``) contains the float ``x``."""
    return (b + arb(0, tol)).contains(arb(x))
