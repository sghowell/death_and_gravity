"""Bivariate truncated series  sum_{n<cap} sum_{k<dm} c[n,k] t^n delta^k  over arb.

Used to carry the dependence on a scalar parameter (delta = V0 - c at the sonic
point, delta = mu - c at the centre) *through* the coefficient recursion as a
truncated Taylor polynomial instead of an interval, which removes the
exponential wrapping blow-up of naive interval propagation.  The recursion is
triangular in the delta-degree, so truncating at delta^{dm} is exact for every
retained coefficient.

Encoding (Kronecker): one ``arb_poly`` in z with c[n,k] at position n*B + k,
block size B = 2*dm - 1, so that a product of two encoded series puts
c[n1,k1] c[n2,k2] at (n1+n2)*B + (k1+k2) with k1+k2 <= 2dm-2 < B (no block
overlap); after each product the entries k >= dm are dropped.  Multiplication
is therefore a single C-level polynomial product.

``s[n]`` returns the delta-coefficient of t^n as a ``Series`` (cap = dm), so
``BiSeries`` can be used wherever ``Series`` is (same duck-typed interface).
"""
from __future__ import annotations

from flint import arb, arb_poly

from .arbseries import Series, to_arb


class BiSeries:
    __slots__ = ("cap", "dm", "p")

    def __init__(self, p, cap, dm, normalize=True):
        self.cap, self.dm = cap, dm
        self.p = self._normalize(p) if normalize else p

    @property
    def B(self):
        return 2 * self.dm - 1

    # -- constructors -------------------------------------------------------
    @classmethod
    def from_blocks(cls, blocks, cap, dm):
        """``blocks``: list over n of Series-in-delta (or arb / int) coefficients."""
        B = 2 * dm - 1
        c = [arb(0)] * (len(blocks) * B)
        for n, b in enumerate(blocks):
            if isinstance(b, Series):
                bc = b.coeffs(dm)
            else:
                bc = [to_arb(b)] + [arb(0)] * (dm - 1)
            for k in range(dm):
                c[n * B + k] = bc[k]
        return cls(arb_poly(c), cap, dm)

    @classmethod
    def const(cls, c, cap, dm):
        return cls.from_blocks([c], cap, dm)

    @classmethod
    def var(cls, cap, dm):
        return cls.from_blocks([arb(0), arb(1)], cap, dm)

    def zero(self):
        return BiSeries(arb_poly([arb(0)]), self.cap, self.dm, normalize=False)

    def scalar(self, c):
        return BiSeries.const(c, self.cap, self.dm)

    # -- internals ----------------------------------------------------------
    def _normalize(self, p):
        """Drop delta-degrees >= dm in every block and t-degrees >= cap."""
        B, dm = self.B, self.dm
        c = p.coeffs()
        nmax = len(c) if self.cap is None else min(len(c), self.cap * B)
        out = []
        for pos in range(nmax):
            k = pos % B
            out.append(c[pos] if k < dm else arb(0))
        return arb_poly(out) if out else arb_poly([arb(0)])

    def _cap(self, o):
        if self.cap is None:
            return o.cap
        if o.cap is None:
            return self.cap
        return min(self.cap, o.cap)

    def _lift(self, o):
        if isinstance(o, BiSeries):
            return o
        return BiSeries.const(o, self.cap, self.dm)

    # -- access ------------------------------------------------------------
    def __len__(self):
        return (self.p.length() + self.B - 1) // self.B

    def __getitem__(self, n):
        """delta-Series coefficient of t^n (cap = dm)."""
        B = self.B
        return Series([self.p[n * B + k] for k in range(self.dm)], self.dm)

    def coeffs(self, n=None):
        m = len(self) if n is None else n
        return [self[i] for i in range(m)]

    def with_cap(self, cap):
        if cap is None or (self.cap is not None and cap >= self.cap):
            return BiSeries(self.p, cap, self.dm, normalize=False)
        return BiSeries(self.p, cap, self.dm)

    # -- arithmetic ---------------------------------------------------------
    def __add__(self, o):
        o = self._lift(o)
        return BiSeries(self.p + o.p, self._cap(o), self.dm)

    __radd__ = __add__

    def __neg__(self):
        return BiSeries(-self.p, self.cap, self.dm, normalize=False)

    def __sub__(self, o):
        o = self._lift(o)
        return BiSeries(self.p - o.p, self._cap(o), self.dm)

    def __rsub__(self, o):
        o = self._lift(o)
        return BiSeries(o.p - self.p, self._cap(o), self.dm)

    def __mul__(self, o):
        if isinstance(o, BiSeries):
            return BiSeries(self.p * o.p, self._cap(o), self.dm)
        if isinstance(o, Series):                      # delta-polynomial scalar
            return self * BiSeries.const(o, self.cap, self.dm)
        return BiSeries(self.p * arb_poly([to_arb(o)]), self.cap, self.dm, normalize=False)

    __rmul__ = __mul__

    def __pow__(self, n):
        n = int(n)
        assert n >= 0
        out = BiSeries.const(1, self.cap, self.dm)
        base = self
        while n:
            if n & 1:
                out = out * base
            n >>= 1
            if n:
                base = base * base
        return out

    def _blockwise(self, factor, shift):
        """New series with block n-shift = factor(n) * block n."""
        B, dm = self.B, self.dm
        c = self.p.coeffs()
        nblocks = (len(c) + B - 1) // B
        out = [arb(0)] * (max(nblocks - shift, 1) * B)
        for n in range(shift, nblocks):
            f = factor(n)
            for k in range(dm):
                pos = n * B + k
                if pos < len(c):
                    out[(n - shift) * B + k] = c[pos] * f
        return BiSeries(arb_poly(out), self.cap, self.dm, normalize=False)

    def deriv(self):
        return self._blockwise(lambda n: n, 1)

    def euler(self):
        return self._blockwise(lambda n: n, 0)

    def __call__(self, x):
        """Horner evaluation in t at the ball ``x``; returns a delta-Series."""
        acc = Series([arb(0)], self.dm)
        for n in reversed(range(len(self))):
            acc = acc * to_arb(x) + self[n]
        return acc

    def __repr__(self):
        return f"BiSeries(len={len(self)}, cap={self.cap}, dm={self.dm})"
