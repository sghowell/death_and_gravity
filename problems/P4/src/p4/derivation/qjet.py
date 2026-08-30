"""Exact truncated Taylor jets over Q in n variables (forward-mode automatic differentiation).

A ``QJet`` is a polynomial  sum_e c_e X^e  in n formal displacement variables, truncated to
total degree <= ``order``, with ``fractions.Fraction`` coefficients.  Products truncate to the
smaller order of the factors; ``d(i)`` differentiates (order drops by one); ``val`` is the
constant term.  Division uses the geometric series of the reciprocal, exact to the order.

This is deliberately *not* sympy: a jet evaluated at a rational point is a finite exact
object, so identities checked with it hold exactly at that point (used for pointwise
identity testing in ``independent_check``).
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product


def _F(v):
    return v if isinstance(v, Fraction) else Fraction(v)


class QJet:
    __slots__ = ("c", "n", "order")

    def __init__(self, n, order, coeffs=None):
        self.n, self.order = n, order
        self.c = {}
        for e, v in (coeffs or {}).items():
            if sum(e) <= order and v != 0:
                self.c[tuple(e)] = _F(v)

    # -- constructors ---------------------------------------------------------------
    @classmethod
    def const(cls, n, order, v):
        return cls(n, order, {(0,) * n: v})

    @classmethod
    def var(cls, n, order, i, v0):
        """The coordinate function x_i with value v0 (jet: v0 + X_i)."""
        e = [0] * n
        e[i] = 1
        return cls(n, order, {(0,) * n: v0, tuple(e): 1})

    # -- accessors ---------------------------------------------------------------------
    @property
    def val(self):
        return self.c.get((0,) * self.n, Fraction(0))

    def coeff(self, e):
        return self.c.get(tuple(e), Fraction(0))

    def d(self, i):
        """Partial derivative w.r.t. displacement variable i (order - 1)."""
        out = {}
        for e, v in self.c.items():
            if e[i] > 0:
                f = list(e)
                f[i] -= 1
                out[tuple(f)] = v * e[i]
        return QJet(self.n, self.order - 1, out)

    # -- arithmetic --------------------------------------------------------------------
    def _lift(self, other):
        if isinstance(other, QJet):
            return other
        return QJet.const(self.n, self.order, other)

    def __add__(self, other):
        o = self._lift(other)
        k = min(self.order, o.order)
        out = dict(self.c)
        for e, v in o.c.items():
            out[e] = out.get(e, Fraction(0)) + v
        return QJet(self.n, k, {e: v for e, v in out.items() if sum(e) <= k})

    __radd__ = __add__

    def __neg__(self):
        return QJet(self.n, self.order, {e: -v for e, v in self.c.items()})

    def __sub__(self, other):
        return self + (-self._lift(other))

    def __rsub__(self, other):
        return self._lift(other) - self

    def __mul__(self, other):
        o = self._lift(other)
        k = min(self.order, o.order)
        out = {}
        for e1, v1 in self.c.items():
            for e2, v2 in o.c.items():
                e = tuple(i + j for i, j in zip(e1, e2))
                if sum(e) <= k:
                    out[e] = out.get(e, Fraction(0)) + v1 * v2
        return QJet(self.n, k, out)

    __rmul__ = __mul__

    def __pow__(self, m):
        assert isinstance(m, int) and m >= 0
        out = QJet.const(self.n, self.order, 1)
        for _ in range(m):
            out = out * self
        return out

    def inv(self):
        """1/self via 1/(c + d) = (1/c) sum_k (-d/c)^k, exact to the order."""
        c0 = self.val
        if c0 == 0:
            raise ZeroDivisionError("jet with zero constant term")
        d = (self - c0) * Fraction(1, 1)
        u = d * Fraction(-1) * (Fraction(1) / c0)
        out = QJet.const(self.n, self.order, 1)
        term = QJet.const(self.n, self.order, 1)
        for _ in range(self.order):
            term = term * u
            out = out + term
        return out * (Fraction(1) / c0)

    def __truediv__(self, other):
        o = self._lift(other)
        return self * o.inv()

    def __rtruediv__(self, other):
        return self._lift(other) * self.inv()

    def __eq__(self, other):
        o = self._lift(other)
        return self.c == o.c

    def __repr__(self):
        return f"QJet({self.n},{self.order},{self.c})"


def exponents(n, order):
    """All exponent tuples with total degree <= order."""
    return [e for e in product(range(order + 1), repeat=n) if sum(e) <= order]
