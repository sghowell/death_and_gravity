"""Truncated power series (TPS) and dual numbers.

Used to (i) generate high-order Taylor expansions of the CSS background at the
sonic point, and (ii) linearize the field equations exactly (forward-mode AD
via ``Dual``) without any symbolic package at runtime.  Only the operations
``+ - * /`` and small integer powers are needed by ``p4.css``; every routine
there is written generically so that it works on floats, complex numbers,
numpy arrays, ``TPS`` and ``Dual`` objects alike.
"""
from __future__ import annotations

import numpy as np

_SCALARS = (int, float, complex, np.integer, np.floating, np.complexfloating, np.ndarray)


class TPS:
    """Power series sum_{n<K} c[n] x^n, truncated at fixed order K."""

    __slots__ = ("c",)

    def __init__(self, c):
        self.c = np.asarray(c)

    @classmethod
    def const(cls, v, K, dtype=float):
        c = np.zeros(K, dtype=np.result_type(dtype, np.asarray(v).dtype))
        c[0] = v
        return cls(c)

    @classmethod
    def var(cls, K, dtype=float):
        c = np.zeros(K, dtype=dtype)
        c[1] = 1.0
        return cls(c)

    @property
    def K(self):
        return len(self.c)

    def _lift(self, o):
        if isinstance(o, TPS):
            return o
        if isinstance(o, _SCALARS):
            return TPS.const(o, self.K, dtype=self.c.dtype)
        return NotImplemented

    def __add__(self, o):
        o = self._lift(o)
        return NotImplemented if o is NotImplemented else TPS(self.c + o.c)

    __radd__ = __add__

    def __neg__(self):
        return TPS(-self.c)

    def __sub__(self, o):
        o = self._lift(o)
        return NotImplemented if o is NotImplemented else TPS(self.c - o.c)

    def __rsub__(self, o):
        o = self._lift(o)
        return NotImplemented if o is NotImplemented else TPS(o.c - self.c)

    def __mul__(self, o):
        if isinstance(o, _SCALARS):
            return TPS(self.c * o)
        if not isinstance(o, TPS):
            return NotImplemented
        K = self.K
        return TPS(np.convolve(self.c, o.c)[:K])

    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, _SCALARS):
            return TPS(self.c / o)
        if not isinstance(o, TPS):
            return NotImplemented
        return TPS(_series_div(self.c, o.c))

    def __rtruediv__(self, o):
        o = self._lift(o)
        return NotImplemented if o is NotImplemented else TPS(_series_div(o.c, self.c))

    def __pow__(self, n):
        assert isinstance(n, (int, np.integer)) and n >= 0
        out = TPS.const(1.0, self.K, dtype=self.c.dtype)
        for _ in range(int(n)):
            out = out * self
        return out

    def deriv(self):
        c = np.zeros_like(self.c)
        c[:-1] = self.c[1:] * np.arange(1, self.K)
        return TPS(c)

    def __call__(self, x):
        return np.polyval(self.c[::-1], x)

    def __repr__(self):
        return f"TPS({self.c!r})"


def _series_div(p, r):
    """Coefficients of p/r (both length-K arrays)."""
    K = len(p)
    q = np.zeros(K, dtype=np.result_type(p.dtype, r.dtype))
    for n in range(K):
        s = p[n] - np.dot(r[1:n + 1], q[n - 1::-1][:n]) if n else p[0]
        q[n] = s / r[0]
    return q


class Dual:
    """a + eps*b with eps^2 = 0.  ``a`` and ``b`` may be scalars, arrays or TPS."""

    __slots__ = ("a", "b")

    def __init__(self, a, b):
        self.a, self.b = a, b

    @staticmethod
    def _parts(o):
        return (o.a, o.b) if isinstance(o, Dual) else (o, 0)

    def __add__(self, o):
        oa, ob = self._parts(o)
        return Dual(self.a + oa, self.b + ob)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.a, -self.b)

    def __sub__(self, o):
        oa, ob = self._parts(o)
        return Dual(self.a - oa, self.b - ob)

    def __rsub__(self, o):
        oa, ob = self._parts(o)
        return Dual(oa - self.a, ob - self.b)

    def __mul__(self, o):
        oa, ob = self._parts(o)
        return Dual(self.a * oa, self.a * ob + self.b * oa)

    __rmul__ = __mul__

    def __truediv__(self, o):
        oa, ob = self._parts(o)
        inv = 1 / oa
        return Dual(self.a * inv, (self.b - self.a * inv * ob) * inv)

    def __rtruediv__(self, o):
        oa, ob = self._parts(o)
        inv = 1 / self.a
        return Dual(oa * inv, (ob - oa * inv * self.b) * inv)

    def __pow__(self, n):
        assert isinstance(n, (int, np.integer)) and n >= 0
        out = Dual(1.0, 0.0)
        for _ in range(int(n)):
            out = out * self
        return out

    def __repr__(self):
        return f"Dual({self.a!r}, {self.b!r})"
