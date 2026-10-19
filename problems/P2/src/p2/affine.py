"""First-order affine forms over two complex noise symbols (S4, regime L).

An AF represents the set  { c0 + c1 dw + c2 dh + E :  |dw| <= Rw, |dh| <= Rh, |E| <= e }  where
c0, c1, c2 are acb balls (any point of each ball may be the true coefficient) and e is an upper bound.
dw is the deviation of omega from the box centre (a complex disc), dh the deviation of h^2 from the
centre of an h^2 sub-interval (the real interval is contained in the disc).  Products keep the linear
terms exactly and put every second-order term into e; there is no rectangle wrapping (the linear
parts are carried as coefficients, not as boxes), which is what the plain acb recurrence suffered.
"""
from __future__ import annotations

from flint import acb, arb


def up(x) -> arb:
    """Exact upper endpoint of an arb (as an arb)."""
    return arb(arb(x).upper())


class AF:
    __slots__ = ("Rh", "Rw", "c0", "c1", "c2", "e")

    def __init__(self, c0, c1=0, c2=0, e=0, Rw=None, Rh=None):
        self.c0, self.c1, self.c2 = acb(c0), acb(c1), acb(c2)
        self.e = up(e)
        self.Rw, self.Rh = Rw, Rh

    def _new(self, c0, c1, c2, e):
        return AF(c0, c1, c2, e, self.Rw, self.Rh)

    def const(self, x):
        return AF(x, 0, 0, 0, self.Rw, self.Rh)

    def lin(self) -> arb:
        return up(abs(self.c1).upper() * self.Rw + abs(self.c2).upper() * self.Rh)

    def dev(self) -> arb:
        return up(self.lin() + self.e)

    def _co(self, o):
        return o if isinstance(o, AF) else self.const(o)

    def __add__(self, o):
        o = self._co(o)
        return self._new(self.c0 + o.c0, self.c1 + o.c1, self.c2 + o.c2, self.e + o.e)

    __radd__ = __add__

    def __neg__(self):
        return self._new(-self.c0, -self.c1, -self.c2, self.e)

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) - self

    def __mul__(self, o):
        if not isinstance(o, AF):
            o = acb(o)
            return self._new(self.c0 * o, self.c1 * o, self.c2 * o, self.e * abs(o).upper())
        e = abs(self.c0).upper() * o.e + abs(o.c0).upper() * self.e + self.dev() * o.dev()
        return self._new(self.c0 * o.c0, self.c0 * o.c1 + o.c0 * self.c1, self.c0 * o.c2 + o.c0 * self.c2, e)

    __rmul__ = __mul__

    def inv(self):
        """1/(c0 + X), |X| <= s < |c0|: 1/c0 - (c1 dw + c2 dh)/c0^2 + E,
        |E| <= e/|c0|^2 + s^2/(|c0|^2 (|c0| - s))."""
        a = abs(self.c0).lower()
        s = self.dev()
        if not (a > s):
            raise ZeroDivisionError("affine form may vanish")
        i0 = 1 / self.c0
        i2 = i0 * i0
        e = self.e / (a * a) + s * s / (a * a * (a - s))
        return self._new(i0, -self.c1 * i2, -self.c2 * i2, e)

    def __truediv__(self, o):
        if not isinstance(o, AF):
            return self * (1 / acb(o))
        return self * o.inv()

    def __rtruediv__(self, o):
        return self.inv() * o

    def absu(self) -> arb:
        return up(abs(self.c0).upper() + self.dev())

    def re_lower(self) -> arb:
        return arb((self.c0.real - self.dev()).lower())

    def re_upper(self) -> arb:
        return arb((self.c0.real + self.dev()).upper())

    def to_acb(self) -> acb:
        d = self.dev()
        return self.c0 + acb(arb(0, d), arb(0, d))
