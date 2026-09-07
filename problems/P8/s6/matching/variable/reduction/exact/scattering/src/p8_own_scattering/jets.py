"""Exact interval Taylor jets in two variables through total degree three.

Coefficients are derivative/factorial intervals, not sampled derivatives or
a finite-step numerical approximation.  Inversion and the positive cube
root use nilpotent binomial identities in the truncated polynomial ring.
"""

from collections.abc import Mapping
from fractions import Fraction
from math import factorial

from p8_exact_stationary.intervals import Interval, as_interval

DEGREE = 3
INDICES = tuple((i, j) for i in range(DEGREE + 1) for j in range(DEGREE + 1 - i))


def _index(i, j):
    if type(i) is not int or type(j) is not int:
        raise TypeError("jet indices must be exact integers")
    if i < 0 or j < 0 or i + j > DEGREE:
        raise ValueError("jet indices must have total degree at most three")
    return i, j


class TaylorJet:
    """Immutable interval coefficients in increments of (v,zeta)."""

    __slots__ = ("_coefficients",)

    def __init__(self, values):
        if not isinstance(values, Mapping):
            raise TypeError("a jet must be constructed from a coefficient mapping")
        for key in values:
            if not isinstance(key, tuple) or len(key) != 2:
                raise TypeError("a jet index is a pair of integers")
            _index(*key)
        object.__setattr__(self, "_coefficients",
                           tuple(as_interval(values.get(index, 0)) for index in INDICES))

    def __setattr__(self, name, value):
        raise AttributeError("TaylorJet coefficients are immutable")

    @classmethod
    def constant(cls, value):
        return cls({(0, 0): as_interval(value)})

    @classmethod
    def variable(cls, axis, value):
        if type(axis) is not int:
            raise TypeError("the variable axis must be an exact integer")
        if axis not in (0, 1):
            raise ValueError("the only axes are v=0 and zeta=1")
        return cls({(0, 0): as_interval(value), (1 - axis, axis): Interval(1)})

    def coefficient(self, i, j):
        return self._coefficients[INDICES.index(_index(i, j))]

    def partial(self, i, j):
        return self.coefficient(i, j)*factorial(i)*factorial(j)

    def coefficients(self):
        return dict(zip(INDICES, self._coefficients, strict=True))

    def __add__(self, other):
        other = as_jet(other)
        return TaylorJet({key: self.coefficient(*key) + other.coefficient(*key) for key in INDICES})

    __radd__ = __add__

    def __neg__(self):
        return TaylorJet({key: -value for key, value in self.coefficients().items()})

    def __sub__(self, other):
        return self + -as_jet(other)

    def __rsub__(self, other):
        return as_jet(other) + -self

    def __mul__(self, other):
        other = as_jet(other)
        result = {}
        for (i, j), left in self.coefficients().items():
            for (p, q), right in other.coefficients().items():
                key = i + p, j + q
                if sum(key) <= DEGREE:
                    result[key] = result.get(key, Interval(0)) + left*right
        return TaylorJet(result)

    __rmul__ = __mul__

    def _binomial(self, exponent, leading):
        base = self.coefficient(0, 0)
        # Remove the constant formally.  Subtracting two equal interval
        # constants would incorrectly introduce a nonzero constant width.
        h = TaylorJet({key: value/base for key, value in self.coefficients().items()
                       if key != (0, 0)})
        total, power, coefficient = TaylorJet.constant(1), TaylorJet.constant(1), Fraction(1)
        for n in range(1, DEGREE + 1):
            coefficient *= (exponent - n + 1)/n
            power = power*h
            total = total + power*coefficient
        return total*leading

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("jet powers require an integer; use cube_root with a proved bracket")
        if exponent < 0:
            return self._binomial(Fraction(exponent), self.coefficient(0, 0)**exponent)
        result = TaylorJet.constant(1)
        for _ in range(exponent):
            result = result*self
        return result

    def inverse(self):
        return self**(-1)

    def cube_root(self, bracket):
        """Positive root, accepted only after an exact cubed-endpoint proof."""
        bracket = as_interval(bracket)
        base = self.coefficient(0, 0)
        if base.lo <= 0 or bracket.lo <= 0:
            raise ValueError("a positive real cube-root jet needs positive base and bracket")
        if bracket.lo**3 > base.lo or bracket.hi**3 < base.hi:
            raise ValueError("the proposed rational cube-root bracket was not proved")
        return self._binomial(Fraction(1, 3), bracket)

    def __truediv__(self, other):
        return self*as_jet(other).inverse()

    def __rtruediv__(self, other):
        return as_jet(other)*self.inverse()


def as_jet(value):
    return value if isinstance(value, TaylorJet) else TaylorJet.constant(value)


def constant(value):
    return TaylorJet.constant(value)


def variable(axis, value):
    return TaylorJet.variable(axis, value)
