"""Exact outward rational enclosures with adaptive positive algebraic roots."""

from dataclasses import dataclass
from fractions import Fraction as F
from functools import cache

import sympy as s


def fraction(value):
    if isinstance(value, (bool, float, s.Float)):
        raise TypeError("Require exact rational input, not boolean or float")
    if not isinstance(value, (int, F, s.Rational)):
        raise TypeError("Require an integer or exact rational value")
    if isinstance(value, F):
        return value
    value = s.Rational(value)
    return F(int(value.p), int(value.q))


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    def __init__(self, lo, hi=None):
        object.__setattr__(self, "lo", fraction(lo))
        object.__setattr__(self, "hi", fraction(lo if hi is None else hi))
        assert self.lo <= self.hi

    def __add__(self, other):
        other = box(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -box(other)

    def __rsub__(self, other):
        return box(other) + -self

    def __mul__(self, other):
        other = box(other)
        values = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return I(min(values), max(values))

    __rmul__ = __mul__

    def inverse(self):
        assert self.lo > 0 or self.hi < 0, "interval inverse crosses zero"
        return I(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        return self * box(other).inverse()

    def __rtruediv__(self, other):
        return box(other) * self.inverse()

    def __pow__(self, power):
        power = fraction(power)
        if power.denominator == 1:
            n = power.numerator
            if n < 0:
                return self.inverse() ** (-n)
            if n == 0:
                return I(1)
            if n % 2 == 0 and self.lo <= 0 <= self.hi:
                return I(0, max(abs(self.lo), abs(self.hi)) ** n)
            values = [self.lo**n, self.hi**n]
            return I(min(values), max(values))
        assert self.lo > 0, "fractional power requires positive interval"
        n = power.denominator

        def root_bound(value):
            exponent = (
                value.numerator.bit_length() - value.denominator.bit_length()
            ) // n
            scale = 2 ** max(100, 100 - exponent)
            scaled = value * scale**n
            lower = int(s.integer_nthroot(scaled.numerator // scaled.denominator, n)[0])
            return F(lower, scale), F(lower + 1, scale)

        left = root_bound(self.lo)[0]
        right = root_bound(self.hi)[1]
        return I(left, right) ** power.numerator

    def inside(self, lo, hi):
        return fraction(lo) <= self.lo <= self.hi <= fraction(hi)

    def bounds(self):
        return [s.Rational(v.numerator, v.denominator) for v in (self.lo, self.hi)]


def box(value):
    return value if isinstance(value, I) else I(value)


def evaluate(expression, bindings):
    @cache
    def visit(expr):
        if expr in bindings:
            return box(bindings[expr])
        if expr.is_Rational:
            return I(expr)
        if expr.is_Add:
            return sum((visit(v) for v in expr.args), I(0))
        if expr.is_Mul:
            result = I(1)
            for term in expr.args:
                result = result * visit(term)
            return result
        if expr.is_Pow and expr.exp.is_Rational:
            return visit(expr.base) ** expr.exp
        raise TypeError("Unsupported exact interval expression: " + str(expr))

    return visit(s.sympify(expression))
