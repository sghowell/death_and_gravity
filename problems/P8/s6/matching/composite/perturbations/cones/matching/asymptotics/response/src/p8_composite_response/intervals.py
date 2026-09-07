"""Independent rational outward intervals and first derivative jets.

Every endpoint is a Fraction. Rounding is OUTWARD to 32 binary places
after every operation; sqrt uses integer isqrt, never binary floating
point. This deliberately simple enclosure engine carries no dynamics
or prior gate conclusions.
"""

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt

BITS = 32
GRID = 1 << BITS


def exact(value):
    if not isinstance(value, bool) and isinstance(value, (int, Fraction)):
        return Fraction(value)
    raise TypeError("Intervals accept exact int/Fraction endpoints only")


def floor_grid(value):
    return Fraction((value.numerator*GRID)//value.denominator, GRID)


def ceil_grid(value):
    return -floor_grid(-value)


@dataclass(frozen=True)
class Interval:
    lo: Fraction
    hi: Fraction

    def __init__(self, lo, hi=None):
        low, high = exact(lo), exact(lo if hi is None else hi)
        if low > high:
            raise ValueError("Reversed interval")
        object.__setattr__(self, "lo", floor_grid(low))
        object.__setattr__(self, "hi", ceil_grid(high))

    @staticmethod
    def lift(value):
        return value if isinstance(value, Interval) else Interval(value)

    def __add__(self, other):
        other = self.lift(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        values = [x*y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return Interval(min(values), max(values))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ValueError("Interval denominator crosses zero")
        return Interval(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        return self*self.lift(other).reciprocal()

    def __rtruediv__(self, other):
        return self.lift(other)*self.reciprocal()

    def __pow__(self, power):
        if isinstance(power, bool) or not isinstance(power, int):
            raise TypeError("Only integer powers are supported")
        if power < 0:
            return self.reciprocal()**(-power)
        result = Interval(1)
        for _ in range(power):
            result = result*self
        return result

    def sqrt(self):
        if self.lo < 0:
            raise ValueError("Negative square root interval")
        low = isqrt((self.lo.numerator*GRID*GRID)//self.lo.denominator)
        high = isqrt((self.hi.numerator*GRID*GRID)//self.hi.denominator)
        if Fraction(high*high, GRID*GRID) < self.hi:
            high += 1
        return Interval(Fraction(low, GRID), Fraction(high, GRID))

    def magnitude(self):
        return max(abs(self.lo), abs(self.hi))

    def contains(self, value):
        return self.lo <= exact(value) <= self.hi


@dataclass(frozen=True)
class Jet:
    value: Interval
    gradient: tuple

    @classmethod
    def constant(cls, value):
        return cls(Interval.lift(value), (Interval(0),)*3)

    @classmethod
    def variable(cls, value, index):
        if isinstance(index, bool) or not isinstance(index, int):
            raise TypeError("A jet variable index must be an integer")
        if not 0 <= index < 3:
            raise ValueError("A jet variable index must be 0, 1 or 2")
        gradient = tuple(Interval(int(i == index)) for i in range(3))
        return cls(Interval.lift(value), gradient)

    @staticmethod
    def lift(value):
        return value if isinstance(value, Jet) else Jet.constant(value)

    def __add__(self, other):
        other = self.lift(other)
        return Jet(self.value+other.value, tuple(a+b for a, b in zip(self.gradient, other.gradient, strict=True)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, tuple(-a for a in self.gradient))

    def __sub__(self, other):
        return self+-self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        return Jet(self.value*other.value, tuple(a*other.value+self.value*b
                   for a, b in zip(self.gradient, other.gradient, strict=True)))

    __rmul__ = __mul__

    def reciprocal(self):
        value = self.value.reciprocal()
        return Jet(value, tuple(-a*value*value for a in self.gradient))

    def __truediv__(self, other):
        return self*self.lift(other).reciprocal()

    def __rtruediv__(self, other):
        return self.lift(other)*self.reciprocal()

    def __pow__(self, power):
        if isinstance(power, bool) or not isinstance(power, int):
            raise TypeError("Only integer powers are supported")
        if power < 0:
            return self.reciprocal()**(-power)
        result = Jet.constant(1)
        for _ in range(power):
            result = result*self
        return result

    def sqrt(self):
        value = self.value.sqrt()
        return Jet(value, tuple(a/(2*value) for a in self.gradient))

    def derivative_norm(self, indices=(0, 1, 2)):
        return sum((self.gradient[i].magnitude() for i in indices), Fraction(0))
