"""Small exact multivariate rational engine for generic symbol identities.

The native polynomial backend is python-flint, already used by the frozen
repository. It avoids treating checks at several rational times as a
generic-in-time proof. Coefficients are rational; no floating inputs or
numerical cancellation tolerance is used.
"""

from fractions import Fraction

import sympy as sp
from flint import fmpq, fmpq_mpoly, fmpq_mpoly_ctx

CONTEXT = fmpq_mpoly_ctx.get(("u", "c", "K"))


def polynomial(value):
    if isinstance(value, fmpq_mpoly):
        return value
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Integer, sp.Rational, fmpq)):
        raise TypeError("Only exact rational polynomial inputs are supported")
    value = Fraction(str(value))
    return CONTEXT.constant(fmpq(value.numerator, value.denominator))


class Rational:
    __slots__ = ("denominator", "numerator")

    def __init__(self, numerator=0, denominator=1):
        if isinstance(numerator, Rational) and denominator == 1:
            self.numerator, self.denominator = numerator.numerator, numerator.denominator
            return
        n, d = polynomial(numerator), polynomial(denominator)
        if not d:
            raise ZeroDivisionError("Zero rational-function denominator")
        if not n:
            self.numerator, self.denominator = polynomial(0), polynomial(1)
            return
        gcd = n.gcd(d)
        n, d = n/gcd, d/gcd
        scale = next(iter(d.to_dict().values()))
        self.numerator, self.denominator = n/scale, d/scale

    def __bool__(self):
        return bool(self.numerator)

    def __eq__(self, other):
        other = as_rational(other)
        return self.numerator*other.denominator == other.numerator*self.denominator

    def __neg__(self):
        return Rational(-self.numerator, self.denominator)

    def __add__(self, other):
        other = as_rational(other)
        if not self:
            return other
        if not other:
            return self
        gcd = self.denominator.gcd(other.denominator)
        left, right = self.denominator/gcd, other.denominator/gcd
        return Rational(self.numerator*right+other.numerator*left, left*other.denominator)

    __radd__ = __add__

    def __sub__(self, other):
        return self+-as_rational(other)

    def __rsub__(self, other):
        return as_rational(other)+-self

    def __mul__(self, other):
        other = as_rational(other)
        if not self or not other:
            return Rational(0)
        first = self.numerator.gcd(other.denominator)
        second = other.numerator.gcd(self.denominator)
        return Rational((self.numerator/first)*(other.numerator/second),
                        (self.denominator/second)*(other.denominator/first))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_rational(other)
        return self*Rational(other.denominator, other.numerator)

    def __rtruediv__(self, other):
        return as_rational(other)/self

    def __pow__(self, exponent):
        if isinstance(exponent, bool) or not isinstance(exponent, int):
            raise TypeError("Only exact integer exponents are supported")
        if exponent < 0:
            return Rational(self.denominator**(-exponent), self.numerator**(-exponent))
        return Rational(self.numerator**exponent, self.denominator**exponent)

    def derivative(self, index=0):
        if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 3:
            raise ValueError("Derivative index must be zero, one or two")
        return Rational(self.numerator.derivative(index)*self.denominator
                        -self.numerator*self.denominator.derivative(index), self.denominator**2)

    def degree(self, index=2):
        if not self:
            return -sp.oo
        return int(max(power[index] for power in self.numerator.to_dict())-max(
            power[index] for power in self.denominator.to_dict()))

    def leading(self, index=2):
        """Coefficient after removing the leading power of one variable."""
        if not self:
            return Rational(0)
        polynomials = []
        for p in (self.numerator, self.denominator):
            degree = max(power[index] for power in p.to_dict())
            terms = {}
            for power, coefficient in p.to_dict().items():
                if power[index] == degree:
                    exponent = list(power)
                    exponent[index] = 0
                    terms[tuple(exponent)] = coefficient
            polynomials.append(CONTEXT.from_dict(terms))
        return Rational(*polynomials)

    def sympy(self):
        variables = sp.symbols("u c K", real=True)
        def convert(p):
            return sum(sp.Rational(str(coefficient))*sp.prod(v**n for v, n in zip(variables, power, strict=True))
                       for power, coefficient in p.to_dict().items())
        return convert(self.numerator)/convert(self.denominator)

    def evaluate(self, values):
        if len(values) != 3 or any(isinstance(v, bool) or not isinstance(v, (int, Fraction)) for v in values):
            raise TypeError("Require three exact rational evaluation coordinates")
        values = tuple(map(Fraction, values))
        def evaluate_polynomial(p):
            total = Fraction(0)
            for power, coefficient in p.to_dict().items():
                term = Fraction(str(coefficient))
                for value, exponent in zip(values, power, strict=True):
                    term *= value**int(exponent)
                total += term
            return total
        return evaluate_polynomial(self.numerator)/evaluate_polynomial(self.denominator)


def as_rational(value):
    return value if isinstance(value, Rational) else Rational(value)


def from_sympy(value, variables):
    """Direct rational expression conversion, without calling cancel/factor."""
    value = sp.sympify(value)
    mapping = dict(zip(variables, (Rational(p) for p in CONTEXT.gens()), strict=True))
    def convert(expr):
        if expr in mapping:
            return mapping[expr]
        if expr.is_Rational:
            return Rational(expr)
        if expr.is_Add:
            return sum(convert(term) for term in expr.args)
        if expr.is_Mul:
            product = Rational(1)
            for term in expr.args:
                product *= convert(term)
            return product
        if expr.is_Pow and expr.exp.is_Integer:
            return convert(expr.base)**int(expr.exp)
        raise TypeError(f"Not a rational expression in the declared variables: {expr}")
    return convert(value)


def zeros(rows, columns):
    return [[Rational(0) for _ in range(columns)] for _ in range(rows)]


def transpose(matrix):
    return [list(row) for row in zip(*matrix, strict=True)]


def multiply(left, right):
    columns = transpose(right)
    return [[sum(x*y for x, y in zip(row, column, strict=True)) for column in columns] for row in left]


def add(left, right):
    return [[x+y for x, y in zip(row, other, strict=True)] for row, other in zip(left, right, strict=True)]


def differentiate(matrix):
    return [[x.derivative() for x in row] for row in matrix]


def inverse(matrix):
    size = len(matrix)
    augmented = [list(row)+[Rational(int(i == j)) for j in range(size)] for i, row in enumerate(matrix)]
    for j in range(size):
        pivot = next((i for i in range(j, size) if augmented[i][j]), None)
        if pivot is None:
            raise ValueError("Singular exact rational matrix")
        augmented[j], augmented[pivot] = augmented[pivot], augmented[j]
        coefficient = augmented[j][j]
        augmented[j] = [x/coefficient for x in augmented[j]]
        for i in range(size):
            if i != j and augmented[i][j]:
                coefficient = augmented[i][j]
                augmented[i] = [x-coefficient*y for x, y in zip(augmented[i], augmented[j], strict=True)]
    return [row[size:] for row in augmented]


def determinant(matrix):
    size = len(matrix)
    work = [list(row) for row in matrix]
    result = Rational(1)
    for j in range(size):
        pivot = next((i for i in range(j, size) if work[i][j]), None)
        if pivot is None:
            return Rational(0)
        if pivot != j:
            work[j], work[pivot] = work[pivot], work[j]
            result = -result
        coefficient = work[j][j]
        result *= coefficient
        for i in range(j+1, size):
            ratio = work[i][j]/coefficient
            for column in range(j+1, size):
                work[i][column] -= ratio*work[j][column]
    return result
