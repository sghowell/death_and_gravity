"""Fraction-only operator composition, endpoint jets and Hermite replay.

No SymPy or production imports. Sparse polynomials use powers (u,K).
The rational ring has the proved denominator factors d=1+u² and e=1-u².
"""

from fractions import Fraction as F
from functools import cache
from math import comb, factorial


def fraction(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError("Independent arithmetic admits exact integers/Fractions only")
    return F(value)


def poly(values):
    result = {}
    for powers, value in values.items():
        if (not isinstance(powers, tuple) or len(powers) != 2
                or any(isinstance(p, bool) or not isinstance(p, int) or p < 0 for p in powers)):
            raise ValueError("Require two nonnegative integer polynomial powers")
        value = fraction(value)
        if value:
            result[powers] = value
    return result


def add(left, right):
    result = dict(left)
    for power, value in right.items():
        result[power] = result.get(power, F(0))+value
    return poly(result)


def scale(value, coefficient):
    coefficient = fraction(coefficient)
    return poly({power: coefficient*item for power, item in value.items()})


def multiply(left, right):
    result = {}
    for (i, j), a in left.items():
        for (k, m), b in right.items():
            key = i+k, j+m
            result[key] = result.get(key, F(0))+a*b
    return poly(result)


def power(value, degree):
    result = {(0, 0): F(1)}
    for _ in range(degree):
        result = multiply(result, value)
    return result


def derivative(value):
    return {(i-1, j): i*coefficient for (i, j), coefficient in value.items() if i}


def divide_quadratic(value, leading):
    """Exact coefficient division by 1+leading*u², or None."""
    quotient = {}
    for k in {j for _, j in value}:
        remainder = {i: coefficient for (i, j), coefficient in value.items() if j == k}
        while remainder and max(remainder) >= 2:
            degree = max(remainder)
            coefficient = remainder[degree]/leading
            quotient[degree-2, k] = coefficient
            del remainder[degree]
            remainder[degree-2] = remainder.get(degree-2, F(0))-coefficient
            remainder = {i: v for i, v in remainder.items() if v}
        if remainder:
            return None
    return poly(quotient)


D_POLY = {(0, 0): F(1), (2, 0): F(1)}
E_POLY = {(0, 0): F(1), (2, 0): F(-1)}


class Rational:
    """Exact numerator divided by d^d_power e^e_power."""

    def __init__(self, numerator, d_power=0, e_power=0):
        if (isinstance(d_power, bool) or isinstance(e_power, bool)
                or not isinstance(d_power, int) or not isinstance(e_power, int)
                or min(d_power, e_power) < 0):
            raise ValueError("Require nonnegative integer denominator powers")
        numerator = poly(numerator)
        if not numerator:
            d_power = e_power = 0
        for name, leading in (("d", F(1)), ("e", F(-1))):
            count = d_power if name == "d" else e_power
            while count:
                quotient = divide_quadratic(numerator, leading)
                if quotient is None:
                    break
                numerator, count = quotient, count-1
            if name == "d":
                d_power = count
            else:
                e_power = count
        self.numerator, self.d_power, self.e_power = numerator, d_power, e_power

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Rational) else Rational({(0, 0): fraction(value)})

    def __add__(self, other):
        other = self.coerce(other)
        dp, ep = max(self.d_power, other.d_power), max(self.e_power, other.e_power)
        left = multiply(self.numerator, multiply(power(D_POLY, dp-self.d_power), power(E_POLY, ep-self.e_power)))
        right = multiply(other.numerator, multiply(power(D_POLY, dp-other.d_power), power(E_POLY, ep-other.e_power)))
        return Rational(add(left, right), dp, ep)

    __radd__ = __add__

    def __neg__(self):
        return Rational(scale(self.numerator, -1), self.d_power, self.e_power)

    def __sub__(self, other):
        return self+-self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other)+-self

    def __mul__(self, other):
        other = self.coerce(other)
        return Rational(multiply(self.numerator, other.numerator),
                        self.d_power+other.d_power, self.e_power+other.e_power)

    __rmul__ = __mul__

    def diff(self, order=1):
        if isinstance(order, bool) or not isinstance(order, int) or order < 0:
            raise ValueError("Require a nonnegative integer derivative order")
        value = self
        u = {(1, 0): F(1)}
        for _ in range(order):
            value = (Rational(derivative(value.numerator), value.d_power, value.e_power)
                     +Rational(scale(multiply(u, value.numerator), -2*value.d_power), value.d_power+1, value.e_power)
                     +Rational(scale(multiply(u, value.numerator), 2*value.e_power), value.d_power, value.e_power+1))
        return value

    def bound(self):
        radius = F(1, 50)
        numerator = sum(abs(v)*radius**i*4**j for (i, j), v in self.numerator.items())
        return numerator/(1-radius)**(2*self.e_power)

    def record(self):
        return {"numerator": self.numerator, "d_power": self.d_power, "e_power": self.e_power}


@cache
def source_coefficients():
    d, momentum = Rational(D_POLY), Rational({(0, 1): F(1)})
    d2, d6, d8 = Rational(power(D_POLY, 2)), Rational(power(D_POLY, 6)), Rational(power(D_POLY, 8))
    kg, kf, gg, gf = d6*F(1, 8), Rational({(0, 0): F(1, 2)}, 6), d2*F(1, 8), Rational({(0, 0): F(1, 2)}, 2)
    inverse = Rational(multiply(d8.numerator, add(power(d.numerator, 4), {(0, 0): F(-1)})), 0, 1)*F(1, 4)
    lf = [gf*momentum, kf.diff(), kf]
    lg = [gg*momentum, kg.diff(), kg]
    transform = [inverse*value for value in lf]
    transform[0] += 1
    # Literal composition Lg*(I+P*Lf)+Lf, independently using Leibniz.
    out = [Rational({}) for _ in range(5)]
    for i, left in enumerate(lg):
        for j, right in enumerate(transform):
            for r in range(i+1):
                out[j+r] += comb(i, r)*left*right.diff(i-r)
    for i, value in enumerate(lf):
        out[i] += value
    return tuple(value*Rational({(0, 0): F(4)}, 6) for value in out)


class Jet:
    """Ordinary Taylor coefficients of a scalar background, through degree5."""

    ORDER = 5

    def __init__(self, coefficients):
        self.values = tuple(fraction(coefficients[j]) if j < len(coefficients) else F(0)
                            for j in range(self.ORDER+1))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Jet) else Jet([value])

    def __add__(self, other):
        other = self.coerce(other)
        return Jet([a+b for a, b in zip(self.values, other.values)])

    __radd__ = __add__

    def __neg__(self):
        return Jet([-v for v in self.values])

    def __sub__(self, other):
        return self+-self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other)+-self

    def __mul__(self, other):
        other = self.coerce(other)
        return Jet([sum(self.values[k]*other.values[n-k] for k in range(n+1)) for n in range(self.ORDER+1)])

    __rmul__ = __mul__

    def inverse(self):
        if not self.values[0]:
            raise ZeroDivisionError("The punctured endpoint denominator is zero")
        out = [1/self.values[0]]
        for n in range(1, self.ORDER+1):
            out.append(-sum(self.values[k]*out[n-k] for k in range(1, n+1))/self.values[0])
        return Jet(out)

    def __truediv__(self, other):
        return self*self.coerce(other).inverse()

    def __rtruediv__(self, other):
        return self.coerce(other)*self.inverse()

    def __pow__(self, exponent):
        if isinstance(exponent, bool) or not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Require a nonnegative integer jet power")
        result = Jet([1])
        for _ in range(exponent):
            result = result*self
        return result

    def diff(self):
        return Jet([(j+1)*self.values[j+1] for j in range(self.ORDER)])


def matrix_add(left, right):
    return [[add(a, b) for a, b in zip(lrow, rrow)] for lrow, rrow in zip(left, right)]


def matrix_multiply(left, right):
    return [[sum_poly(multiply(left[i][k], right[k][j]) for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


def sum_poly(values):
    result = {}
    for value in values:
        result = add(result, value)
    return result


@cache
def endpoint_jets():
    ell = F(1, 100)
    u = Jet([-ell, ell])
    d = 1+u*u
    kg, kf, gg, gf = d**6/8, 1/(2*d**6), d**2/8, 1/(2*d**2)
    spring = 4*(1-u*u)/(d**8*(d**4-1))
    zero = Jet([0])
    pairs = [[(zero, zero) for _ in range(4)] for _ in range(4)]
    pairs[0][1] = pairs[2][3] = (Jet([1]), zero)
    pairs[1][0], pairs[1][1], pairs[1][2] = (-ell**2*spring/kg, -ell**2*gg/kg), (-kg.diff()/kg, zero), (ell**2*spring/kg, zero)
    pairs[3][0], pairs[3][2], pairs[3][3] = (ell**2*spring/kf, zero), (-ell**2*spring/kf, -ell**2*gf/kf), (-kf.diff()/kf, zero)
    derivatives = []
    for n in range(5):
        derivatives.append([[poly({(0, 0): factorial(n)*pair[0].values[n],
                                   (0, 1): factorial(n)*pair[1].values[n]}) for pair in row] for row in pairs])
    states = [[[{(0, 0): F(int(i == j))} if i == j else {} for j in range(4)] for i in range(4)]]
    for n in range(5):
        new = [[{} for _ in range(4)] for _ in range(4)]
        for j in range(n+1):
            product = matrix_multiply(derivatives[j], states[n-j])
            new = matrix_add(new, [[scale(entry, comb(n, j)) for entry in row] for row in product])
        states.append(new)
    return tuple(tuple(state[2]) for state in states)


def inverse_matrix(matrix):
    size = len(matrix)
    work = [[fraction(value) for value in row]+[F(int(i == j)) for j in range(size)]
            for i, row in enumerate(matrix)]
    for j in range(size):
        pivot = next(i for i in range(j, size) if work[i][j])
        work[j], work[pivot] = work[pivot], work[j]
        divisor = work[j][j]
        work[j] = [entry/divisor for entry in work[j]]
        for i in range(size):
            if i != j:
                factor = work[i][j]
                work[i] = [a-factor*b for a, b in zip(work[i], work[j])]
    return [row[size:] for row in work]


@cache
def hermite():
    # Solve the six right-end power-basis equations; the six left zeros
    # were imposed by f(x)=x^6 sum_(n=0)^5 a_n x^n.
    jets = endpoint_jets()
    matrix = [[F(factorial(6+n), factorial(6+n-j)) for n in range(6)] for j in range(6)]
    inverse = inverse_matrix(matrix)
    powers = [[{} for _ in range(4)] for _ in range(12)]
    for n in range(6):
        powers[6+n] = [sum_poly(scale(jets[j][column], inverse[n][j]) for j in range(6)) for column in range(4)]
    derivative_rows = []
    for order in range(7):
        degree = 11-order
        differentiated = [[scale(entry, F(factorial(n+order), factorial(n))) for entry in powers[n+order]]
                          for n in range(degree+1)]
        converted = [[sum_poly(scale(differentiated[n][column], F(comb(i, n), comb(degree, n)))
                               for n in range(i+1)) for column in range(4)] for i in range(degree+1)]
        derivative_rows.append(tuple(tuple(row) for row in converted))
    return {"power_rows": tuple(tuple(row) for row in powers), "derivative_rows": tuple(derivative_rows)}


def ceiling(value):
    return -((-value.numerator)//value.denominator)


@cache
def calibration():
    ell = F(1, 100)
    source = tuple(tuple(ceiling(value.diff(order).bound()/ell**j) for j, value in enumerate(source_coefficients()))
                   for order in range(3))
    f_jets = tuple(ceiling(max(sum(sum(abs(co)*4**j for (_, j), co in entry.items()) for entry in row)
                               for row in rows)) for rows in hermite()["derivative_rows"])
    cs = sum(source[0][j]*f_jets[j] for j in range(5))
    c2 = sum(comb(2, r)*source[2-r][j]*ell**(-r)*f_jets[j+r] for j in range(5) for r in range(3))
    return {"source_coefficient_bounds": source, "hermite_derivative_bounds": f_jets,
            "source_sup_per_scaled_endpoint_linf": cs, "source_uu_sup_per_scaled_endpoint_linf": c2,
            "source_L2_raw": F(4, 10)*cs, "source_uu_L2_raw": F(4, 10)*c2}
