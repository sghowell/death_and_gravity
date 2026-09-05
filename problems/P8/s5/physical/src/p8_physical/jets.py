"""Exact multilinear Fourier jets in labelled external legs.

epsilon_i^2=0 extracts a vertex with one copy of each labelled leg; it does
not drop repeated-field vertices (give those fields distinct leg labels).
No collinear/single-wave restriction. Coefficients live in QQ(i) or QQ(i)(x).
"""

from itertools import product

import sympy as sp
from sympy.polys.domains import QQ_I


class Context:
    def __init__(self, momenta, parameter=None):
        self.momenta = tuple(tuple(map(sp.Rational, k)) for k in momenta)
        if not self.momenta or any(len(k) != 3 for k in self.momenta):
            raise ValueError("Supply one or more three-dimensional momenta")
        self.n = len(self.momenta)
        self.full = (1 << self.n)-1
        self.domain = QQ_I if parameter is None else QQ_I.frac_field(parameter)
        self.zero_coefficient = self.domain.zero
        self.one_coefficient = self.domain.one
        self.wave = {mask: tuple(sum(self.momenta[j][i] for j in range(self.n) if mask & (1 << j))
                                for i in range(3)) for mask in range(self.full+1)}

    def convert(self, value):
        if self.domain.of_type(value):
            return value
        return self.domain.from_sympy(sp.sympify(value))

    def jet(self, value=0):
        return Jet(self, {0: self.convert(value)})

    def leg(self, index, value=1):
        return Jet(self, {1 << index: self.convert(value)})


class Jet:
    def __init__(self, context, data):
        self.context = context
        self.data = {k: v for k, v in data.items() if v != context.zero_coefficient}

    def _coerce(self, other):
        if isinstance(other, Jet):
            if other.context is not self.context:
                raise ValueError("Cannot mix Fourier contexts")
            return other
        return self.context.jet(other)

    def __add__(self, other):
        other = self._coerce(other)
        out = self.data.copy()
        for k, value in other.data.items():
            out[k] = out.get(k, self.context.zero_coefficient)+value
        return Jet(self.context, out)

    __radd__ = __add__

    def __neg__(self):
        return Jet(self.context, {k: -v for k, v in self.data.items()})

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) + (-self)

    def __mul__(self, other):
        other = self._coerce(other)
        out = {}
        zero = self.context.zero_coefficient
        for i, left in self.data.items():
            for j, right in other.data.items():
                if not i & j:
                    key = i | j
                    out[key] = out.get(key, zero)+left*right
        return Jet(self.context, out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Jet):
            return self*other.power(-1)
        return Jet(self.context, {k: v/self.context.convert(other) for k, v in self.data.items()})

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            return self.power(exponent)
        result = self.context.jet(1)
        for _ in range(exponent):
            result = result*self
        return result

    def power(self, exponent):
        constant = self.coefficient(0)
        if constant != 1:
            raise ValueError("Analytic powers currently require constant term one")
        out = self.context.jet(1)
        term = self.context.jet(1)
        remainder = self-1
        for k in range(1, self.context.n+1):
            term *= remainder
            out += sp.binomial(exponent, k)*term
        return out

    def derivative(self, axis):
        return Jet(self.context, {mask: value*self.context.convert(sp.I*self.context.wave[mask][axis])
                                  for mask, value in self.data.items()})

    def homogeneous(self, order):
        return Jet(self.context, {mask: value for mask, value in self.data.items() if mask.bit_count() == order})

    def coefficient(self, mask):
        return self.context.domain.to_sympy(self.data.get(mask, self.context.zero_coefficient))

    def is_zero(self):
        return not self.data


def zeros(context):
    return [[context.jet() for _ in range(3)] for _ in range(3)]


def identity(context):
    return [[context.jet(int(i == j)) for j in range(3)] for i in range(3)]


def madd(left, right):
    return [[left[i][j]+right[i][j] for j in range(3)] for i in range(3)]


def mscale(matrix, value):
    return [[entry*value for entry in row] for row in matrix]


def mmul(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def trace(matrix):
    return sum(matrix[i][i] for i in range(3))


def contract(left, right):
    return sum(left[i][j]*right[i][j] for i, j in product(range(3), repeat=2))


def determinant(matrix):
    a, b, c = matrix
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
