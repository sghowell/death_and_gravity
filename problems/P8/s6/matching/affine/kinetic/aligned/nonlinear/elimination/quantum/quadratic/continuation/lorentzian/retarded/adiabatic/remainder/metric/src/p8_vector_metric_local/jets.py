"""Exact physical lapse and log-scale source jets."""
from functools import cache

import sympy as sp
from p8_affine_kinetic import scalar
from p8_vector_mass_adiabatic import jets as parent

D, z = parent.D, parent.z
H = parent.H+sp.symbols("H_seventh H_eighth H_ninth", real=True)
alpha = parent.alpha+sp.symbols("alpha7:10", real=True)
beta = parent.beta+sp.symbols("beta7:10", real=True)
alpha2 = sp.symbols("a_NN0:10", real=True)
beta2 = sp.symbols("b_NN0:10", real=True)
n = parent.n+sp.symbols("source7:10", real=True)
v = sp.symbols("logscale_source0:10", real=True)
u = scalar.u
kind, order = parent.kind, parent.order


def output(value):
    if type(value) is not str or value not in ("N", "Z"):
        raise ValueError("Require physical lapse N or log-scale Z current")
    return value


def time(value):
    result = -2*H[0]*z*(1-z)*sp.diff(value, z)
    for row in (H, alpha, beta, alpha2, beta2, n, v):
        result += sum(sp.diff(value, row[j])*row[j+1] for j in range(len(row)-1))
    return sp.expand(result)


def weighted(value):
    return sp.expand(time(value)+D*H[0]*value)


@cache
def actual_mapping():
    r = 1+u**2
    functions = ((H, 4*u/r), (alpha, 4/(9*r**3)), (beta, 28/(81*r**3)),
                 (alpha2, -4/(3*r**3)-8/(9*r**6)),
                 (beta2, -28/(27*r**3)+152/(729*r**6)))
    return {field: sp.diff(function, u, j) for row, function in functions
            for j, field in enumerate(row)}


def actual(value):
    mapping = actual_mapping()
    fields = n+v
    polynomial = sp.Poly(sp.expand(value), *fields)
    return sp.Add(*(sp.factor(coefficient.subs(mapping))*sp.prod(field**power for field, power in zip(fields, powers))
                    for powers, coefficient in polynomial.terms()))


def linear_clean(value):
    fields = n+v
    expanded = sp.expand(value)
    if expanded == 0:
        return sp.Integer(0)
    polynomial = sp.Poly(expanded, *fields)
    if any(sum(powers) != 1 for powers, _ in polynomial.terms()):
        raise ValueError("Require a source-linear response")
    return sp.Add(*(sp.factor(expanded.coeff(field))*field for field in fields))
