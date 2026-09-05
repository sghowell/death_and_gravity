"""Finite ordinary Taylor jets with exact SymPy polynomial coefficients."""

import sympy as sp

ORDER = 4


def constant(value):
    return (sp.sympify(value),)+(sp.Integer(0),)*ORDER


def add(*rows):
    return tuple(sp.expand(sum(row[k] for row in rows)) for k in range(ORDER+1))


def scale(row, value):
    return tuple(sp.expand(value*entry) for entry in row)


def mul(*rows):
    result = constant(1)
    for row in rows:
        result = tuple(sp.expand(sum(result[j]*row[k-j] for j in range(k+1))) for k in range(ORDER+1))
    return result


def power(row, exponent):
    if row[0] != 1:
        raise ValueError("This binomial jet requires constant term exactly one")
    h = (sp.Integer(0), *row[1:])
    result, term = constant(1), constant(1)
    for j in range(1, ORDER+1):
        term = mul(term, h)
        result = add(result, scale(term, sp.binomial(exponent, j)))
    return result


def Npower(exponent):
    return tuple(sp.binomial(exponent, k) for k in range(ORDER+1))


def derivative_coefficients(row):
    return tuple(sp.factor(value*sp.factorial(k)) for k, value in enumerate(row))


def identities():
    n, exponent = sp.symbols("n exponent", real=True)
    ts = sp.symbols("t1:5", real=True)
    row = (sp.Integer(1), *ts)
    p = power(row, exponent)
    T = sum(value*n**k for k, value in enumerate(row))
    F = sum(value*n**k for k, value in enumerate(p))
    residual = sp.Poly(sp.expand(T*sp.diff(F, n)-exponent*sp.diff(T, n)*F), n)
    return {f"binomial_differential_order_{k}": sp.simplify(residual.nth(k)) for k in range(ORDER)}
