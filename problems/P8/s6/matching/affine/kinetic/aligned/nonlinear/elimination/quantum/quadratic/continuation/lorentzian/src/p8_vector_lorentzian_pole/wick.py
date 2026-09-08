"""Local analytic continuation of jets, not a global Euclidean state."""
from functools import lru_cache

import sympy as sp
from p8_vector_quadratic import tensors
from p8_vector_quadratic_matching import kernel


def exact_expression(value):
    if type(value) is int:
        value = sp.Integer(value)
    if not isinstance(value, sp.Expr) or value.atoms(sp.Float) or value.has(sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("Require a finite exact symbolic local expression")
    return value


def continuation(value):
    value = exact_expression(value)
    mapping = {field: (-sp.I)**(order+1)*field for order, field in enumerate(tensors.H)}
    for mode in (tensors.plus, tensors.minus):
        for fields in (mode.temporal, mode.spatial):
            mapping.update({field: (-sp.I)**order*field for order, field in enumerate(fields)})
    return sp.expand(value.xreplace(mapping))


@lru_cache(maxsize=None, typed=True)
def loop_pole(order, jet=0):
    if type(order) is not int or order not in (2, 4) or type(jet) is not int or jet not in (0, 1):
        raise ValueError("Require derivative order 2/4 and dimensional Taylor jet 0/1")
    return -continuation(kernel.dimension_jet(order, jet))


def graded_actual(value, order):
    if type(order) is not int or order not in (0, 2, 4):
        raise ValueError("Require a supported even covariant derivative order")
    return sp.factor(-(-1)**(order//2)*exact_expression(value).subs(tensors.k, sp.I*tensors.k))
