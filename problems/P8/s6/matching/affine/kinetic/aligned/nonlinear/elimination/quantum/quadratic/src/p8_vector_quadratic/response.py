"""Actual clock-lapse insertion and compact-support self-adjoint local form."""
from functools import cache, lru_cache

import sympy as sp
from p8_affine_kinetic import scalar
from p8_vector_variation import readouts

from . import kernel, tensors

u = scalar.u
n = sp.symbols("clock_lapse0:9", real=True)


def time(value):
    flow = sp.diff(value, u)-scalar.background()["H"]*tensors.k*sp.diff(value, tensors.k)
    for order in range(len(n)-1):
        flow += n[order+1]*sp.diff(value, n[order])
    return sp.factor(flow)


def weighted_time(value):
    return sp.factor(time(value)+3*scalar.background()["H"]*value)


@cache
def mass_coefficients():
    h = scalar.background()["h"]
    return 4/(9*h), 28/(81*h)


@lru_cache(maxsize=None, typed=True)
def raw(order):
    if type(order) is not int or order not in (0, 2, 4):
        raise ValueError("Require native covariant derivative order 0,2 or 4")
    alpha, beta = mass_coefficients()
    mapping = {value: sp.diff(scalar.background()["H"], u, j) for j, value in enumerate(tensors.H)}
    for mode in (tensors.plus, tensors.minus):
        for data, coefficient in ((mode.temporal, alpha), (mode.spatial, beta)):
            for j, value in enumerate(data):
                mapping[value] = tensors.mass**2*sum(sp.binomial(j, r)*sp.diff(coefficient, u, j-r)*n[r] for r in range(j+1))
    # Use polynomial coefficient collection before rational simplification.
    substituted = sp.Poly(sp.expand(kernel.pole(order).subs(mapping)), *n[:5])
    return sp.Add(*(sp.factor(coefficient)*sp.prod(value**power for value, power in zip(n[:5], monomial))
                    for monomial, coefficient in substituted.terms()))


@lru_cache(maxsize=None, typed=True)
def local_form(order):
    if type(order) is not int or order not in (0, 2, 4):
        raise ValueError("Require native covariant derivative order 0,2 or 4")
    polynomial = sp.Poly(raw(order), *n[:3])
    c = lambda i, j, k: polynomial.coeff_monomial(n[0]**i*n[1]**j*n[2]**k)
    second = sp.factor(c(0, 0, 2))
    first = sp.factor(c(0, 2, 0)-c(1, 0, 1)-weighted_time(c(0, 1, 1))/2)
    zero = sp.factor(c(2, 0, 0)-weighted_time(c(1, 1, 0))/2
                     +weighted_time(weighted_time(c(1, 0, 1)))/2)
    return {"second_time_derivative_squared": second,
            "first_time_derivative_squared": first, "lapse_squared": zero}


@lru_cache(maxsize=None, typed=True)
def checks(order):
    if type(order) is not int or order not in (0, 2, 4):
        raise ValueError("Require native covariant derivative order 0,2 or 4")
    data = local_form(order)
    diagonal = data["second_time_derivative_squared"]*n[2]**2+data["first_time_derivative_squared"]*n[1]**2+data["lapse_squared"]*n[0]**2
    residual = 0
    for j in range(3):
        term = sp.diff(raw(order)-diagonal, n[j])
        for _ in range(j):
            term = -weighted_time(term)
        residual += term
    out = {"actual_lapse_compact_action_order_"+str(order): sp.factor(residual)}
    for key, value in closed_form(order).items():
        out["actual_lapse_closed_coefficient_"+str(order)+"_"+key] = sp.factor(data[key]-value)
    return out


@lru_cache(maxsize=None, typed=True)
def closed_form(order):
    if type(order) is not int or order not in (0, 2, 4):
        raise ValueError("Require native covariant derivative order 0,2 or 4")
    r, q, m = 1+u**2, tensors.k**2, tensors.mass
    if order == 0:
        coefficients = (0, 0, -sp.Rational(452, 2187)*m**4/r**6)
    elif order == 2:
        coefficients = (0, -sp.Rational(832, 6561)*m**2/r**6,
                        -sp.Rational(16, 2187)*m**2*(13*q*r**2+258*u**2+78)/r**8)
    else:
        coefficients = (-sp.Rational(128, 6561)/r**6,
                        -sp.Rational(64, 19683)*(11*q*r**2+612*u**2+90)/r**8,
                        -sp.Rational(16, 98415)*(101*q**2*r**4+10*q*r**2*(265*u**2+283)
                                                  +238230*u**4+172980*u**2-4770)/r**10)
    return dict(zip(("second_time_derivative_squared", "first_time_derivative_squared", "lapse_squared"),
                    map(sp.sympify, coefficients)))


@cache
def actual_mass_checks():
    alpha, beta = mass_coefficients()
    # S6.58's variation coefficients are the same lapse insertion used here.
    h = sp.Symbol("h", positive=True)
    data = readouts.actual_mass_jets()
    return {"temporal_profile_ratio": sp.factor(beta-sp.Rational(7, 9)*alpha),
            "frozen_temporal_lapse_mass_derivative": sp.factor(alpha-data["a"]["N_first"].subs(h, scalar.background()["h"])),
            "frozen_spatial_lapse_mass_derivative": sp.factor(beta-data["b"]["N_first"].subs(h, scalar.background()["h"]))}
