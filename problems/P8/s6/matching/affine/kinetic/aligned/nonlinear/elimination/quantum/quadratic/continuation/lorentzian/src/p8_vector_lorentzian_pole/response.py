"""Physical clock projection and variation of the Lorentzian local pole."""
from functools import lru_cache

import sympy as sp
from p8_vector_quadratic import response as clock
from p8_vector_quadratic import tensors
from p8_vector_quadratic_matching import response as euclidean
from p8_vector_quadratic_matching.geometry import dimension

from . import bimetric, frame, wick

n, u, H = euclidean.n, euclidean.u, euclidean.H


def valid(order, jet):
    euclidean.valid(order, jet)


@lru_cache(maxsize=None, typed=True)
def raw(order, jet):
    valid(order, jet)
    generic = frame.second_loop_pole() if order == 2 else bimetric.loop_pole()
    value = sp.expand(sp.diff(generic, dimension, jet).subs(dimension, 3))
    fields = tensors.plus.temporal+tensors.plus.spatial+tensors.minus.temporal+tensors.minus.spatial
    mapping = {}
    for mode in (tensors.plus, tensors.minus):
        for data, coefficient in zip((mode.temporal, mode.spatial), clock.mass_coefficients()):
            for j, field in enumerate(data):
                mapping[field] = tensors.mass**2*sum(sp.binomial(j, r)*sp.diff(coefficient, u, j-r)*n[r] for r in range(j+1))
    hmap = {field: sp.diff(H, u, j) for j, field in enumerate(tensors.H)}
    coefficients = {}
    for powers, coefficient in sp.Poly(value, *fields).terms():
        factors = [field for field, power in zip(fields, powers) for _ in range(power)]
        if len(factors) != 2:
            raise ValueError("Require a quadratic mass insertion")
        projected = sp.Poly(sp.expand(mapping[factors[0]]*mapping[factors[1]]), *n[:3])
        for monomial, field_coefficient in projected.terms():
            coefficients.setdefault(monomial, []).append(coefficient.subs(hmap)*field_coefficient)
    return sp.Add(*(sp.factor(sum(terms))*sp.prod(field**power for field, power in zip(n[:3], monomial))
                    for monomial, terms in coefficients.items()))


@lru_cache(maxsize=None, typed=True)
def compact(order, jet):
    valid(order, jet)
    old, current = raw(order, 0), raw(order, jet)
    c = lambda i, j, k: euclidean.coefficient(current, i, j, k)
    second = c(0, 0, 2)
    first = c(0, 2, 0)-c(1, 0, 1)-clock.weighted_time(c(0, 1, 1))/2
    zero = c(2, 0, 0)-clock.weighted_time(c(1, 1, 0))/2+clock.weighted_time(clock.weighted_time(c(1, 0, 1)))/2
    if jet == 1:
        first -= H*euclidean.coefficient(old, 0, 1, 1)/2
        mixed = euclidean.coefficient(old, 1, 0, 1)
        zero += (-H*euclidean.coefficient(old, 1, 1, 0)+clock.weighted_time(H*mixed)+H*clock.weighted_time(mixed))/2
    return {key: sp.factor(value) for key, value in
            zip(("second_time_derivative_squared", "first_time_derivative_squared", "lapse_squared"),
                (second, first, zero))}


@lru_cache(maxsize=None, typed=True)
def operator(order, jet):
    valid(order, jet)
    value = euclidean.ordinary_operator(compact(order, jet))
    if jet == 1:
        old = compact(order, 0)
        value += (-2*H*old["first_time_derivative_squared"]*n[1]
                  +2*clock.weighted_time(H*old["second_time_derivative_squared"]*n[2])
                  +2*H*clock.weighted_time(old["second_time_derivative_squared"]*n[2]))
    return sp.factor(value)


@lru_cache(maxsize=None, typed=True)
def checks(order):
    valid(order, 0)
    out = {}
    for jet in (0, 1):
        out["Lorentzian_raw_"+str(order)+"_"+str(jet)] = sp.factor(raw(order, jet)-wick.graded_actual(euclidean.raw_jet(order, jet), order))
        for key, value in compact(order, jet).items():
            out["Lorentzian_compact_"+str(order)+"_"+str(jet)+"_"+key] = sp.factor(
                value-wick.graded_actual(euclidean.compact_jet(order, jet)[key], order))
        out["Lorentzian_operator_"+str(order)+"_"+str(jet)] = sp.factor(
            operator(order, jet)-wick.graded_actual(euclidean.operator_jet(order, jet), order))
    return out


@lru_cache(maxsize=None, typed=True)
def finite_counterterm(order):
    valid(order, 1)
    return 2*operator(order, 1)
