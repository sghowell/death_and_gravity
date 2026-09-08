"""Actual lapse counterterm variation before the dimensional limit."""
from functools import lru_cache

import sympy as sp
from p8_affine_kinetic import scalar
from p8_vector_quadratic import response as four
from p8_vector_quadratic import tensors

from . import kernel

n, u = four.n, four.u
H = scalar.background()["H"]


def valid(order, jet):
    if type(order) is not int or order not in (2, 4) or type(jet) is not int or jet not in (0, 1):
        raise ValueError("Require derivative order 2/4 and dimensional Taylor jet 0/1")


@lru_cache(maxsize=None, typed=True)
def raw_jet(order, jet):
    valid(order, jet)
    mass_fields = tensors.plus.temporal+tensors.plus.spatial+tensors.minus.temporal+tensors.minus.spatial
    mapping = {}
    for mode in (tensors.plus, tensors.minus):
        for fields, coefficient in zip((mode.temporal, mode.spatial), four.mass_coefficients()):
            for j, value in enumerate(fields):
                mapping[value] = tensors.mass**2*sum(sp.binomial(j, r)*sp.diff(coefficient, u, j-r)*n[r] for r in range(j+1))
    hmap = {value: sp.diff(H, u, j) for j, value in enumerate(tensors.H)}
    coefficients = {}
    for powers, coefficient in sp.Poly(kernel.dimension_jet(order, jet), *mass_fields).terms():
        factors = [field for field, power in zip(mass_fields, powers) for _ in range(power)]
        if len(factors) != 2:
            raise ValueError("Require a quadratic mass insertion")
        field_product = sp.Poly(sp.expand(mapping[factors[0]]*mapping[factors[1]]), *n[:3])
        for monomial, field_coefficient in field_product.terms():
            coefficients.setdefault(monomial, []).append(coefficient.subs(hmap)*field_coefficient)
    return sp.Add(*(sp.factor(sum(terms))*sp.prod(value**power for value, power in zip(n[:3], monomial))
                    for monomial, terms in coefficients.items()))


def coefficient(raw, i, j, k):
    return sp.Poly(raw, *n[:3]).coeff_monomial(n[0]**i*n[1]**j*n[2]**k)


@lru_cache(maxsize=None, typed=True)
def compact_jet(order, jet):
    valid(order, jet)
    old, current = raw_jet(order, 0), raw_jet(order, jet)
    c = lambda i, j, k: coefficient(current, i, j, k)
    second = c(0, 0, 2)
    first = c(0, 2, 0)-c(1, 0, 1)-four.weighted_time(c(0, 1, 1))/2
    zero = c(2, 0, 0)-four.weighted_time(c(1, 1, 0))/2+four.weighted_time(four.weighted_time(c(1, 0, 1)))/2
    if jet == 1:
        first -= H*coefficient(old, 0, 1, 1)/2
        mixed = coefficient(old, 1, 0, 1)
        zero += (-H*coefficient(old, 1, 1, 0)+four.weighted_time(H*mixed)+H*four.weighted_time(mixed))/2
    return {key: sp.factor(value) for key, value in
            zip(("second_time_derivative_squared", "first_time_derivative_squared", "lapse_squared"),
                (second, first, zero))}


def ordinary_operator(coefficients):
    return (2*coefficients["lapse_squared"]*n[0]
            -2*four.weighted_time(coefficients["first_time_derivative_squared"]*n[1])
            +2*four.weighted_time(four.weighted_time(coefficients["second_time_derivative_squared"]*n[2])))


@lru_cache(maxsize=None, typed=True)
def operator_jet(order, jet):
    valid(order, jet)
    value = ordinary_operator(compact_jet(order, jet))
    if jet == 1:
        old = compact_jet(order, 0)
        value += (-2*H*old["first_time_derivative_squared"]*n[1]
                  +2*four.weighted_time(H*old["second_time_derivative_squared"]*n[2])
                  +2*H*four.weighted_time(old["second_time_derivative_squared"]*n[2]))
    return sp.factor(value)


@lru_cache(maxsize=None, typed=True)
def checks(order):
    valid(order, 0)
    out = {"actual_dimension_zero_raw_"+str(order): sp.factor(raw_jet(order, 0)-four.raw(order))}
    for key, value in compact_jet(order, 0).items():
        out["actual_dimension_zero_compact_"+str(order)+"_"+key] = sp.factor(value-four.local_form(order)[key])
    for jet in (0, 1):
        current, old = raw_jet(order, jet), raw_jet(order, 0)
        direct = sp.diff(current, n[0])-four.weighted_time(sp.diff(current, n[1]))
        direct += four.weighted_time(four.weighted_time(sp.diff(current, n[2])))
        if jet == 1:
            direct += (-H*sp.diff(old, n[1])+four.weighted_time(H*sp.diff(old, n[2]))
                       +H*four.weighted_time(sp.diff(old, n[2])))
        out["vary_before_dimension_limit_order_"+str(order)+"_jet_"+str(jet)] = sp.factor(direct-operator_jet(order, jet))
    out.update({"actual_first_dimension_jet_"+str(order)+"_"+key: sp.factor(compact_jet(order, 1)[key]-value)
                for key, value in closed_first_jet(order).items()})
    return out


@lru_cache(maxsize=None, typed=True)
def closed_first_jet(order):
    valid(order, 1)
    r, q, m = 1+u**2, tensors.k**2, tensors.mass
    if order == 2:
        values = (0, -sp.Rational(574, 19683)*m**2/r**6,
                  -sp.Rational(2, 19683)*m**2*(259*q*r**2+32994*u**2+2106)/r**8)
    else:
        values = (-sp.Rational(2896, 98415)/r**6,
                  -sp.Rational(4, 98415)*(1389*q*r**2+94370*u**2+8070)/r**8,
                  -sp.Rational(4, 295245)*(1995*q**2*r**4+q*r**2*(149764*u**2+49984)
                                           +6795480*u**4+3319008*u**2-107352)/r**10)
    return dict(zip(("second_time_derivative_squared", "first_time_derivative_squared", "lapse_squared"),
                    map(sp.sympify, values)))


@lru_cache(maxsize=None, typed=True)
def finite_counterterm_operator(order):
    valid(order, 1)
    # The counterterm is -Q_D/(32*pi^2*epsilon), D=3-2epsilon.
    # Its normalized lapse variation is -E0/(32*pi^2*epsilon)
    # +2*(partial_D E_D)_3/(32*pi^2). This is only its
    # evanescent component, not the combined finite quantum kernel.
    return 2*operator_jet(order, 1)


@lru_cache(maxsize=None, typed=True)
def measure_control(order):
    valid(order, 1)
    naive = ordinary_operator(compact_jet(order, 1))
    omitted = sp.factor(operator_jet(order, 1)-naive)
    fixture = {value: 0 for value in n}
    fixture.update({u: 0, tensors.k: 0, tensors.mass: 1000, n[2]: 1})
    return {"omitted_measure_variation": omitted,
            "bounce_fixture": sp.factor(omitted.subs(fixture))}
