"""Unintegrated time response to full radial Laurent and counterterm jets."""
from functools import cache

import sympy as sp
from p8_affine_kinetic import scalar
from p8_vector_dimensional import local
from p8_vector_lorentzian_pole import response as counterterm
from p8_vector_quadratic import tensors

from . import jets, variation

mass = tensors.mass
u = scalar.u


def laurent(order):
    order = jets.order(order)
    return _laurent(order)


@cache
def _laurent(order):
    coefficient = variation.combined(order)
    radial = local.radial_polynomial(coefficient, order)
    factor = sp.Rational((-1)**(2-order), sp.factorial(2-order))
    pole = sp.expand(factor*radial.subs(jets.D, 3))
    evanescent = sp.expand(-2*factor*sp.diff(radial, jets.D).subs(jets.D, 3))
    return {"pole": pole, "evanescent": evanescent,
            "finite_MSbar_mu_m": sp.expand(sp.harmonic(2-order)*pole+evanescent)}


def actual(value):
    H = scalar.background()["H"]
    alpha, beta = 4/(9*(1+u**2)**3), 28/(81*(1+u**2)**3)
    mapping = {}
    for row, function in ((jets.H, H), (jets.alpha, alpha), (jets.beta, beta)):
        mapping.update({field: sp.diff(function, u, j) for j, field in enumerate(row)})
    polynomial = sp.Poly(sp.expand(value), *jets.n)
    return sp.Add(*(sp.factor(coefficient.subs(mapping))*sp.prod(field**power for field, power in zip(jets.n, powers))
                    for powers, coefficient in polynomial.terms()))


def matched(order):
    order = jets.order(order)
    return _matched(order)


@cache
def _matched(order):
    data = laurent(order)
    finite = mass**(4-2*order)*actual(data["finite_MSbar_mu_m"])
    if order:
        mapping = {field: jets.n[j] for j, field in enumerate(counterterm.n) if j < len(jets.n)}
        extra = counterterm.finite_counterterm(2*order).subs(mapping).subs(tensors.k, 0)
        finite += 2*extra
    expanded = sp.expand(finite)
    return sp.Add(*(sp.factor(expanded.coeff(field))*field for field in jets.n))


@cache
def pole_checks():
    out = {}
    for j in (0, 1, 2):
        actual_pole = mass**(4-2*j)*actual(laurent(j)["pole"])
        if j:
            mapping = {field: jets.n[k] for k, field in enumerate(counterterm.n) if k < len(jets.n)}
            target = 2*counterterm.operator(2*j, 0).subs(mapping).subs(tensors.k, 0)
        else:
            target = sp.Rational(1808, 2187)*mass**4*jets.n[0]/(1+u**2)**6
        out["homogeneous_curved_response_pole_"+str(2*j)] = sp.factor(actual_pole-target)
    return out
