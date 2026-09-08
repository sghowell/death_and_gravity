"""Variational consistency of the full two-current local response."""
from functools import cache

import sympy as sp
from p8_aligned_quantum import potential
from p8_vector_dimensional import local

from . import counterterms, jets, radial


def W(value):
    return sp.factor(sp.diff(value, jets.u)+12*jets.u*value/(1+jets.u**2))


def coefficients(output, order, source):
    output, order, source = jets.output(output), jets.order(order), jets.output(source)
    value = sp.expand(radial.matched(output, order))
    row = jets.n if source == "N" else jets.v
    return tuple(sp.factor(value.coeff(field)) for field in row[:5])


def adjoint(row):
    result = [sp.Integer(0) for _ in row]
    for j, coefficient in enumerate(row):
        for k in range(j+1):
            value = coefficient
            for _ in range(j-k):
                value = W(value)
            result[k] += (-1)**j*sp.binomial(j, k)*value
    return tuple(sp.factor(value) for value in result)


@cache
def self_adjoint_checks():
    out = {}
    for order in (0, 1, 2):
        for output, source in (("N", "N"), ("N", "Z"), ("Z", "Z")):
            left = coefficients(output, order, source)
            right = adjoint(coefficients(source, order, output))
            for j in range(5):
                out[output+source+"_finite_adjoint_"+str(2*order)+"_"+str(j)] = sp.factor(left[j]-right[j])
        for output in ("N", "Z"):
            value = radial.matched(output, order)
            out[output+"_no_higher_source_jet_"+str(2*order)] = sp.expand(
                value-sum(sp.diff(value, field)*field for field in jets.n[:5]+jets.v[:5]))
    return out


@cache
def potential_checks():
    e, n, v = counterterms.e, jets.n[0], jets.v[0]
    am = 1+e*jets.alpha[0]*n+e**2*jets.alpha2[0]*n**2/2
    bm = 1+e*jets.beta[0]*n+e**2*jets.beta2[0]*n**2/2
    B = potential.coefficients()["finite_weight"].subs(potential.Lscale, 0)
    F = -(1+e*n)*sp.exp(3*e*v)*B.subs({potential.a: am, potential.b: bm}, simultaneous=True)
    quadratic = sp.diff(F, e, 2).subs(e, 0)/2
    return {output+"_full_finite_potential": sp.factor(radial.matched(output, 0)
                -counterterms.mass**4*jets.actual(sp.diff(quadratic, field)))
            for output, field in (("N", n), ("Z", v))}


@cache
def ordinary_checks():
    out = {}
    zero = {field: 0 for row in (jets.alpha, jets.beta, jets.alpha2, jets.beta2) for field in row}
    e, D = counterterms.e, jets.D
    volume = 1+e*(jets.n[0]+D*jets.v[0])+e**2*(D*jets.n[0]*jets.v[0]+D**2*jets.v[0]**2/2)
    geometry = counterterms.geometry()
    for order in (0, 1, 2):
        invariant = {0: sp.Rational(5, 2), 1: 5*geometry["R"]/3, 2: -4*geometry["a4sc"]}[order]
        density = counterterms.keep(volume*invariant).coeff(e, 2).subs(D, 3)
        for output, row in (("N", jets.n), ("Z", jets.v)):
            target = sp.Integer(0)
            for j, field in enumerate(row):
                term = sp.diff(density, field)
                for _ in range(j):
                    term = -jets.weighted(term).subs(D, 3)
                target += term
            raw = radial.laurent(output, order)["finite_MSbar_mu_m"].subs(zero)
            extra = 2*sp.diff(counterterms.generic_operator(output, order), D).subs(D, 3).subs(zero)
            out[output+"_ordinary_independent_finite_heat_action_"+str(2*order)] = sp.expand(raw+extra-target)
    return out


@cache
def scale_checks():
    out = {}
    finite = counterterms.first.finite_coefficients()
    for order in (0, 1, 2):
        for output, component, factor in (("N", "energy", -3), ("Z", "pressure", 9)):
            target = factor*counterterms.mass**(4-2*order)*jets.actual(finite[order][component].subs(local.ell, 0))
            out[output+"_constant_spatial_rescaling_"+str(2*order)] = sp.factor(
                coefficients(output, order, "Z")[0]-target)
    return out
