"""Finite homogeneous local response in a self-adjoint quadratic form."""
from functools import cache

import sympy as sp
from p8_aligned_quantum import potential
from p8_vector_dimensional import local

from . import jets, radial


def W(value):
    return sp.factor(sp.diff(value, radial.u)+3*radial.scalar.background()["H"]*value)


def coefficients(order):
    return _coefficients(jets.order(order))


@cache
def _coefficients(order):
    operator = sp.expand(radial.matched(order))
    c = [sp.factor(operator.coeff(field)) for field in jets.n]
    A = c[4]/2
    B = W(W(A))-c[2]/2
    C = c[0]/2
    return {"second_time_square": sp.factor(A), "first_time_square": sp.factor(B), "source_square": sp.factor(C)}


def reconstructed(order):
    data = coefficients(order)
    A, B, C = [data[key] for key in ("second_time_square", "first_time_square", "source_square")]
    H, Hprime = radial.scalar.background()["H"], sp.diff(radial.scalar.background()["H"], radial.u)
    u, n = radial.u, jets.n
    return (2*C*n[0]-2*W(B)*n[1]-2*B*n[2]
            +2*(sp.diff(A, u, 2)+6*H*sp.diff(A, u)+(3*Hprime+9*H**2)*A)*n[2]
            +2*(2*sp.diff(A, u)+6*H*A)*n[3]+2*A*n[4])


@cache
def checks():
    out = {"finite_self_adjoint_local_form_"+str(2*j): sp.factor(radial.matched(j)-reconstructed(j))
           for j in (0, 1, 2)}
    source = sp.Symbol("flat_source", real=True)
    mapping = {potential.a: 1+jets.alpha[0]*source, potential.b: 1+jets.beta[0]*source}
    finite_potential = potential.coefficients()["finite_weight"].subs(potential.Lscale, 0).subs(mapping)
    target = -sp.diff(finite_potential, source, 2).subs(source, 0)*jets.n[0]
    out["finite_static_response_matches_frozen_full_potential"] = sp.factor(radial.laurent(0)["finite_MSbar_mu_m"]-target)
    for j in (0, 1, 2):
        angular = local.radial_polynomial(radial.variation.combined(j), j)
        factor = sp.Rational((-1)**(2-j), sp.factorial(2-j))
        # An early D=3 limit of the coefficient drops its contribution to the finite term.
        frozen = local.radial_polynomial(radial.variation.combined(j).subs(jets.D, 3), j)
        loss = sp.expand(-2*factor*sp.diff(angular-frozen, jets.D).subs(jets.D, 3))
        out["early_dimension_loss_formula_"+str(2*j)] = sp.expand(
            loss+2*factor*local.radial_polynomial(sp.diff(radial.variation.combined(j), jets.D), j).subs(jets.D, 3))
    return out


@cache
def controls():
    out = {}
    for order in (1, 2):
        raw = sp.expand(radial.mass**(4-2*order)*radial.actual(radial.laurent(order)["finite_MSbar_mu_m"]))
        top = 2*order
        defect = sp.factor(raw.coeff(jets.n[top-1])-order*W(raw.coeff(jets.n[top])))
        out["omitted_counterterm_finite_adjoint_defect_"+str(top)] = defect
        out["half_time_unit_mass_adjoint_fixture_"+str(top)] = defect.subs({radial.u: sp.Rational(1, 2), radial.mass: 1})
    return out
