"""Continuous C4-to-C0 bound for the displayed finite local component only."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact

from . import compact, jets, radial

LEFT, RIGHT = -sp.Rational(1, 2), sp.Rational(1, 2)


def rational_bound(value):
    if not isinstance(value, sp.Expr) or value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("Require an exact finite rational function")
    if value.free_symbols-{radial.u}:
        raise ValueError("Only the fixed clock time may remain in a normalized coefficient")
    numerator, denominator = sp.fraction(sp.factor(value))
    try:
        p, d = sp.Poly(numerator, radial.u), sp.Poly(denominator, radial.u)
    except sp.PolynomialError as error:
        raise ValueError("Require a rational clock coefficient") from error
    if d.degree() % 2:
        raise ValueError("Require a positive constant times a power of 1+u^2")
    power, constant = int(d.degree())//2, d.eval(0)
    if constant <= 0 or sp.expand(denominator-constant*(1+radial.u**2)**power) != 0:
        raise ValueError("Unsupported denominator for this continuous bound")
    return sum(abs(coefficient)*RIGHT**exponent[0] for exponent, coefficient in p.terms())/constant


@cache
def coefficient_envelopes():
    out = {}
    for order in (0, 1, 2):
        operator = sp.expand(radial.matched(order)/radial.mass**(4-2*order))
        out[order] = {j: rational_bound(sp.factor(operator.coeff(jets.n[j]))) for j in range(5)}
    return out


def operator_bound(planck_time_product, mass_time_product):
    L, R = exact(planck_time_product, "M_tau"), exact(mass_time_product, "m0_tau")
    if not bool(L > 0) or not bool(R >= 1000):
        raise ValueError("Require M*tau>0 and m0*tau>=1000")
    # pi^2>9; all loop components are normalized by 64*pi^2*L^2.
    pieces = {j: sp.factor(R**(4-2*j)*sum(row.values())/(576*L**2))
              for j, row in coefficient_envelopes().items()}
    return {"by_adiabatic_half_order": pieces, "C4_to_C0_upper_bound": sum(pieces.values()),
            "source_norm": "max over time derivatives zero through four on [-1/2,1/2]",
            "scope": "Homogeneous finite local subtraction component only; not the nonlocal response or a feedback contraction"}


@cache
def closed_coefficients():
    u, m, r = radial.u, radial.mass, 1+radial.u**2
    return {0: {"second_time_square": sp.Integer(0), "first_time_square": sp.Integer(0),
                "source_square": -sp.Rational(1652, 6561)*m**4/r**6},
            1: {"second_time_square": sp.Integer(0), "first_time_square": -sp.Rational(800, 6561)*m**2/r**6,
                "source_square": -sp.Rational(128, 6561)*m**2*(917*u**2+164)/r**8},
            2: {"second_time_square": sp.Rational(628, 6561)/r**6,
                "first_time_square": sp.Rational(80, 19683)*(1483*u**2+459)/r**8,
                "source_square": sp.Rational(64, 6561)*(3651*u**4+4285*u**2+338)/r**10}}


@cache
def checks():
    return {"closed_finite_"+str(order)+"_"+name: sp.factor(compact.coefficients(order)[name]-target)
            for order, row in closed_coefficients().items() for name, target in row.items()}
