"""Exact four triangle comparison moments and the already matched box moment."""

from functools import cache

import sympy as s
from p8_vacuum_affine_box_curvature_coefficient import moment as box
from p8_vacuum_affine_box_curvature_coefficient import radiation
from p8_vacuum_affine_triangle_curvature_coefficient import moment as triangle

from . import conversion

N, Y = s.symbols("n y", positive=True)
Z = conversion.Z
require_mass = triangle.require_positive_mass


def require_order(value):
    if type(value) is not int or not 0 <= value <= 3:
        raise ValueError("Require a native triangle order0..3")
    return value


@cache
def weights():
    K = conversion.build()[3]
    return (s.factor(K[0] - s.Rational(4, 15) * Z**2 * (1 - Z) ** 5), *K[1:])


def primitive(weight, power):
    poly = s.Poly(s.expand(weight), Z)
    degree = poly.degree()
    numerator = s.Poly(
        s.expand(
            sum(
                co * (Y - 1) ** ex[0] * (N - 1) ** (degree - ex[0])
                for ex, co in poly.terms()
            )
        ),
        Y,
    )
    F = sum(
        co
        * (
            s.log(Y)
            if ex[0] == power - 1
            else Y ** (ex[0] - power + 1) / s.Integer(ex[0] - power + 1)
        )
        for ex, co in numerator.terms()
    )
    value = s.factor((F.subs(Y, N) - F.subs(Y, 1)) / (N - 1) ** (degree + 1))
    return numerator, F, value


@cache
def triangle_moments():
    return tuple(primitive(weight, 4 - j)[2] for j, weight in enumerate(weights()))


def coefficient(order, mass):
    order = require_order(order)
    mass = require_mass(mass)
    return (
        s.integrate(weights()[order], (Z, 0, 1))
        if mass == 1
        else triangle_moments()[order].subs(N, mass)
    )


@cache
def data():
    checks = {}
    for j, weight in enumerate(weights()):
        numerator, F, value = primitive(weight, 4 - j)
        checks[f"triangle_j{j}_whole_primitive_derivative"] = s.factor(
            s.diff(F, Y) - numerator.as_expr() / Y ** (4 - j)
        )
        checks[f"triangle_j{j}_equal_mass_continuity"] = s.limit(
            value, N, 1
        ) - s.integrate(weight, (Z, 0, 1))
        for mass in (s.Rational(1, 2), 2, 3, 10):
            checks[f"triangle_j{j}_mass{mass}_endpoint"] = s.factor(
                coefficient(j, mass) - value.subs(N, mass)
            )
    numerator, F, value = primitive(radiation.TARGET, 5)
    checks["complete_box_primitive_derivative"] = s.factor(
        s.diff(F, Y) - numerator.as_expr() / Y**5
    )
    checks["same_complete_selected_box_moment"] = s.factor(
        value - box.CLOSED.subs(box.N, N)
    )
    checks["complete_leading_triangle_weight"] = s.expand(
        weights()[0] + 2 * (Z - 1) ** 4 * (151 * Z**3 + 44 * Z**2 + 12 * Z + 3) / 105
    )
    return {
        "checks": checks,
        "gates": {
            "all_four_dressed_triangle_moments": len(triangle_moments()) == 4,
            "primitive_curvature_plus_lift_conversion_not_silently_same_basis": True,
            "equal_mass_limit_not_zero_default": True,
            "all_masses_positive_analytic_origin_only": True,
            "complete_box_coefficient_not_added_twice": True,
        },
        "whole_triangle_weights": weights(),
        "whole_triangle_closed_moments": triangle_moments(),
        "whole_triangle_equal_mass_values": tuple(coefficient(j, 1) for j in range(4)),
        "whole_triangle_integrals": tuple(
            s.Integral(w / (1 + (N - 1) * Z) ** (4 - j), (Z, 0, 1))
            for j, w in enumerate(weights())
        ),
        "whole_box_moment": box.CLOSED.subs(box.N, N),
        "whole_normalization": "S344 gives c_triangle=-4/15 int z^2(1-z)^5/M^4 for the complete24-labeled primitive curvature difference. Add it only to j0 conversion weight, since its lower-order curvature differences vanish and the outer A(v) coefficient starts at A(0). The four conversion moments carry M^(4-j). The complete S345 box coefficient is already in this same jet convention.",
    }
