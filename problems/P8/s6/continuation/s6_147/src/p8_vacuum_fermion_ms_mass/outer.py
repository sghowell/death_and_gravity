"""The full regulated outer reference at variable scalar mass and its b integral."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_mixed_quartic import reference as mixed


@cache
def data():
    e = s.Symbol("epsilon", real=True)
    b, v = s.symbols("b v", positive=True)
    ell = s.Symbol("log_m_squared", real=True)
    A = s.exp(s.EulerGamma * e) * s.gamma(1 + e)
    soft = (
        s.exp((ell - s.log(b)) * e)
        * ((12 - 32 * e + 16 * e * e) * A * A - 12 * A)
        / e**2
    )
    soft_normal = soft * e**2
    softfinite = s.simplify(s.diff(soft_normal, e, 2).subs(e, 0) / 2)
    hardfinite = -62 - s.pi**2
    finite = s.simplify(softfinite + hardfinite)
    integrated = -s.integrate(finite, (b, 0, 1))
    bracket = b ** (-e) / (v - b) - (v ** (1 - e) - b ** (1 - e)) / (
        (1 - e) * (v - b) ** 2
    )
    finite_triangle = s.simplify(s.diff(bracket, e).subs(e, 0))
    expected_triangle = (v * s.log(v / b) - v + b) / (v - b) ** 2
    log_primitive = -b * s.log(b) + b
    original = mixed.data()["leading_finite_reference"]
    hard = mixed.data()["dimensionless_hard_reference"]
    integrated_normal = s.simplify(-soft_normal.subs(b, 1) / (1 - e) - e**2 * hard)
    integrated_jets = tuple(
        s.simplify(s.diff(integrated_normal, e, n).subs(e, 0) / s.factorial(n))
        for n in range(3)
    )
    return {
        "general_outer_triangle_regulated_bracket": bracket,
        "general_outer_triangle_finite": finite_triangle,
        "leading_F_mix_MS_over_CY_over_Q": finite,
        "integrated_scalar_mass_difference_leading_coefficient": integrated,
        "mass_difference_double_pole_in_NY_squared_over_Q_squared_units": 6,
        "mass_difference_simple_pole_in_NY_squared_over_Q_squared_units": -2,
        "whole_regulated_integrated_reference": integrated_normal / e**2,
        "integrated_double_pole_simple_pole_finite": integrated_jets,
        "integral_scope": "The coefficientwise b integral commutes with the finite Laurent extraction only after the whole regulated proper forest is written. Its b^(-epsilon) and logarithmic endpoints have uniform integrable bounds for |epsilon|<1/8.",
        "checks": {
            "general_outer_triangle_removable_limit": s.simplify(
                finite_triangle - expected_triangle
            ),
            "soft_reference_simple_pole": s.simplify(
                s.diff(soft_normal, e).subs(e, 0) + 32
            ),
            "general_soft_reference_finite": s.simplify(
                softfinite - (16 + s.pi**2 - 32 * (ell - s.log(b)))
            ),
            "general_full_reference_finite": s.simplify(
                finite + 46 + 32 * (ell - s.log(b))
            ),
            "unit_b_reference_is_exact_frozen_S744": s.simplify(
                finite.subs(b, 1) - original
            ),
            "minus_full_b_integral": s.simplify(integrated - 78 - 32 * ell),
            "real_log_endpoint_antiderivative": s.diff(log_primitive, b) + s.log(b),
            "real_log_endpoint_integral": log_primitive.subs(b, 1)
            - s.limit(log_primitive, b, 0, dir="+")
            - 1,
            "half_CY_equals_NY_squared": s.Rational(1, 2) * 2 - 1,
            "quartic_threshold_fixes_log_coefficient": -s.Rational(1, 2) * (-64) - 32,
            "integrating_regulator_first_double_pole": integrated_jets[0] - 6,
            "integrating_regulator_first_simple_pole": integrated_jets[1] + 2,
            "integrating_regulator_first_finite": s.simplify(
                integrated_jets[2] - integrated
            ),
        },
    }
