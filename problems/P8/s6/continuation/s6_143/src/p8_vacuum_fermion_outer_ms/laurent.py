"""Exact leading gamma-function integral and its MS Laurent coefficients."""

from functools import cache

import sympy as sp


@cache
def data():
    e = sp.Symbol("epsilon", real=True)
    normal = (
        sp.exp(2 * sp.EulerGamma * e)
        * 4 ** (-e)
        * sp.sqrt(sp.pi)
        / 2
        * (sp.Rational(3, 2) - e)
        / (1 - e)
        * sp.gamma(1 + e)
        * sp.gamma(1 + 2 * e)
        / sp.gamma(sp.Rational(5, 2) + e)
    )
    L = sp.log(normal)
    first = sp.simplify(sp.diff(L, e).subs(e, 0))
    # Apply two exact shifts of psi_1(x+1)=psi_1(x)-1/x^2.
    trigamma_five_halves = sp.polygamma(1, sp.Rational(1, 2)) - 4 - sp.Rational(4, 9)
    second = sp.simplify(
        sp.diff(L, e, 2)
        .subs(e, 0)
        .xreplace({sp.polygamma(1, sp.Rational(5, 2)): trigamma_five_halves})
    )
    finite = sp.simplify((second + first**2) / 4)
    H = (
        sp.exp(2 * sp.EulerGamma * e)
        * 4 ** (-e)
        * sp.sqrt(sp.pi)
        / (2 * sp.gamma(sp.Rational(3, 2) - e))
    )
    prefactor_finite = sp.simplify(sp.diff(H, e).subs(e, 0) - sp.EulerGamma)
    return {
        "normalized_leading_gamma_function": normal,
        "dimensionless_leading_regulated_F": normal / (2 * e**2),
        "double_pole_coefficient": sp.Rational(1, 2),
        "simple_pole_coefficient": first / 2,
        "leading_finite_MS_coefficient": finite,
        "difference_prefactor_finite_coefficient": prefactor_finite,
        "scope": "The r=0 reference integral at mu=mF. The exact r=1/(4mF^2) correction is not omitted.",
        "checks": {
            "normalized_gamma_value": sp.simplify(normal.subs(e, 0) - 1),
            "first_logarithmic_derivative": first + sp.Rational(7, 3),
            "second_logarithmic_derivative": second - 5 - sp.pi**2 / 3,
            "leading_finite_coefficient": finite - sp.Rational(47, 18) - sp.pi**2 / 12,
            "difference_prefactor_finite_coefficient": prefactor_finite
            - (2 - 4 * sp.log(2)),
            "digamma_half_integer_recurrence": sp.polygamma(0, sp.Rational(5, 2))
            + sp.EulerGamma
            + 2 * sp.log(2)
            - sp.Rational(8, 3),
            "trigamma_two_shift_reference": trigamma_five_halves
            - sp.pi**2 / 2
            + sp.Rational(40, 9),
            "Gamma_function_shift_underlying_polygamma_recurrence": sp.gammasimp(
                sp.gamma(e + 1) - e * sp.gamma(e)
            ),
            "trigamma_half_reference": sp.polygamma(1, sp.Rational(1, 2))
            - sp.pi**2 / 2,
            "leading_pi_upper_bound": sp.Rational(47, 18)
            + sp.Rational(16, 12)
            - sp.Rational(71, 18),
            "leading_strict_upper_four_slack": 4
            - sp.Rational(71, 18)
            - sp.Rational(1, 18),
        },
    }
