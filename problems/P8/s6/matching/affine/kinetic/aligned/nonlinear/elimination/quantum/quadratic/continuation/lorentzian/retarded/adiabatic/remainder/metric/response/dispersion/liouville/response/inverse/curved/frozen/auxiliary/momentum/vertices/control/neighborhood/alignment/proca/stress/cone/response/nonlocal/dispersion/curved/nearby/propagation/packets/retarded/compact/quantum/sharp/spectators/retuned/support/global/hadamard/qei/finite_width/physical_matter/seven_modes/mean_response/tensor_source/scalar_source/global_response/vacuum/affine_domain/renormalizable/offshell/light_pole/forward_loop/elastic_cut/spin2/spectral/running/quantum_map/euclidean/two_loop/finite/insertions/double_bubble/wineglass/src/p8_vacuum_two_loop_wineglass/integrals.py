"""Anchored all-radius integrals for the subtracted wineglass majorants."""

from functools import cache

import sympy as sp
from p8_vacuum_two_loop_double_bubble import triangle
from p8_vacuum_two_loop_insertions import radial


@cache
def data():
    y, t, M = sp.symbols(
        "positive_radial positive_sqrt_radial positive_heavy_mass_squared",
        positive=True,
    )
    delta = sp.Rational(1, 4)
    primitive_half = (
        sp.Rational(3, 2) * sp.atan(2 * t) - t * (20 * t**2 + 3) / (4 * t**2 + 1) ** 2
    )
    primitive_integer = -1 / (y + delta) + delta / (2 * (y + delta) ** 2)
    beta_half = 3 * sp.pi / 4
    beta_integer = sp.Integer(2)
    first = 8 * beta_half + 4 * beta_integer
    second = sp.Rational(33, 16) * beta_half
    ell = sp.Symbol("positive_heavy_mass_logarithm_upper", positive=True)
    oldtriangle, oldradial = triangle.data(), radial.data()
    J = oldtriangle["gap_one_quarter_radial_integral"].subs(oldtriangle["M"], M)
    K = oldradial["actual_gap_one_quarter_integral"].subs(oldradial["M"], M)
    loggap = t - sp.log(1 + t**2)
    checks = {
        "half_power_radial_substitution_primitive": sp.factor(
            sp.diff(primitive_half, t) - 2 * t**4 / (t**2 + delta) ** 3
        ),
        "half_power_radial_zero_anchor": primitive_half.subs(t, 0),
        "half_power_radial_infinity_anchor": sp.limit(primitive_half, t, sp.oo)
        - beta_half,
        "integer_radial_primitive_derivative": sp.factor(
            sp.diff(primitive_integer, y) - y / (y + delta) ** 3
        ),
        "integer_radial_zero_anchor": -primitive_integer.subs(y, 0) - beta_integer,
        "integer_radial_infinity_anchor": sp.limit(primitive_integer, y, sp.oo),
        "logarithm_sqrt_majorant_positive_derivative": sp.factor(
            sp.diff(loggap, t) - (t - 1) ** 2 / (1 + t**2)
        ),
        "logarithm_sqrt_majorant_zero_anchor": loggap.subs(t, 0),
        "light_inverse_square_difference_numerator_bound": sp.expand(
            sp.Rational(11, 4) * (y + 1)
            - (2 * y + sp.Rational(11, 4))
            - sp.Rational(3, 4) * y
        ),
        "light_inverse_square_difference_prefactor": sp.Rational(3, 4)
        * sp.Rational(11, 4)
        - sp.Rational(33, 16),
        "first_subtracted_integral_bound": first - (6 * sp.pi + 8),
        "second_subtracted_integral_bound": second - sp.Rational(99, 64) * sp.pi,
        "subtracted_radial_bound_below_forty_margin": 40
        - 32
        - sp.Rational(99, 16)
        - sp.Rational(29, 16),
        "heavy_log_product_argument_gap": sp.expand(
            (1 + M) * (1 + y / M) - (1 + y) - M - y / M
        ),
        "same_complete_triangle_integral": sp.simplify(
            J - M * sp.log(4 * M) / (M - delta) ** 2 + 1 / (M - delta)
        ),
        "same_complete_decaying_log_integral": sp.simplify(
            K - sp.log(4 * M) / (M - delta)
        ),
        "all_radius_heavy_logarithm_bound": sp.expand(
            (ell + 3) * 2 * ell / M + 2 * ell / M - 2 * ell * (ell + 4) / M
        ),
        "heavy_numerator_twenty_prefactor": 4 + 3 * sp.Rational(16, 3) - 20,
        "complete_channel_subtracted_and_heavy_prefactors": 3 * 2 * 40 - 240,
        "complete_channel_heavy_coupling_prefactor": 3 * 20 * 2 * sp.Rational(1, 3)
        - 40,
    }
    return {
        "half_power_radial_primitive_after_sqrt_substitution": primitive_half,
        "integer_radial_primitive": primitive_integer,
        "half_power_radial_integral": beta_half,
        "integer_radial_integral": beta_integer,
        "subtracted_shift_part_radial_upper": first,
        "subtracted_light_gap_part_radial_upper": second,
        "complete_subtracted_wineglass_radial_upper": sp.Integer(40),
        "same_full_triangle_radial_integral": J,
        "same_full_heavy_logarithm_radial_integral": K,
        "complete_heavy_weighted_logarithm_upper": 2 * ell * (ell + 4) / M,
        "complete_family_quartic_cubed_prefactor_upper": (240 + 40 * ell * (ell + 4))
        / (16 * sp.pi**2) ** 2,
        "scope": "Every radial integral extends to infinity. The subtracted local core uses the decaying logarithm difference and an inverse-square denominator difference. Heavy terms use full unexpanded propagators and inherited exact integrals. Cauchy at unit radius applies to the complete integrated family.",
        "checks": checks,
    }
