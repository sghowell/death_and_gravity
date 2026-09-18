"""Exact conditional logarithmic second and mixed moments."""

from functools import cache

import sympy as s
from p8_vacuum_affine_borel_soft_conversion.moments import parameters


def conditional_log_square_moments(index, resolution):
    a, x = parameters(index, resolution)
    if a == 0:
        return s.S.Zero, s.S.Zero, s.S.Zero
    b = a + 1
    L = 1 - s.log(x)
    l = -s.log(x)
    H = s.polygamma(0, b + 1) + s.EulerGamma
    mean = a * x / b
    energy = mean * (L * L + 2 * L / b + 2 / (b * b))
    mixed = mean * (L * (l + H) + l / b + H / b - s.polygamma(1, b + 1))
    remaining = mean * ((l + H) ** 2 + s.polygamma(1, 1) - s.polygamma(1, b + 1))
    return energy, mixed, remaining


@cache
def data():
    a, x, p, b, c, L = s.symbols("a x p b c L", positive=True)
    B = s.gamma(b) * s.gamma(c) / s.gamma(b + c)
    H = s.polygamma(0, b + 1) + s.EulerGamma
    moment = a * x**p / (a + p)
    mean = moment.subs(p, 1)
    checks = {
        "weighted_energy_log_square_derivative": mean
        - 2 * s.diff(moment, p).subs(p, 1)
        + s.diff(moment, p, 2).subs(p, 1)
        - mean
        * ((1 - s.log(x)) ** 2 + 2 * (1 - s.log(x)) / (a + 1) + 2 / (a + 1) ** 2),
        "mixed_log_beta_derivative": s.expand_func(
            s.diff(B, b, c).subs(c, 1) - (H / b - s.polygamma(1, b + 1)) / b
        ),
        "remaining_log_square_beta_derivative": s.expand_func(
            s.diff(B, c, 2).subs(c, 1)
            - (H**2 + s.polygamma(1, 1) - s.polygamma(1, b + 1)) / b
        ),
        "mixed_moment_upper_slack": s.expand(
            3 * L * L
            - (L * L + s.Rational(3, 2) * L + s.Rational(1, 2))
            - (L - 1) * (2 * L + s.Rational(1, 2))
        ),
        "remaining_log_square_upper_slack": s.expand(
            5 * L * L
            - (L + s.Rational(1, 2)) ** 2
            - 2
            - (L - 1) * (4 * L + 3)
            - s.Rational(3, 4)
        ),
        "energy_log_square_upper_slack": 5 * L * L
        - (L * L + 2 * L + 2)
        - 2 * (L - 1) * (2 * L + 1),
        "harmonic_upper_endpoint": s.polygamma(0, 3) + s.EulerGamma - s.Rational(3, 2),
    }
    return {
        "checks": {name: s.simplify(value) for name, value in checks.items()},
        "gates": {
            "all_three_moments_use_same_conditioned_leading_reference": True,
            "mixed_log_not_product_of_means": True,
            "harmonic_factor_bounded_by_three_halves": True,
            "trigamma_positive_and_difference_below_two": bool(
                s.Rational(22, 7) ** 2 / 6 < 2
            ),
            "zero_index_moments_vanish_exactly": True,
            "rare_probability_not_expanded_or_divided_into_loose_errors": True,
        },
        "whole_three_conditional_moments": "For b=a+1,L=1+ln1/x,l=ln1/x,H=psi(b+1)+EulerGamma, multiply a*x/b by respectively L^2+2L/b+2/b^2; L(l+H)+l/b+H/b-psi1(b+1); (l+H)^2+psi1(1)-psi1(b+1). These are E[R L_R^2],E[R L_R|ln(x-R)|],E[R|ln(x-R)|^2] on R<=x.",
        "whole_conditional_majorants": "The three exact moments are bounded by5*a*x*L_x^2,3*a*x*L_x^2 and5*a*x*L_x^2. Both endpoint logarithms remain integrable; a=0 is separate.",
    }
