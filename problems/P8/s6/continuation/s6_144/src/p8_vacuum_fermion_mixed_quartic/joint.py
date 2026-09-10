"""An absolutely integrable joint bound for the full nonzero-soft box remainder."""

from functools import cache

import sympy as sp

SUNSET_GAMMA = sp.Rational(7, 4)


def sunset_constant(alpha, beta, gamma=SUNSET_GAMMA):
    return (
        sp.gamma(2 - gamma)
        * sp.gamma(alpha + beta + gamma - 4)
        * sp.gamma(alpha + gamma - 2)
        * sp.gamma(beta + gamma - 2)
        / (sp.gamma(alpha) * sp.gamma(beta) * sp.gamma(alpha + beta + 2 * gamma - 4))
    )


@cache
def data():
    ab = tuple(
        (sp.Rational(a + 1, 2), sp.Rational(b, 2)) for a, b in ((1, 3), (2, 2), (3, 1))
    ) + tuple(
        (sp.Rational(a, 2), sp.Rational(b + 1, 2)) for a, b in ((1, 3), (2, 2), (3, 1))
    )
    g = sp.Rational(7, 4)
    checks = {
        "fractional_majorant_high_exponent_gap": 2
        - sp.Rational(7, 4)
        - sp.Rational(1, 4),
        "Euler_gamma_upper_reference": 4 + sp.factorial(3) - 10,
        "three_denominator_gamma_lower_bound": sp.Rational(1, 2) ** 3
        - sp.Rational(1, 8),
        "individual_sunset_constant_upper": 10**4 * 8 - 80000,
        "two_minimum_radius_terms": 2 * 80000 - 160000,
        "one_channel_assignment_remainder_prefactor": 6 * 128 * 360 * 16 * 160000
        - 707788800000,
        "all_channel_remainder_prefactor": 3 * 707788800000 - 2123366400000,
        "mass_suppression": 4 - sp.Rational(5, 2) - g + sp.Rational(1, 4),
    }
    for i, (a, b) in enumerate(ab):
        checks[f"sunset_sum_{i}"] = a + b - sp.Rational(5, 2)
        checks[f"sunset_beta_denominator_{i}"] = a + b + 2 * g - 4 - 2
    return {
        "Schwinger_integral": "Integral d4k d4l /Q-normalized angular measures of (m^2+k^2)^(-alpha)(m^2+l^2)^(-beta)((k-l)^2)^(-gamma) = (m^2)^(4-alpha-beta-gamma) sunset_constant/Q^2.",
        "exponents": ab,
        "gamma": g,
        "individual_exact_constants": tuple(
            sp.simplify(sunset_constant(a, b)) for a, b in ab
        ),
        "individual_constant_upper": 80000,
        "all_channel_remainder_prefactor": 2123366400000,
        "bound": "E_delta <=2123366400000 N Y^2 V/(Q^2 sqrt(m)), V=L+g/(M-3)+4g/M.",
        "method": "Cauchy in fermion soft momenta only. Keep the light lines and the whole W vertex at the physical external momenta. No large-radius contour is assigned to a light or heavy line.",
        "checks": checks,
    }
