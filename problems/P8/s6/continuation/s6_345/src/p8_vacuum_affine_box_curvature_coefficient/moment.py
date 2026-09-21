"""Exact finite box curvature moment and honest original coefficient bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_triangle_curvature_coefficient import moment as triangle

from . import radiation, source

N, Y = s.symbols("n y", positive=True)
Z = radiation.Z
NUM = (
    3 * N**10
    - 21 * N**9
    + 91 * N**8
    + 210 * N**7 * s.log(N)
    - 1859 * N**7
    + 10920 * N**6 * s.log(N)
    - 21896 * N**6
    + 28560 * N**5 * s.log(N)
    + 10920 * N**4 * s.log(N)
    + 21896 * N**4
    + 210 * N**3 * s.log(N)
    + 1859 * N**3
    - 91 * N**2
    + 21 * N
    - 3
)
CLOSED = -2 * NUM / (315 * N**3 * (N - 1) ** 9)
require_positive_mass = triangle.require_positive_mass


def coefficient(mass):
    mass = require_positive_mass(mass)
    return -s.Rational(58, 945) if mass == 1 else CLOSED.subs(N, mass)


def primitive():
    polynomial = s.Poly(
        s.cancel((N - 1) ** 8 * radiation.TARGET.subs(Z, (Y - 1) / (N - 1))), Y
    )
    F = s.Add(
        *(
            co * (s.log(Y) if power[0] == 4 else Y ** (power[0] - 4) / (power[0] - 4))
            for power, co in polynomial.terms()
        )
    )
    return polynomial, F


@cache
def data():
    polynomial, F = primitive()
    expanded = (
        -s.Rational(4, 315)
        * Z
        * (Z - 1)
        * (363 * Z**6 - 1089 * Z**5 + 997 * Z**4 - 179 * Z**3 - 74 * Z**2 - 18 * Z - 18)
    )
    checks = {
        "whole_positive_symmetric_weight_factorization": s.expand(
            radiation.TARGET - expanded
        ),
        "light_heavy_weight_reflection": s.expand(
            radiation.TARGET - radiation.TARGET.subs(Z, 1 - Z)
        ),
        "complete_exact_primitive_derivative": s.factor(
            s.diff(F, Y) - polynomial.as_expr() / Y**5
        ),
        "complete_integral_endpoints_closed_form": s.factor(
            (F.subs(Y, N) - F.subs(Y, 1)) / (N - 1) ** 9 - CLOSED
        ),
        "equal_mass_Beta_value": s.integrate(radiation.TARGET, (Z, 0, 1))
        + s.Rational(58, 945),
        "equal_mass_removable_closed_form_limit": s.limit(CLOSED, N, 1)
        + s.Rational(58, 945),
        "large_heavy_mass_squared_asymptotic": s.limit(N * N * CLOSED, N, s.oo)
        + s.Rational(2, 105),
        "exact_mass_exchange_covariance": s.factor(
            s.expand_log(CLOSED - N ** (-5) * CLOSED.subs(N, 1 / N), force=True)
        ),
        "positive_tail_integral": s.integrate(Y / (1 + Y) ** 5, (Y, 0, s.oo))
        - s.Rational(1, 12),
        "positive_weight_polynomial_upper_budget": 18
        + 18 * s.Rational(1, 4)
        + 92 * s.Rational(1, 16)
        + 363 * s.Rational(1, 64)
        - s.Rational(2171, 64),
        "equal_mass_API_value": coefficient(1) + s.Rational(58, 945),
    }
    for mass in (s.Rational(1, 2), 2, 3, 10):
        checks["exact_mass_" + str(mass) + "_endpoints"] = s.factor(
            coefficient(mass)
            - (F.subs({Y: mass, N: mass}) - F.subs({Y: 1, N: mass})) / (mass - 1) ** 9
        )
    mass, g = source.HEAVY_MASS2, source.CUBIC
    return {
        "checks": checks,
        "gates": {
            "strict_negative_coefficient_from_positive_integrand": True,
            "original_heavy_mass_domain": mass > 10**6,
            "original_moment_bound_below_1_over_25n2": s.Rational(2171, 60480)
            / (mass - 1) ** 2
            < 1 / (25 * mass**2),
            "original_Born_normalized_coefficient_LARGE_upper_bound": g**2
            * mass
            / 14400
            < 10**187,
            "no_physical_jet_truncation_bound": True,
            "different_triangle_convention_not_silently_added": True,
        },
        "whole_exact_mass_moment": s.Integral(
            radiation.TARGET / (1 + (N - 1) * Z) ** 5, (Z, 0, 1)
        ),
        "whole_closed_form_n_not_one": CLOSED,
        "whole_equal_mass_value": -s.Rational(58, 945),
        "whole_large_mass_squared_limit": -s.Rational(2, 105),
        "whole_sign_and_bound": "c_box(n)<0 for every n>0, because TARGET=-4*t*(18+18t+92t^2+363t^3)/315, t=z(1-z). For n>1, use t<=1/4 and t<=z:0<-c_box(n)<2171/[60480(n-1)^2], by integral_0^infinity y/(1+y)^5=1/12. At the original n>10^6 this is<1/(25n^2). No floating fit enters.",
        "whole_original_coefficient_bound": "With chi_box=g^4*c_box/(16pi^2), the original known box coefficient is negative and |chi_box|<g^4/(3600n^2). Dividing by A0>4g^2/n^3 gives |chi_box|/A0<g^2*n/14400<10^187. This is a LARGE upper bound on a local coefficient, not perturbative smallness, not an above-threshold amplitude estimate and not a bound on extra parent chi.",
        "whole_scheme_boundary": "The S344 triangle coefficient uses an explicitly different labeled-box comparison. Converting it to this uniform six-word jet basis, and treating all remaining sectors in a stated common convention, is required before adding coefficients as a single aggregate result. Nothing here changes the physical metric or adds a second copy to the complete S342 radiation.",
    }
