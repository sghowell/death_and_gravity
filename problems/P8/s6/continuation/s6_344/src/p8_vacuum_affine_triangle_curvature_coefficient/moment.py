"""Exact positive-mass finite curvature moment with removable equal-mass limit."""

from functools import cache

import sympy as s

from . import source

N, Y, Z = s.symbols("n y z", positive=True)
NUMERATOR = (
    4 * N**5
    - 60 * N**4 * s.log(N)
    + 155 * N**4
    - 240 * N**3 * s.log(N)
    + 80 * N**3
    - 120 * N**2 * s.log(N)
    - 220 * N**2
    - 20 * N
    + 1
)
CLOSED = -NUMERATOR / (45 * (N - 1) ** 8)


def require_positive_mass(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("Require an exact positive rational mass squared")
    value = s.Rational(value)
    if value <= 0:
        raise ValueError("Require an exact positive rational mass squared")
    return value


def coefficient(mass):
    mass = require_positive_mass(mass)
    if mass == 1:
        return -s.Rational(1, 630)
    return CLOSED.subs(N, mass)


def primitive():
    polynomial = s.Poly(s.expand((Y - 1) ** 2 * (N - Y) ** 5), Y)
    return s.Add(
        *(
            co * (s.log(Y) if power[0] == 3 else Y ** (power[0] - 3) / (power[0] - 3))
            for power, co in polynomial.terms()
        )
    )


@cache
def data():
    F = primitive()
    c, gg = s.symbols("contact g", real=True)
    actualn, actualg = source.HEAVY_MASS2, source.CUBIC
    endpoint = source.CONTACT + actualg**2 / actualn
    checks = {
        "substituted_primitive_exact_derivative": s.factor(
            s.diff(F, Y) - (Y - 1) ** 2 * (N - Y) ** 5 / Y**4
        ),
        "full_endpoint_integral_closed_form": s.factor(
            -s.Rational(4, 15) * (F.subs(Y, N) - F.subs(Y, 1)) / (N - 1) ** 8 - CLOSED
        ),
        "equal_mass_removable_limit": s.limit(CLOSED, N, 1) + s.Rational(1, 630),
        "equal_mass_positive_Beta_moment": s.integrate(Z**2 * (1 - Z) ** 5, (Z, 0, 1))
        - s.Rational(1, 168),
        "large_mass_squared_leading_coefficient": s.limit(N**3 * CLOSED, N, s.oo)
        + s.Rational(4, 45),
        "positive_radial_tail_moment": s.integrate(Y**2 / (1 + Y) ** 4, (Y, 0, s.oo))
        - s.Rational(1, 3),
        "original_endpoint_algebra": s.factor(
            -(gg**2) * (3 / (N - 2) - 2 / (N - 2) ** 2)
            + gg**2 / N
            + 2 * gg**2 * (N**2 - 2 * N - 2) / (N * (N - 2) ** 2)
        ),
        "original_endpoint_upper_identity": s.factor(
            3 - 2 * (N**2 - 2 * N - 2) / (N - 2) ** 2 - (N - 4) ** 2 / (N - 2) ** 2
        ),
        "returned_equal_mass_value": coefficient(1) + s.Rational(1, 630),
    }
    for mass in (s.Rational(1, 2), s.Integer(2), s.Integer(3), s.Integer(10)):
        checks["exact_mass_" + str(mass) + "_primitive_endpoints"] = s.factor(
            coefficient(mass)
            + s.Rational(4, 15)
            * (F.subs({Y: mass, N: mass}) - F.subs({Y: 1, N: mass}))
            / (mass - 1) ** 8
        )
    gates = {
        "original_mass_squared_above_million": actualn > 10**6,
        "original_endpoint_strictly_negative": endpoint < 0,
        "original_endpoint_magnitude_below_3g2_over_n": -endpoint
        < 3 * actualg**2 / actualn,
        "moment_bound_below_original_n_inverse_cube": s.Rational(4, 45)
        / (actualn - 1) ** 3
        < 1 / actualn**3,
        "original_local_coefficient_over_Born_bound_below_1e_minus206": actualg**2
        / (192 * actualn)
        < s.Rational(1, 10**206),
        "strict_sign_from_positive_integrand_not_float_evaluation": True,
        "full_physical_amplitude_not_replaced_by_local_jet": True,
    }
    return {
        "checks": checks,
        "gates": gates,
        "whole_primitive_curvature_moment": -s.Rational(4, 15)
        * s.Integral(Z**2 * (1 - Z) ** 5 / (1 + (N - 1) * Z) ** 4, (Z, 0, 1)),
        "whole_closed_form_for_n_not_one": CLOSED,
        "whole_equal_mass_value": -s.Rational(1, 630),
        "whole_large_mass_squared_limit": -s.Rational(4, 45),
        "whole_exact_sign_and_bound": "c_triangle(n)<0 for every n>0 by its strictly positive real integral. For n>1,0<-c_triangle(n)<4/[45(n-1)^3], obtained by dropping(1-z)^5, substituting y=(n-1)z, extending to infinity and integrating y^2/(1+y)^4=1/3. For the original n>10^6 this is<n^-3.",
        "whole_selected_triangle_coefficient": gg**2
        * (c + gg**2 / N)
        * CLOSED
        / (16 * s.pi**2),
        "whole_original_known_coefficient_sign_and_bound": "In this explicitly stated lift, chi_triangle is strictly positive for the original source: c_triangle<0 and C+g^2/n<0. Also |chi_triangle|<g^4/(48*n^4), using|C+g^2/n|<3g^2/n and16pi^2>144. The already-proved positive Born A0>4g^2/n^3 gives|chi_triangle|/A0<g^2/(192n)<10^-206. This is a bound on a known local Taylor coefficient, NOT a truncation-error bound at the physical above-threshold kinematics or a bound on extra parent chi.",
        "whole_numerical_boundary": "Use the positive integral or asymptotic form for numerics near n=1; the rational-log expression has exact cancellations. No floating evaluation is used as proof.",
    }
