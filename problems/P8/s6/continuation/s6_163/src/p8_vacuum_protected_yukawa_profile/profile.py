"""A globally real-analytic bounded Yukawa argument, not an entire function."""

from fractions import Fraction
from functools import cache

import sympy as s

Z = s.Symbol("z", real=True)
R, M, Y = s.symbols("R m Y", positive=True)
ETA = s.Rational(1, 100)


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, s.Rational)):
        raise TypeError("Require a finite exact rational")
    return s.Rational(value)


def coefficient(k):
    if type(k) is not int or k < 0:
        raise TypeError("Require a nonnegative native binomial index")
    return s.binomial(-s.Rational(1, 8), k)


def expression(z=Z, radius=R):
    return z / (1 + (z / radius) ** 8) ** s.Rational(1, 8)


def gap(mass, yukawa_squared_upper, radius, eta=ETA):
    m, y2, r, e = map(exact, (mass, yukawa_squared_upper, radius, eta))
    if min(m, y2, r) <= 0 or not 0 < e < 1:
        raise ValueError(
            "Require positive scales and a relative mass allowance in (0,1)"
        )
    ratio = y2 * r * r / (m * m)
    if not ratio < e * e:
        raise ValueError("The strict paired-mass gap condition fails")
    return {
        "relative_mass_shift_squared_upper": ratio,
        "both_mass_strict_lower": (1 - e) * m,
        "both_mass_strict_upper": (1 + e) * m,
    }


@cache
def data():
    f = expression()
    x = Z / R
    d1 = (1 + x**8) ** (-s.Rational(9, 8))
    d2 = -9 * x**7 / (R * (1 + x**8) ** s.Rational(17, 8))
    germ = sum(coefficient(k) * Z ** (8 * k + 1) / R ** (8 * k) for k in range(4))
    q2 = s.Symbol("q_squared", nonnegative=True)
    first = Z - Z**9 / (8 * R**8)
    pair = (q2 + M * M + Y * first**2) ** 2 - 4 * M * M * Y * first**2
    oldpair = (q2 + M * M + Y * Z**2) ** 2 - 4 * M * M * Y * Z**2
    log10 = s.expand(pair - oldpair).coeff(Z, 10) / (q2 + M * M) ** 2
    return {
        "candidate": "GY14-SAT8: active masses m +/- y f_R(Phi), R=10^300; two active and twelve inert SU(3) fundamental Dirac flavors. All scalar/gauge terms and the lower-field reference prescription are unchanged.",
        "bounded_argument": f,
        "argument_first_derivative": d1,
        "argument_second_derivative": d2,
        "vacuum_germ_through_degree_twenty_five": germ,
        "first_new_scalar_degree": 9,
        "paired_log_determinant_first_difference_coefficient": log10,
        "real_domain": "All real z, with R>0 and the unique positive real eighth root. |f_R(z)|<R, f'_R in (0,1], and |f''_R|<=9/R. This pointwise mass profile has no new field or kinetic derivative.",
        "complex_domain": "The vacuum Taylor branch is analytic for |z|<R. The eight branch points z/R=exp(i(2j+1)pi/8) lie on |z|=R. The function is real analytic on all of the real axis, not entire in complex z.",
        "checks": {
            "odd_real_profile": s.simplify(expression(-Z) + f),
            "derivative_first": s.simplify(s.diff(f, Z) - d1),
            "derivative_second": s.simplify(s.diff(f, Z, 2) - d2),
            "zeroth_binomial": coefficient(0) - 1,
            "first_new_binomial": coefficient(1) + s.Rational(1, 8),
            "second_new_binomial": coefficient(2) - s.Rational(9, 128),
            "third_new_binomial": coefficient(3) + s.Rational(51, 1024),
            "no_scalar_degree_two_through_eight_change": sum(
                s.expand(germ - Z).coeff(Z, j) for j in range(2, 9)
            ),
            "paired_logdet_detects_higher_field_change": s.factor(
                log10 + Y * (q2 - M * M) / (2 * R**8 * (q2 + M * M) ** 2)
            ),
            "nonzero_zero_momentum_higher_field_control": s.factor(
                log10.subs(q2, 0) - Y / (2 * R**8 * M * M)
            ),
            "unprotected_positive_crossing": M - s.sqrt(Y) * (M / s.sqrt(Y)),
            "saturation_eighth_power_identity": s.simplify(
                f**8 / R**8 - x**8 / (1 + x**8)
            ),
        },
    }
