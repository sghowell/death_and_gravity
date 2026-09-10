"""Cancellation-free exact moment series for one-loop on-shell anchors."""

from fractions import Fraction
from functools import cache

import sympy as sp


def moment(n):
    if type(n) is not int:
        raise TypeError("Require a native integer moment index")
    if n < 0 or n > 64:
        raise ValueError("Require a moment index between zero and 64")
    return sp.factorial(n) ** 2 / sp.factorial(2 * n + 1)


def enclosure(mass, terms=4):
    if isinstance(mass, bool) or not isinstance(mass, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational fermion mass")
    m = sp.Rational(mass)
    if m < 36:
        raise ValueError("Require mF>=36")
    if type(terms) is not int:
        raise TypeError("Require a native integer anchor truncation")
    if not 1 <= terms <= 16:
        raise ValueError("Require one to sixteen retained anchor terms")
    coefficients = {n: 3 * moment(n) / (n * (2 * n + 3)) for n in range(1, terms + 1)}
    S = sum(coefficients[n] / m ** (2 * n) for n in coefficients)
    P = sum((n + 1) * coefficients[n] / m ** (2 * n) for n in coefficients)
    R = sum(n * coefficients[n] / m ** (2 * n) for n in coefficients)
    ratio = 1 / (4 * m * m)
    tail = sp.Rational(3, 5) * ratio ** (terms + 1) / (1 - ratio)
    return {
        "mass": m,
        "retained_terms": terms,
        "coefficients": coefficients,
        "mass_increment_over_2NY_div_Q_lower": sp.Rational(2, 3) - S - tail,
        "mass_increment_over_2NY_div_Q_upper": sp.Rational(2, 3) - S,
        "slope_over_2NY_div_Q_lower": sp.Rational(2, 3) - P - 2 * tail,
        "slope_over_2NY_div_Q_upper": sp.Rational(2, 3) - P,
        "OS_zero_remainder_over_2NY_div_Q_lower": -R - tail,
        "OS_zero_remainder_over_2NY_div_Q_upper": -R,
        "common_positive_series_tail_upper": tail,
    }


@cache
def data():
    n = sp.symbols("positive_moment_index", integer=True, positive=True)
    b = sp.symbols("positive_moment_b", positive=True)
    ratio = (n + 1) / (2 * (2 * n + 3))
    c = 3 * b / (n * (2 * n + 3))
    x = sp.symbols("unit_interval_x", real=True)
    checks = {
        "moment_successive_ratio": sp.simplify(
            sp.factorial(n + 1) ** 2
            / sp.factorial(2 * n + 3)
            / (sp.factorial(n) ** 2 / sp.factorial(2 * n + 1))
            - ratio
        ),
        "mass_anchor_grouped_coefficient": sp.factor(
            4 * ratio * b / (n + 1) - b / n + c
        ),
        "slope_anchor_grouped_coefficient": sp.factor(
            4 * ratio * b - (1 + 1 / n) * b + (n + 1) * c
        ),
        "OS_zero_anchor_grouped_coefficient": sp.factor(
            4 * ratio * b * n / (n + 1) - b + n * c
        ),
        "mass_series_weight_bound_residual": sp.factor(
            sp.Rational(3, 5)
            - 3 / (n * (2 * n + 3))
            - 3 * (n - 1) * (2 * n + 5) / (5 * n * (2 * n + 3))
        ),
        "slope_series_weight_bound_residual": sp.factor(
            sp.Rational(6, 5)
            - 3 * (n + 1) / (n * (2 * n + 3))
            - 3 * (n - 1) * (4 * n + 5) / (5 * n * (2 * n + 3))
        ),
        "OS_zero_series_weight_bound_residual": sp.factor(
            sp.Rational(3, 5) - 3 / (2 * n + 3) - 6 * (n - 1) / (5 * (2 * n + 3))
        ),
    }
    for degree in range(9):
        checks[f"independent_beta_integral_moment_{degree}"] = sp.integrate(
            (x * (1 - x)) ** degree, (x, 0, 1)
        ) - moment(degree)
    return {
        "moments": "b_n=integral_0^1 [x(1-x)]^n dx=(n!)^2/(2n+1)!",
        "positive_grouped_coefficients": "a_n=3 b_n/[n(2n+3)]",
        "mass_anchor_series": "f(1)=f(0)+(2NY/Q)[2/3-sum_{n>=1} a_n/mF^(2n)]",
        "slope_anchor_series": "f'(1)=(2NY/Q)[2/3-sum_{n>=1}(n+1)a_n/mF^(2n)]",
        "OS_zero_remainder_series": "f_R(0)=-(2NY/Q)sum_{n>=1}n a_n/mF^(2n)",
        "tail_proof": "b_n<=4^-n, a_n/b_n<=3/5, (n+1)a_n/b_n<=6/5, n a_n/b_n<=3/5. Sum the omitted geometric tail; bounds cover every remaining term.",
        "scope": "No subtraction of two huge decimal masses is used to recover a tiny on-shell remainder. The anchor interfaces return normalized rational enclosures; physical coupling and loop factors are retained separately.",
        "checks": checks,
    }
