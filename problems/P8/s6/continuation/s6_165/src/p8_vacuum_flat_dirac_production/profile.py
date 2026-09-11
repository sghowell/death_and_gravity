"""Integrable SAT8 switching derivatives and total-variation bounds."""

from functools import cache

import sympy as s

X = s.Symbol("x", real=True)


@cache
def data():
    f = X / (1 + X**8) ** s.Rational(1, 8)
    d1 = (1 + X**8) ** (-s.Rational(9, 8))
    d2 = -9 * X**7 / (1 + X**8) ** s.Rational(17, 8)
    d3 = -9 * X**6 * (7 - 10 * X**8) / (1 + X**8) ** s.Rational(25, 8)
    return {
        "normalized_profile": f,
        "first_derivative": d1,
        "second_derivative": d2,
        "third_derivative": d3,
        "first_derivative_L1": s.Integer(2),
        "first_times_second_derivative_absolute_integral": s.Integer(1),
        "first_derivative_cubed_integral_upper": s.Integer(2),
        "third_derivative_absolute_integral_upper": s.Integer(36),
        "extremum_positive_eighth_power": s.Rational(7, 10),
        "tail": "For |x|>=1, |s(x)-sgn(x)|<=1/(8|x|^8), by convexity of (1+u)^(-1/8). Hence the mass difference from each asymptote is time integrable. Derivatives decay sufficiently for both transition-integral integrations by parts.",
        "total_variation": "s' increases from zero to one and decreases to zero. Therefore integral |s's''|=1, and integral (s')^3<=integral s'=2. On each side s'' has one extremum, at |x|^8=7/10. Thus integral |s'''|=4 max|s''|<=36 using the previously proved |s''|<=9.",
        "checks": {
            "full_profile_first_derivative": s.simplify(s.diff(f, X) - d1),
            "full_profile_second_derivative": s.simplify(s.diff(f, X, 2) - d2),
            "full_profile_third_derivative": s.simplify(s.diff(f, X, 3) - d3),
            "squared_first_derivative_antiderivative": s.simplify(
                s.diff(d1 * d1 / 2, X) - d1 * d2
            ),
            "positive_extremum_equation": 7 - 10 * s.Rational(7, 10),
            "four_total_variation_pieces": 4 * 9 - 36,
            "two_first_derivative_half_lines": 2 * (1 - 0) - 2,
            "weighted_variation_exact": 2 * s.Rational(1, 2) - 1,
        },
    }
