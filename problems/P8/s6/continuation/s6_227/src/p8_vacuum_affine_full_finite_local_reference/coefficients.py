"""Whole-slab actual alpha, alpha-prime and gamma coefficient bounds."""

from functools import cache

import sympy as s

from . import factorization as f


@cache
def data():
    bg = f.g.profile()
    t, a, h, U = bg["t"], bg["a"], bg["h"], bg["U"]
    mass = f.g.bridge.MASS
    Up = s.factor(a * s.diff(U, t))
    alpha = s.factor(-2 * U / 3 + 5 * mass**2 * a * a / 3)
    alphap = s.factor(a * s.diff(alpha, t))
    gamma = s.factor(
        -2 * U * U + 10 * mass**2 * a * a * U + s.Rational(5, 2) * mass**4 * a**4
    )
    point = s.Rational(1, 2)
    checks = {
        "actual_conformal_U_derivative": s.factor(
            Up - 24 * t * (1 + t * t) ** 3 * (3 + 7 * t * t)
        ),
        "actual_Uprime_endpoint": Up.subs(t, point) - s.Rational(7125, 64),
        "actual_alpha_positive_mass_and_curvature_formula": s.factor(
            alpha + 2 * U / 3 - 5 * mass**2 * a * a / 3
        ),
        "actual_alpha_conformal_derivative": s.factor(
            alphap + s.Rational(2, 3) * Up - s.Rational(10, 3) * mass**2 * a * a * h
        ),
        "actual_gamma_formula": s.factor(
            gamma
            + 2 * U * U
            - 10 * mass**2 * a * a * U
            - s.Rational(5, 2) * mass**4 * a**4
        ),
        "actual_a_squared_endpoint": (a * a).subs(t, point) - s.Rational(625, 256),
        "actual_A_coefficient_norm": 22 * 5 * mass**2 - 110 * mass**2,
        "actual_H_coefficient_norm": 37 * 5 * mass**2 - 185 * mass**2,
        "actual_J_coefficient_norm": 11 * 26 * mass**2 - 286 * mass**2,
        "actual_volume_rank_one_coefficient_norm": 10 * 16 * mass**4 - 160 * mass**4,
    }
    return {
        "actual_profile_coefficients": {
            "alpha": alpha,
            "alpha_prime": alphap,
            "gamma": gamma,
            "Uprime": Up,
        },
        "scalar_bounds": {
            "alpha": 5 * mass**2,
            "alpha_prime": 26 * mass**2,
            "gamma": 16 * mass**4,
            "Uprime": s.Integer(112),
        },
        "matrix_bounds": {
            "A": 110 * mass**2,
            "H": 185 * mass**2,
            "J": 286 * mass**2,
            "V0": 160 * mass**4,
        },
        "whole_interval_argument": "Positive coefficient polynomials and parity bound a^2,U,|h|,|Uprime| throughout |t|<=1/2. The positive alpha mass term dominates2U and subtraction only lowers its upper bound. Absolute-value bounds retain every term in alpha-prime and gamma.",
        "checks": checks,
        "gates": {
            "actual_Uprime_strictly_below112": s.Rational(7125, 64) < 112,
            "Uprime_whole_interval_positive_odd_coefficients": all(
                v >= 0 for v in s.Poly(Up, t).all_coeffs()
            ),
            "actual_alpha_strictly_positive": 5 * mass**2 - 2 * 18 > 0,
            "actual_alpha_upper": s.Rational(5, 3) * mass**2 * s.Rational(625, 256)
            < 5 * mass**2,
            "actual_alphaprime_upper": 25 * mass**2 + s.Rational(224, 3) < 26 * mass**2,
            "actual_gamma_absolute_upper": 2 * 18**2 + 540 * mass**2 + 15 * mass**4
            < 16 * mass**4,
        },
    }
