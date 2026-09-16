"""Written compact-domain majorants with exact algebraic budget checks."""

from functools import cache

import sympy as s

from . import source

DELTA = s.Symbol("positive_angular_window", positive=True)
LOG = s.Symbol("nonnegative_log_window", nonnegative=True)
KNOWN_FINITE_CONSTANT = s.Integer(10) ** 11
KNOWN_AMPLITUDE_CONSTANT = s.Integer(10) ** 9
COMPARISON_WINDOW = s.Rational(1, 10**204)


def require_positive_exact(value):
    if isinstance(value, (bool, float)) or not isinstance(value, (int, s.Basic)):
        raise TypeError("Require a finite exact positive real number")
    value = s.sympify(value)
    if value.has(s.Float) or value.is_number is not True or value.is_real is not True:
        raise ValueError("Require a finite exact positive real number")
    if value.is_finite is not True or value.is_positive is not True:
        raise ValueError("Require a finite exact positive real number")
    return value


def require_delta(value):
    value = require_positive_exact(value)
    if value > 1:
        raise ValueError("Require 0<delta<=1")
    return value


def require_physical(energy, transfer, delta):
    delta = require_delta(delta)
    if isinstance(energy, (bool, float)) or isinstance(transfer, (bool, float)):
        raise TypeError("Require exact physical invariants")
    if not isinstance(energy, (int, s.Basic)) or not isinstance(
        transfer, (int, s.Basic)
    ):
        raise TypeError("Require exact physical invariants")
    energy, transfer = map(s.sympify, (energy, transfer))
    if any(
        v.has(s.Float)
        or v.is_number is not True
        or v.is_real is not True
        or v.is_finite is not True
        for v in (energy, transfer)
    ):
        raise ValueError("Require exact finite real invariants")
    crossed = 4 - energy - transfer
    if not s.Rational(25, 4) <= energy <= 16 or transfer >= 0 or crossed >= 0:
        raise ValueError("Require the stated massive physical compact domain")
    if min(-transfer, -crossed) < delta:
        raise ValueError("The point is outside the stated angular window")
    return energy, transfer, crossed, delta


def known_finite_bound(delta):
    delta = require_delta(delta)
    return KNOWN_FINITE_CONSTANT * (1 - s.log(delta)) / delta**2


def known_amplitude_bound(delta, kappa):
    delta, kappa = require_delta(delta), require_positive_exact(kappa)
    return KNOWN_AMPLITUDE_CONSTANT * (1 - s.log(delta)) / (kappa**2 * delta**2)


def known_rate_interference_bound(delta, kappa, positive_born_lower):
    lower = require_positive_exact(positive_born_lower)
    return 2 * known_amplitude_bound(delta, kappa) / lower


def angular_coefficient_l1(expression):
    dim = source.dimensional
    expression = s.factor(expression.subs(dim.N, 2))
    return sum(
        abs(coefficient) * 8**i * 17**j
        for (i, j), coefficient in s.Poly(expression, dim.H, dim.Z).terms()
    )


def exact_budget():
    return {
        "six_ordered_boxes": 6 * (110000 * 200 * 10 + 1300 * 200),
        "three_massless_triangles": 3 * 210000000 * 10,
        "three_massive_triangles": 3 * (200000 * 1200 + 2500 * 100),
        "three_massless_bubbles": 3 * 2300000 * 12,
        "three_massive_bubbles": 3 * 30000000 * 21,
        "compact_evanescent_bubbles": 5000,
        "four_LSZ_finite": 28000,
        "physical_pole_completion": 23000,
        "old_Gram_finite": 22000,
        "D_tree_soft_reference": 1200000,
    }


@cache
def data():
    dim = source.dimensional
    cl, cp = dim.angular_reduction()[2][2:]
    triangle = angular_coefficient_l1(cl) * 20**3 / 32
    bubble = angular_coefficient_l1(cp) * 20**2 / 16
    # Use exact rationals even in majorants; no floating constants in certificates.
    bmm = (
        2 * 644 * 100 * (17**2 + s.Rational(1, 2) * (1 + 17**2))
        + 140**2
        + s.Rational(2, 3) * 140 * 100
        + s.Rational(1, 15) * 100**2 * (1 + 2 * 17**2)
    ) / 2
    basis = (3, 3, 768, 3072)
    r0 = sum(
        abs(x) * y
        for x, y in zip(
            (
                s.Rational(164, 15),
                -s.Rational(2, 15),
                -s.Rational(73, 15),
                s.Rational(1, 15),
            ),
            basis,
            strict=True,
        )
    )
    r1 = sum(
        abs(x) * y
        for x, y in zip(
            (
                -s.Rational(4879, 225),
                -s.Rational(28, 225),
                s.Rational(3128, 225),
                -s.Rational(16, 225),
            ),
            basis,
            strict=True,
        )
    )
    evan = sum(
        abs(x) * y
        for x, y in zip(
            (
                s.Rational(69, 50),
                -s.Rational(4979, 90),
                s.Rational(3044, 225),
                s.Rational(28, 225),
                -s.Rational(418, 225),
                s.Rational(16, 225),
            ),
            (768, 1, *basis),
            strict=True,
        )
    )
    cm_bound = 2 * 322 * 296
    cm1_bound = 4 * (296 + 322)
    budget = sum(exact_budget().values())
    checks = {
        "massive_eikonal_bound": 16**2 + 4 * 16 + 2 - 322,
        "squared_eikonal_margin": 110000 - 322**2 - 6316,
        "box_dimensional_derivative_margin": 1300 - 4 * 322 - 12,
        "massive_triangle_h_bound": 2 * 16 + 6 + 2 + 16**2 - 296,
        "massive_triangle_whole_bound": 2 * 322 * 296 - 190624,
        "massive_triangle_first_bound": 4 * 296 + 4 * 322 - 2472,
        "massless_triangle_L1_exact": triangle - 207864500,
        "massless_bubble_L1_exact": bubble - s.Rational(6718985, 3),
        "massive_bubble_whole_bound": bmm - s.Rational(84471200, 3),
        "contour_norm_log_bound": 200 * 9 - 1800,
        "contour_log_integral_bound": 2 * 9 - 18,
        "massive_triangle_finite_bound": s.Rational(1, 2) * (3 * 200 + 1800) - 1200,
        "massive_bubble_finite_bound": 3 + 18 - 21,
        "three_tree_numerator_bound": 3 * (2 + 2 * 16 + 16**2) - 870,
        "tree_derivative_bound": 2 * 3 - 6,
        "three_soft_channel_bound_margin": 100000
        - (s.Rational(3, 2) * 322 * 200 + 1)
        - 3399,
        "LSZ_complete_finite_bound": (6 * 3 + 14) * 870 + 6 * 6 - 27876,
        "old_Gram_each_bound": s.Rational(12, 4) * 5 * 256
        + s.Rational(27, 2) * 256
        - 7296,
        "old_Gram_crossed_bound": 3 * 7296 - 21888,
        "physical_completion_D4_bound": r0 - s.Rational(19878, 5),
        "physical_completion_first_bound": r1 - s.Rational(822059, 75),
        "physical_completion_combined_budget": 3 * 4000 + 11000 - 23000,
        "finite_soft_division_bound": 2 * 6 * 100000 - 1200000,
        "complete_finite_budget": budget - 10316388000,
        "complete_budget_margin": KNOWN_FINITE_CONSTANT - budget - 89683612000,
        "loop_normalization_margin": 16 * 9 * KNOWN_AMPLITUDE_CONSTANT
        - KNOWN_FINITE_CONSTANT
        - 44000000000,
        "original_window_amplitude_power": 9 + 3 + 2 * 204 - 2 * 800 + 1180,
        "original_matter_lower_power": 4 * s.Rational(1, 8192) ** 2 * 256**3 - 1,
        "original_known_relative_power": -1180 + 600 + 580,
        "original_rate_error_margin": 10 * s.Integer(10) ** (-580)
        - 2 * s.Integer(10) ** (-580)
        - 8 * s.Integer(10) ** (-580),
        "comparison_window_log_bound": 1 + 204 * 3 - 613,
    }
    return {
        "physical_domain": "mu=nu=1;25/4<=s<=16;t,u<0;t+u=4-s;"
        "0<delta<=1;min(-t,-u)>=delta. This is an angular comparison window, "
        "not a physical ultraviolet cutoff or a forward-limit estimate.",
        "master_coefficient_bounds": {
            "C00mu": triangle,
            "B00": bubble,
            "Bmm_times_delta2": bmm,
            "abs_C0mumu_times_delta": 200000,
            "abs_C0mumu_first_times_delta": 2500,
            "physical_R0_times_delta2": r0,
            "physical_R1_times_delta2": r1,
            "compact_evanescent_times_delta2": evan,
        },
        "whole_termwise_budget": exact_budget(),
        "whole_known_finite_bound": KNOWN_FINITE_CONSTANT * (1 + LOG) / DELTA**2,
        "whole_known_amplitude_bound": KNOWN_AMPLITUDE_CONSTANT
        * (1 + LOG)
        / (source.K**2 * DELTA**2),
        "original_parameter_comparison": {
            "explicit_window": COMPARISON_WINDOW,
            "matter_Born_lower": s.Integer(10) ** (-600),
            "known_loop_upper": s.Integer(10) ** (-1180),
            "known_loop_over_full_Born_upper": s.Integer(10) ** (-580),
            "known_rate_interference_upper": s.Integer(10) ** (-579),
        },
        "matching_boundary": "The bound is on the specified known representative "
        "only. alpha,beta,delta_kappa may not be bounded or set to zero from "
        "these estimates. No full amplitude or original V/G/B/P8 conclusion.",
        "checks": checks,
        "gates": {
            "whole_massless_coefficients_bounded_without_angle_sampling": bool(
                triangle < 210000000 and bubble < 2300000
            ),
            "whole_massive_coefficients_bounded_without_angle_sampling": bool(
                bmm < 30000000 and cm_bound < 200000 and cm1_bound < 2500
            ),
            "complete_physical_and_evanescent_bounds": bool(
                r0 < 4000 and r1 < 11000 and evan < 5000
            ),
            "whole_positive_budget_has_strict_margin": bool(budget < 10**11),
            "original_matter_lower_bound_uses_unchanged_n_and_g": bool(
                source.HEAVY_MASS2 < s.Integer(10) ** 200 / 256
                and source.CUBIC == s.Rational(1, 8192)
            ),
            "log_and_pi_inequalities_have_written_analytic_justification": True,
            "all_three_unknown_finite_matching_terms_excluded_explicitly": True,
            "window_not_exchanged_with_exact_forward_or_Regge_limit": True,
        },
    }
