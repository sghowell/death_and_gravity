"""Actual double-bubble MS enclosure with all three terms retained."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_two_loop_double_bubble import calibration as previous

from . import bounds


@cache
def data():
    p = model.data()["actual_parameters"]
    L, g, M, lam = (
        p[k]
        for k in (
            "bare_polynomial_quartic",
            "cubic_coupling_squared",
            "heavy_mass_squared",
            "lambda",
        )
    )
    d = bounds.enclosure(L, 4 * lam)
    return {
        "actual_parameters": p,
        "actual_scale_logarithm": 400 * s.log(10),
        "actual_double_bubble_MS_enclosure": d,
        "checks": {
            "same_old_assigned_family_bound": d["old_complete_double_bubble_upper"]
            - previous.data()["actual_double_bubble_group_b2_absolute_upper"],
            "full_MS_family_three_term_sum": d["complete_MS_double_bubble_upper"]
            - d["old_complete_double_bubble_upper"]
            - d["full_linear_scale_term_upper"]
            - d["quadratic_scale_b2_absolute_upper"],
            "same_actual_tree": 4 * lam - previous.data()["actual_tree_b2"],
        },
        "bounds": {
            "actual_heavy_mass_above_thirty_two": bool(M > 32),
            "actual_external_heavy_gap": bool(M - 3 > M / 2),
            "actual_cubic_over_mass_below_quartic_third": bool(g / M < L / 3),
            "actual_forward_cubic_over_D_below_quartic_half": bool(g / (M - 2) < L / 2),
            "actual_logarithm_argument_below_ten_to_two_hundred": bool(4 * M < 10**200),
            "quadratic_scale_relative_below_one_e_minus_405": bool(
                d["quadratic_scale_b2_relative_upper"] < s.Rational(1, 10**405)
            ),
            "linear_scale_below_seven_e_minus_612": bool(
                d["full_linear_scale_term_upper"] < s.Rational(7, 10**612)
            ),
            "complete_MS_double_bubble_below_seven_e_minus_612": bool(
                d["complete_MS_double_bubble_upper"] < s.Rational(7, 10**612)
            ),
            "complete_MS_relative_below_two_e_minus_12": bool(
                d["complete_MS_double_bubble_relative_upper"] < s.Rational(2, 10**12)
            ),
        },
        "scope": "The complete 24-refinement assigned double-bubble interaction forest in MS at mu=mF, not the complete scalar/GY14 amplitude or its finite field/parameter-map assembly.",
    }
