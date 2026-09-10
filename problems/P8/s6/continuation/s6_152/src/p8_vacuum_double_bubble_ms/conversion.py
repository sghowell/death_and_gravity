"""Full finite scale shift of all selected factorizing refinements."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_double_bubble import forests as previous


@cache
def data():
    C, IR, TA, TB, ell, Q = s.symbols("C finite_I_R full_T_A full_T_B ell Q")
    d = ell / Q
    old = (C**3 * IR**2 + C**2 * IR * (TA + TB)) / 4
    new = (C**3 * (IR + d) ** 2 + C**2 * (IR + d) * (TA + TB)) / 4
    linear = (2 * d * C**3 * IR + d * C**2 * (TA + TB)) / 4
    square = d**2 * C**3 / 4
    left, right = s.symbols("left_triangle_sum right_triangle_sum")
    raw = previous.data()["selected_raw_amplitude_per_channel"]
    return {
        "old_selected_channel": old,
        "MS_selected_channel": new,
        "complete_linear_scale_term": linear,
        "complete_quadratic_scale_term": square,
        "same_parent_raw_selected_channel": raw,
        "checks": {
            "full_MS_minus_old_two_distinct_scale_terms": s.expand(
                new - old - linear - square
            ),
            "complete_left_right_factorization": s.expand(
                C * ((C * (IR + d) + left) * (C * (IR + d) + right) - left * right) / 4
                - new.subs({TA: left, TB: right})
            ),
            "left_full_triangle_coefficient": s.diff(new, TA) - C**2 * (IR + d) / 4,
            "right_full_triangle_coefficient": s.diff(new, TB) - C**2 * (IR + d) / 4,
            "excluded_finite_triangle_product_unchanged_by_scale": s.diff(
                C * TA * TB / 4, ell
            ),
            "zero_scale_change_recovers_old": new.subs(ell, 0) - old,
            "linear_scale_derivative": s.diff(new, ell).subs(ell, 0) * ell - linear,
            "quadratic_scale_derivative": s.diff(new, ell, 2) * ell**2 / 2 - square,
        },
        "scope": "The exact old family with I_R replaced by I_R+ell/Q after the complete MS forest cancellation. Both heavy triangles and full external heavy poles are retained. The separately owned finite triangle-triangle terms are not counted again. No local reference or first-order product is added a second time.",
    }
