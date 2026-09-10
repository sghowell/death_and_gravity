"""Actual necessary matter-cone relation and analytic full-domain tensor bounds."""

from functools import cache

import sympy as sp

from . import family


@cache
def data():
    X = family.X
    d = family.data()
    B = family.B
    R, RX = d["R"], d["RX"]
    h = d["h"]
    theta = sp.Symbol("nonzero_regular_Theta", nonzero=True, real=True)
    Y = sp.Symbol("positive_rolling_matter_Y", positive=True)

    def mismatch(a3):
        delta = X * (2 * RX - X * a3) / (2 * R)
        f = -(R - 2 * X * RX + R * delta) / theta
        g = -R * (1 - 3 * delta) / theta
        return sp.factor(f - g)

    good = mismatch(d["functions_on_transition_or_tube"]["A3"])
    naive = family.B / (h * X)
    wrong = mismatch(naive)
    determinant = -Y * wrong**2 / 4
    z = sp.Symbol("transition_argument", real=True)
    r1, r2 = sp.symbols("positive_left_bump positive_right_bump", positive=True)
    derivative = (
        r1
        * r2
        / (r1 + r2) ** 2
        * ((z - sp.Rational(1, 4)) ** -2 + (sp.Rational(3, 4) - z) ** -2)
    )
    left = sp.exp(-1 / (z - sp.Rational(1, 4)))
    right = sp.exp(-1 / (sp.Rational(3, 4) - z))
    literal_switch = left / (left + right)

    # Substituting X and the switch simultaneously would leave the old
    # switch argument in derivative terms; replace switch jets first.
    def at_center(value):
        return sp.factor(
            value.subs(
                {family.B: sp.Rational(1, 2), sp.diff(B, X): 8}, simultaneous=True
            ).subs({family.u: 0, X: sp.Rational(1, 2)})
        )

    return {
        "regular_chart_Theta": theta,
        "rolling_matter_Y": Y,
        "correct_exceptional_mismatch": good,
        "naive_independently_tapered_A3": naive,
        "naive_exceptional_mismatch": wrong,
        "naive_actual_free_M1_difference_determinant": sp.factor(determinant),
        "positive_switch_derivative_formula": derivative,
        "transition_center_R": at_center(R),
        "transition_center_RX": at_center(RX),
        "transition_center_correct_A3": at_center(
            d["functions_on_transition_or_tube"]["A3"]
        ),
        "transition_center_naive_A3": at_center(naive),
        "transition_center_naive_difference_determinant": at_center(determinant),
        "all_real_X_and_u_tensor_factor_lower": sp.Rational(83, 164),
        "all_real_X_and_u_Ia_denominator_absolute_lower": sp.Rational(83, 328),
        "three_region_lower_bounds": (
            sp.Rational(5, 8),
            sp.Rational(83, 164),
            sp.Rational(3, 5),
        ),
        "middle_switch_odds_exponent": sp.Rational(80, 21),
        "middle_switch_odds_upper": sp.Integer(81),
        "middle_switch_upper": sp.Rational(81, 82),
        "checks": {
            "literal_bump_switch_derivative": sp.factor(
                sp.diff(literal_switch, z) - derivative.subs({r1: left, r2: right})
            ),
            "literal_bump_switch_center": literal_switch.subs(z, sp.Rational(1, 2))
            - sp.Rational(1, 2),
            "literal_bump_switch_center_derivative": sp.diff(literal_switch, z).subs(
                z, sp.Rational(1, 2)
            )
            - 8,
            "correct_relation_zero": good,
            "naive_source_relation_keeps_curvature_switch_derivative": sp.factor(
                wrong + 2 * X * (X - 1) * sp.diff(B, X) / (h * theta)
            ),
            "naive_actual_free_M1_determinant_is_strict_negative_square": sp.factor(
                determinant
                + Y * X**2 * (X - 1) ** 2 * sp.diff(B, X) ** 2 / (h * h * theta * theta)
            ),
            "exact_middle_exponent_from_switch_definition": sp.Rational(1, 1)
            / (sp.Rational(3, 4) - sp.Rational(3, 5))
            - 1 / (sp.Rational(3, 5) - sp.Rational(1, 4))
            - sp.Rational(80, 21),
            "exact_first_region_lower_bound": 1
            - sp.Rational(1, 2) * sp.Rational(3, 4)
            - sp.Rational(5, 8),
            "exact_middle_region_lower_bound": 1
            - sp.Rational(81, 82) * sp.Rational(1, 2)
            - sp.Rational(83, 164),
            "exact_midpoint_R": at_center(R) - sp.Rational(3, 4),
            "exact_midpoint_RX": at_center(RX) + sp.Rational(7, 2),
            "exact_midpoint_A3": at_center(d["functions_on_transition_or_tube"]["A3"])
            + 7,
            "exact_midpoint_naive_A3": at_center(naive) - 1,
            "exact_midpoint_naive_determinant": at_center(determinant)
            + 4 * Y / theta**2,
        },
        "bounds": {
            "middle_exponent_below_four": sp.Rational(80, 21) < 4,
            "global_floor_is_below_each_region": all(
                sp.Rational(83, 164) <= v
                for v in (sp.Rational(5, 8), sp.Rational(83, 164), sp.Rational(3, 5))
            ),
            "tensor_factor_floor_is_strictly_positive": sp.Rational(83, 164) > 0,
            "Ia_denominator_floor_is_positive": sp.Rational(83, 328) > 0,
        },
    }
