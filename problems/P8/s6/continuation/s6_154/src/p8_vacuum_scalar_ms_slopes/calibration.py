"""Actual scalar MS slope and converted fixed-order OS pole enclosures."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_two_loop_light_pole import calibration as previous

from . import bounds


@cache
def data():
    p = model.data()["actual_parameters"]
    L, g, M = [
        p[k]
        for k in (
            "bare_polynomial_quartic",
            "cubic_coupling_squared",
            "heavy_mass_squared",
        )
    ]
    enclosure = bounds.enclosure(L, g, M)
    old = previous.data()
    B1 = old["inherited_one_loop_quadratic_coefficient_upper"]
    B2 = old["actual_two_loop_quadratic_coefficient_upper"]
    converted = B2 + enclosure["additional_OS_quadratic_upper"]
    total = B1 + converted
    return {
        "actual_parameters": p,
        "actual_scalar_MS_slope_enclosure": enclosure,
        "unchanged_one_loop_OS_quadratic_upper": B1,
        "old_two_loop_OS_quadratic_upper": B2,
        "converted_two_loop_OS_quadratic_upper": converted,
        "converted_one_plus_two_loop_OS_quadratic_upper": total,
        "unit_disc_factored_inverse_gap_lower": 1 - total,
        "checks": {
            "same_full_old_two_loop_bound": B2
            - sum(old["actual_two_loop_group_quadratic_coefficient_uppers"].values()),
            "same_inner_alpha_bound": enclosure["same_inner_alpha_upper"]
            - old["inherited_asymptotic_multiplier_upper"],
            "converted_bound_retains_all_old_groups": converted
            - B2
            - enclosure["additional_OS_quadratic_upper"],
            "one_plus_two_loop_sum": total - B1 - converted,
            "seven_disjoint_slope_groups": len(
                enclosure["finite_MS_slope_group_uppers"]
            )
            - 7,
            "canonical_factored_gap": 1 - total - (1 - B1 - converted),
        },
        "bounds": {
            "actual_scale_and_mass_domain": bool(32 < M < 10**400),
            "all_slope_groups_positive": all(
                v > 0 for v in enclosure["finite_MS_slope_group_uppers"].values()
            ),
            "local_MS_slope_below_three_e_minus_412": bool(
                enclosure["finite_MS_slope_group_uppers"]["local_sunset_MS"]
                < s.Rational(3, 10**412)
            ),
            "proper_MS_slope_conversion_below_six_e_minus_412": bool(
                enclosure["finite_MS_slope_group_uppers"][
                    "proper_interaction_MS_conversion"
                ]
                < s.Rational(6, 10**412)
            ),
            "complete_scalar_finite_MS_slope_below_one_e_minus_18": bool(
                enclosure["finite_MS_slope_absolute_upper"] < s.Rational(1, 10**18)
            ),
            "additional_OS_conversion_below_two_e_minus_609": bool(
                enclosure["additional_OS_quadratic_upper"] < s.Rational(2, 10**609)
            ),
            "converted_two_loop_OS_upper_below_one_e_minus_18": bool(
                converted < s.Rational(1, 10**18)
            ),
            "converted_one_plus_two_loop_OS_upper_below_one_e_minus_18": bool(
                total < s.Rational(1, 10**18)
            ),
            "strict_unit_disc_fixed_order_inverse_gap": bool(1 - total > 0),
        },
        "scope": "Complete scalar quadratic intermediate interaction-MS reference with inner physical OS fixed. The finite MS slope is bounded before outer OS. Unit residue follows only upon imposing that outer condition, not for the unnormalized MS field. No GY14-wide finite map, global pole count or physical higher-loop remainder is inferred.",
    }
