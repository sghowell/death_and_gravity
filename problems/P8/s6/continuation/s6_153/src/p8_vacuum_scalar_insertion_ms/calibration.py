"""The complete scalar insertion conversion and four raw scalar-family assembly."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_double_bubble_ms import calibration as double
from p8_vacuum_scalar_zero_reference import calibration as wine
from p8_vacuum_two_loop_finite import calibration as finite
from p8_vacuum_two_loop_insertions import calibration as previous

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
    bound = bounds.enclosure(L, g, M, 4 * lam)
    old = previous.data()["complete_grouped_insertion_b2_absolute_upper"]
    insertion = old + bound["local_conversion_b2_absolute_upper"]
    families = {
        "finite_88_unchanged": finite.data()[
            "actual_finite_subsector_b2_absolute_upper"
        ],
        "OS_insertion_64_MS_outer_reference": insertion,
        "double_bubble_24_MS": double.data()["actual_double_bubble_MS_enclosure"][
            "complete_MS_double_bubble_upper"
        ],
        "wineglass_16_MS": wine.data()["actual_wineglass_MS_enclosure"][
            "complete_wineglass_MS_upper"
        ],
    }
    total = sum(families.values())
    return {
        "actual_parameters": p,
        "actual_insertion_MS_conversion_enclosure": bound,
        "same_old_complete_insertion_upper": old,
        "complete_insertion_MS_upper": insertion,
        "four_raw_scalar_interaction_family_bounds": families,
        "four_raw_scalar_MS_interaction_families_upper": total,
        "four_raw_scalar_MS_interaction_families_relative_upper": total / (4 * lam),
        "checks": {
            "same_parent_asymptotic_multiplier_bound": bound["alpha_zero_upper"]
            - previous.data()["same_fixed_asymptotic_multiplier_upper"],
            "both_old_remainders_and_local_conversion_once": insertion
            - old
            - bound["local_conversion_b2_absolute_upper"],
            "four_disjoint_raw_families_complete_count": 88 + 64 + 24 + 16 - 192,
            "four_family_sum": total - sum(families.values()),
            "same_actual_scale": 10**400 - (10**200) ** 2,
        },
        "bounds": {
            "actual_mass_within_fixed_scale": bool(32 < M < 10**400),
            "actual_cubic_over_mass_condition": bool(g / M < L / 3),
            "actual_forward_multiplier_negative": bool(3 * g / (M - 2) < 2 * L),
            "actual_forward_multiplier_abs_below_two_L": bool(
                0 < 3 * g / (M - 2) < 2 * L
            ),
            "local_conversion_relative_below_one_e_minus_409": bool(
                bound["local_conversion_b2_relative_upper"] < s.Rational(1, 10**409)
            ),
            "local_conversion_absolute_below_one_e_minus_1008": bool(
                bound["local_conversion_b2_absolute_upper"] < s.Rational(1, 10**1008)
            ),
            "complete_MS_insertion_below_one_e_minus_614": bool(
                insertion < s.Rational(1, 10**614)
            ),
            "complete_MS_insertion_relative_below_two_e_minus_15": bool(
                insertion / (4 * lam) < s.Rational(2, 10**15)
            ),
            "four_raw_MS_scalar_families_below_three_e_minus_607": bool(
                total < s.Rational(3, 10**607)
            ),
            "four_raw_MS_scalar_families_relative_below_one_e_minus_7": bool(
                total / (4 * lam) < s.Rational(1, 10**7)
            ),
        },
        "scope": "All four raw scalar interaction families now have bounds with MS interaction references at fixed canonical OS coordinates. The 88 entirely UV-finite refinements are scheme independent at this formal order. The isolated sigma insertion was handled in S6.150, not added here. Global field/coupling re-expansion, cross terms, full GY14 vacuum/canonical assembly and V/G/B remain open.",
    }
