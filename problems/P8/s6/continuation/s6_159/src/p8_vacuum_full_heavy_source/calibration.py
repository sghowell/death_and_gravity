"""Complete stationary-source bound at the unchanged GY14 boundary."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_insertion_ms import calibration as fermion
from p8_vacuum_finite_field_covariance import calibration as field
from p8_vacuum_full_one_loop_matching import calibration as first
from p8_vacuum_full_vacuum_reference import calibration as vacuum

from . import bounds


@cache
def data():
    p = first.data()["actual_MS_interaction_parameters"]
    L, g, M, m = p["L"], p["g"], p["M"], p["mF"]
    G = s.Rational(1, 8192)
    scalar = bounds.scalar_source_enclosure(m, M, L, G, 144, 10**25, 10**13)
    fs = field.data()["actual_first_finite_field_enclosure"]
    ell = first.data()["exact_scale_log_rational_enclosure"]["scale_log_upper"]
    cross = bounds.field_reexpression_upper(
        G,
        fs["first_normalization_zero_absolute_upper"],
        fs["first_normalization_epsilon_absolute_upper"],
        ell,
        144,
    )
    F = fermion.data()["assigned_H_source_absolute_upper"]
    total = scalar["complete_scalar_and_cubic_source_absolute_upper"] + F + cross
    shift = total / M
    mass = G * shift
    return {
        "actual_MS_interaction_parameters": p,
        "actual_fundamental_cubic": G,
        "complete_scalar_and_cubic_source_enclosure": scalar,
        "unchanged_complete_fermion_source_absolute_upper": F,
        "regulated_first_Phi_map_source_absolute_upper": cross,
        "complete_second_H_source_absolute_upper": total,
        "complete_second_H_source_over_M_absolute_upper": shift,
        "source_induced_Phi_mass_reference_absolute_upper": mass,
        "unchanged_complete_second_vacuum_reference_upper": vacuum.data()[
            "complete_matched_two_loop_vacuum_absolute_upper"
        ],
        "checks": {
            "actual_fundamental_cubic_square": G * G - g,
            "same_actual_fermion_mass": m - 10**200,
            "three_complete_source_groups_once": total
            - scalar["complete_scalar_and_cubic_source_absolute_upper"]
            - F
            - cross,
            "fixed_H_shift_reference": M * shift - total,
            "source_induced_Phi_mass_reference": mass - G * total / M,
            "unchanged_old_fermion_source_piece": F
            - fermion.data()["assigned_H_source_absolute_upper"],
            "source2_does_not_change_order2_vacuum_reference": vacuum.data()[
                "complete_matched_two_loop_vacuum_absolute_upper"
            ]
            - vacuum.data()["finite_vacuum_counterterm_second_absolute_upper"],
        },
        "bounds": {
            "actual_positive_reference_parameters": bool(min(L, g, G, M, m) > 0),
            "all_scalar_source_groups_positive": all(
                v > 0 for v in scalar["group_uppers"].values()
            ),
            "scalar_and_cubic_source_bound_below_one_e_80": bool(
                0 < scalar["complete_scalar_and_cubic_source_absolute_upper"] < 10**80
            ),
            "finite_Phi_map_source_bound_below_one_e_minus_209": bool(
                0 < cross < s.Rational(1, 10**209)
            ),
            "complete_second_source_bound_below_one_e_189": bool(0 < total < 10**189),
            "complete_second_source_over_M_below_one_e_minus_8": bool(
                0 < shift < s.Rational(1, 10**8)
            ),
            "source_induced_mass_reference_below_one_e_minus_12": bool(
                0 < mass < s.Rational(1, 10**12)
            ),
            "complete_fixed_order_vacuum_bound_unchanged": bool(
                vacuum.data()["complete_matched_two_loop_vacuum_absolute_upper"]
                < 10**595
            ),
        },
        "scope": "Complete fixed-order stationary-H source reference together with the unchanged full Phi pole, forward coefficient and vacuum reference. These are renormalization-reference values, not a claim that the physical H expectation shifts away from its fixed zero condition. Physical truncation, other coordinate dictionaries and V/G/B remain open.",
    }
