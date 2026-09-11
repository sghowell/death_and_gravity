"""Actual complete canonical GY14 Phi forward coefficient through two loops."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_four_scalar import calibration as box
from p8_vacuum_fermion_mixed_quartic import bounds as mixed
from p8_vacuum_fermion_outer_ms import calibration as inserted
from p8_vacuum_fermion_vertex_chord import calibration as direct
from p8_vacuum_finite_field_covariance import calibration as field
from p8_vacuum_full_one_loop_matching import calibration as first
from p8_vacuum_full_phi_normalization import calibration as phi
from p8_vacuum_scalar_insertion_ms import calibration as scalar2

from . import bounds


@cache
def data():
    d = first.data()
    p = d["actual_MS_interaction_parameters"]
    L, g, M, lam = [p[k] for k in ("L", "g", "M", "lambda")]
    tree = 4 * lam
    ell = d["exact_scale_log_rational_enclosure"]["scale_log_upper"]
    scalar_first = d["old_scalar_complete_one_loop_absolute_error"] + tree * abs(
        -2 * L + 3 * g / (M - 2)
    ) * ell / (2 * 144)
    fermion_first = tree * box.data()["rational_box_second_coefficient_relative_upper"]
    k = field.data()["actual_first_finite_field_enclosure"][
        "first_normalization_zero_absolute_upper"
    ]
    groups = {
        "all_192_scalar_MS_interaction_refinements": scalar2.data()[
            "four_raw_scalar_MS_interaction_families_upper"
        ],
        "both_direct_fermionic_quartic_rows": direct.data()[
            "both_complete_quartic_primitive_upper"
        ],
        "complete_fermionic_W2_F0_row": inserted.data()[
            "full_family_b2_upper_in_common_MS_interaction_scheme"
        ],
        "complete_fermionic_W1_F2_row": mixed.data()["enclosure"][
            "complete_paired_mixed_b2_upper"
        ],
    }
    raw = sum(groups.values())
    cross = bounds.first_variation_upper(k, scalar_first, fermion_first)
    pole = phi.data()
    tree_field = pole["finite_field_map_tree_b2_second_relative_upper"]
    E1 = d["complete_one_loop_absolute_error_upper"]
    enclosure = bounds.assemble(tree, E1, raw, cross, tree_field)
    B = pole["complete_canonical_one_plus_two_loop_OS_coefficient_upper"]
    return {
        "actual_MS_interaction_parameters": p,
        "tree_b2": tree,
        "complete_canonical_one_loop_absolute_upper": E1,
        "finite_hybrid_one_loop_group_uppers": {
            "scalar": scalar_first,
            "fermion": fermion_first,
        },
        "fixed_hybrid_two_loop_group_uppers": groups,
        "fixed_hybrid_two_loop_absolute_upper": raw,
        "first_coordinate_cross_absolute_upper": cross,
        "second_tree_field_map_relative_upper": tree_field,
        "complete_canonical_two_loop_enclosure": enclosure,
        "unchanged_complete_canonical_unit_disc_OS_coefficient_upper": B,
        "checks": {
            "same_positive_tree_reference": tree - 2 * g / (M - 2) ** 3,
            "all_raw_groups_once": raw - sum(groups.values()),
            "both_first_hybrid_groups_with_correct_weights": cross
            - k * (4 * scalar_first + 2 * fermion_first),
            "complete_second_absolute_assembly": enclosure["second_absolute"]
            - raw
            - cross
            - tree * tree_field,
            "complete_relative_assembly": enclosure["total_relative"]
            - (E1 + enclosure["second_absolute"]) / tree,
            "no_extra_first_canonical_field_factor": E1
            - d["complete_one_loop_absolute_error_upper"],
            "unchanged_independent_complete_Phi_pole": B
            - pole["complete_canonical_one_plus_two_loop_OS_coefficient_upper"],
            "forward_scalar_MS_reference_factor": -2 * L
            + 3 * g / (M - 2)
            + g * (3 * (M - 2) - 4) / (M - 2) ** 2,
        },
        "bounds": {
            "all_raw_group_uppers_strictly_positive": all(
                v > 0 for v in groups.values()
            ),
            "hybrid_first_group_uppers_strictly_positive": bool(
                scalar_first > 0 and fermion_first > 0
            ),
            "tree_coefficient_strictly_positive": bool(tree > 0),
            "first_coordinate_cross_relative_below_one_e_minus_212": bool(
                cross / tree < s.Rational(1, 10**212)
            ),
            "second_tree_field_relative_below_one_e_minus_18": bool(
                tree_field < s.Rational(1, 10**18)
            ),
            "full_two_loop_relative_below_one_e_minus_7": bool(
                enclosure["second_relative"] < s.Rational(1, 10**7)
            ),
            "full_one_plus_two_loop_relative_below_one_e_minus_6": bool(
                enclosure["total_relative"] < s.Rational(1, 10**6)
            ),
            "formal_through_two_loop_b2_strictly_positive": bool(
                enclosure["formal_uniform_lower"] > 0
            ),
            "same_complete_unit_disc_Phi_pole_gap": bool(B < s.Rational(1, 10**18)),
        },
        "scope": "The complete formal one-plus-two-loop Phi forward coefficient in the named nongravitational GY14 reference model is positive. Physical truncation, vacuum/source references, V contours/cuts, G and B are not inferred.",
    }
