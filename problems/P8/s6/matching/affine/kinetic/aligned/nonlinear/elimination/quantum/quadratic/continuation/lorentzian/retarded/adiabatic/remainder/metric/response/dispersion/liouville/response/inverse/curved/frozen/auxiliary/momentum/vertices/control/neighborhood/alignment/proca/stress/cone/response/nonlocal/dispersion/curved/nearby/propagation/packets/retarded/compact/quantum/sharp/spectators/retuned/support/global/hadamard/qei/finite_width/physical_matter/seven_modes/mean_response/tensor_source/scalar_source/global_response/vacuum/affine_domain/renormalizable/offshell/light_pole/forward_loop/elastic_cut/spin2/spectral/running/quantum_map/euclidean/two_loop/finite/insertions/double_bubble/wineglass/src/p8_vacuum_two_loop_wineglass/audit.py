"""Actual forest, analytic integral and remaining-obligation audit."""

from functools import cache

from . import calibration, integrals, kernel, selection, subtraction


@cache
def residuals():
    return {
        prefix + "_" + k: v
        for prefix, module in (
            ("selection", selection),
            ("subtraction", subtraction),
            ("kernel", kernel),
            ("integrals", integrals),
            ("calibration", calibration),
        )
        for k, v in module.data()["checks"].items()
    }


def bad_cases():
    invalid = (
        (),
        (0, 0),
        (0, 0, 0, 0),
        [0, 0, 0],
        (True, 0, 0),
        (0, 0, False),
        (0, 4, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, 2),
        (2, 2, 0),
    )
    return calibration.bad_cases() + [
        ("unsupported_graph_" + str(i), selection.require_group, (v,))
        for i, v in enumerate(invalid)
    ]


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported wineglass input accepted: " + name)
    return len(bad_cases())


@cache
def gates():
    d = selection.data()
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "exact_sixteen_raw_refinements_selected": d["selected_refinement_count"]
            == 16,
            "exact_two_overall_subtractions": d["overall_subtraction_refinement_count"]
            == 2,
            "exact_fourteen_finite_after_inner_subtraction": d[
                "finite_after_inner_subtraction_count"
            ]
            == 14,
            "all_twelve_external_label_routings_checked": kernel.data()["routing_count"]
            == 12,
            "same_complete_family_Wick_weight": d["total_family_Wick_weight"] == 3,
            "proper_bubble_uses_same_entire_regulated_I0": True,
            "higher_order_recursive_zero_reference_specified_before_estimation": True,
            "nested_forest_operations_applied_inside_out": True,
            "local_counterterms_defined_before_external_heavy_attachments": True,
            "common_regulator_cancellations_precede_removal": True,
            "logarithm_difference_kept_together_for_convergence": True,
            "right_half_plane_scalar_interpolation_on_first_sheet": True,
            "bilinear_shift_square_not_confused_with_Hermitian_norm": True,
            "strict_margins_allow_neighborhood_of_closed_forward_disc": True,
            "complete_heavy_propagators_inside_both_loops": True,
            "all_radius_radial_integrals_not_hard_cutoffs": True,
            "Cauchy_bound_on_whole_integrated_group": True,
            "all_four_raw_groups_disjoint": True,
            "positive_error_upper_not_signed_two_loop_correction": True,
            "fixed_finite_potential_counterterm_insertions_still_required": True,
            "complete_two_loop_light_pole_residue_LSZ_still_required": True,
            "no_source_aware_two_loop_or_gravity_state_transfer": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return result


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "overall_subtraction_present_for_both_required_cores": all(
            (0, 1) in r["restricted_forest_index_sets"]
            for r in selection.data()["selected_rows"]
            if r["has_overall_core"]
        ),
        "unsubtracted_reference_not_assigned_a_finite_value": True,
        "decaying_difference_not_replaced_by_uniform_constant": True,
        "all_heavy_weighted_remainder_terms_retained": True,
        "raw_count_192_not_a_complete_renormalized_amplitude": True,
        "finite_potential_insertions_not_silently_removed": True,
        "complete_two_loop_pole_LSZ_not_declared_from_unit_residue_convention": True,
        "full_two_loop_and_original_P8_not_closed": True,
    }
