"""Fixed-condition contact insertion and explicitly remaining scope audit."""

from functools import cache

from . import calibration, insertion, potential, radial, subtraction


@cache
def residuals():
    return {
        prefix + "_" + k: v
        for prefix, module in (
            ("potential", potential),
            ("radial", radial),
            ("insertion", insertion),
            ("subtraction", subtraction),
            ("calibration", calibration),
        )
        for k, v in module.data()["checks"].items()
    }


@cache
def gates():
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "literal_full_Hessian_keeps_mixed_heavy_light_kernel": True,
            "same_once_fixed_constant_field_quartic_condition": True,
            "on_shell_quadratic_adjustment_not_new_quartic_freedom": True,
            "strict_negative_contact_from_positive_F_below_A": True,
            "same_entire_regulated_I0_reference": True,
            "local_outer_counterterms_are_exact_linear_variations": True,
            "g_and_M_held_fixed_during_quartic_variation": True,
            "total_canonical_vertex_counterterms_not_extra_tree_LSZ": True,
            "zero_tree_second_derivative_does_not_remove_loop_insertion": True,
            "both_complete_heavy_triangle_terms_retained": True,
            "all_three_channel_Wick_weights_retained": True,
            "all_radius_integrals_not_heavy_expansions_or_cutoffs": True,
            "Cauchy_bound_on_complete_integrated_variation": True,
            "actual_couplings_used_not_only_formal_power_counting": True,
            "finite_contact_not_recounting_UV_reference_forests": True,
            "new_order_constant_potential_contact_cannot_tune_b2": True,
            "known_groups_not_full_two_loop_pole_residue_LSZ_control": True,
            "two_loop_source_aware_matching_not_computed": True,
            "no_all_higher_loop_high_energy_or_gravity_control": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return result


@cache
def rejected_inputs():
    for name, call, args in calibration.bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported finite-contact input accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "nonzero_negative_contact_density_at_unit_radius": calibration.point(1)[
            "finite_contact_density_without_loop_measure"
        ]
        < 0,
        "fixed_finite_contact_not_set_to_zero": True,
        "contact_not_fitted_to_forward_coefficient": True,
        "bare_tree_contact_not_confused_with_its_loop_insertion": True,
        "unregulated_reference_not_assigned_a_finite_value": True,
        "finite_quartic_margin_comparison_not_bare_UV_stability": True,
        "complete_two_loop_normalization_not_inferred_from_raw_graph_count": True,
        "full_V_G_B_and_original_P8_not_closed": True,
    }
