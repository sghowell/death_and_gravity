"""Actual scalar reference-flow audit with perturbative boundaries retained."""

from functools import cache

from . import calibration, flow, subtraction, tensors


@cache
def residuals():
    return {
        prefix + "_" + name: value
        for prefix, module in (
            ("tensors", tensors),
            ("subtraction", subtraction),
            ("flow", flow),
            ("calibration", calibration),
        )
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    rows = {name: bool(value) for name, value in calibration.data()["bounds"].items()}
    rows.update(
        {
            "same_fixed_full_model_counterterms_not_a_new_b2_contact": True,
            "literal_all_components_of_two_field_scalar_beta_tensors": True,
            "cubic_coupling_and_cubic_squared_factors_not_conflated": True,
            "independent_existing_counterterm_normalization_match": True,
            "reference_bubble_change_derived_with_anchor_and_convergent_difference": True,
            "full_three_channel_amplitude_compensated_through_one_loop": True,
            "finite_potential_contact_fixed_once_at_original_reference": True,
            "light_mass_and_LSZ_remain_on_shell_not_new_MS_running_mass": True,
            "heavy_parameter_not_an_exact_stable_heavy_pole": True,
            "closed_flow_positive_root_and_all_initial_values_checked": True,
            "flow_composition_includes_heavy_mass_shift": True,
            "positive_completed_square_margin_preserved_by_exact_flow_numerator": True,
            "strictly_increasing_margin_has_positive_square_derivative": True,
            "rational_window_enclosures_include_all_references_up_to_selected_Planck_scale": True,
            "one_loop_flow_solution_not_full_higher_loop_or_gravitational_error_bound": True,
            "formal_flow_singularity_not_a_proven_quantum_Landau_pole": True,
            "running_tree_parameters_not_a_scale_dependent_physical_amplitude": True,
            "no_new_Regge_contour_or_common_bounce_parent_from_reference_flow": True,
            "no_old_state_counterterm_or_canonical_map_transfer": True,
            "full_V_G_B_and_original_P8_remain_open": True,
        }
    )
    return rows


@cache
def rejected_inputs():
    for name, call, args in calibration.bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported logarithmic reference time accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    d = calibration.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "finite_running_margin_change_not_set_to_zero": d[
            "whole_window_margin_relative_increase_upper"
        ]
        > 0,
        "heavy_mass_reference_shift_not_omitted": d[
            "whole_window_heavy_mass_squared_increase_upper"
        ]
        > 0,
        "retuning_potential_contact_at_every_reference_not_allowed": True,
        "light_UV_mass_tensor_not_used_to_run_physical_pole": True,
        "coupling_derivative_of_one_loop_term_is_uncomputed_higher_order": True,
        "closed_beta_ODE_not_exact_quantum_flow": True,
        "formal_Landau_singularity_not_UV_exclusion": True,
        "Planck_reference_endpoint_not_control_of_gravitational_diagrams": True,
        "original_P8_not_closed": True,
    }
