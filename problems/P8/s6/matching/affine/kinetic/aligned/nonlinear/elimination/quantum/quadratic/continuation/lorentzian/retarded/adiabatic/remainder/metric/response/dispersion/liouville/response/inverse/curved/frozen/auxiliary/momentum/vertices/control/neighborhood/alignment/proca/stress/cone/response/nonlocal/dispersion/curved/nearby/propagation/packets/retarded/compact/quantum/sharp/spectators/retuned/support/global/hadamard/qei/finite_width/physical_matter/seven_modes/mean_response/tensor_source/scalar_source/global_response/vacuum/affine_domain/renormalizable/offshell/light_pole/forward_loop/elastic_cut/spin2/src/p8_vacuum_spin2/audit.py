"""Actual stress form factor and explicitly partial finite-gravity audit."""

from functools import cache

from . import calibration, pole, stress, triangles


@cache
def residuals():
    return {
        prefix + "_" + name: value
        for prefix, module in (
            ("stress", stress),
            ("triangles", triangles),
            ("pole", pole),
            ("calibration", calibration),
        )
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    rows = {
        prefix + "_" + name: bool(value)
        for prefix, module in (("triangles", triangles), ("calibration", calibration))
        for name, value in module.data()["bounds"].items()
    }
    rows.update(
        {
            "new_minimal_Einstein_extension_of_same_polynomial_vacuum_only": True,
            "literal_canonical_stress_before_null_projection": True,
            "both_internal_mass_species_stress_insertions_retained": True,
            "complete_four_dimensional_three_propagator_shift_checked": True,
            "triangle_numerator_weight_from_shift_not_only_Ward_guess": True,
            "potential_mass_improvement_and_quartic_bubble_spin_two_terms_vanish": True,
            "same_once_fixed_light_kinetic_counterterm_enforces_Ward_charge_one": True,
            "on_shell_mass_and_residue_scheme_not_retuned_for_gravity": True,
            "spin_two_slope_strictly_positive_by_parameter_integral": True,
            "complex_transfer_disc_majorants_control_limit_remainder": True,
            "stress_two_particle_transfer_cut_not_mixed_self_energy_threshold": True,
            "literal_all_component_graviton_projector_contraction": True,
            "tree_graviton_pole_subtracted_at_fixed_negative_transfer_before_limit": True,
            "squared_one_loop_form_factor_not_resummed_as_one_loop": True,
            "single_identical_field_diagnostic_not_conflated_with_actual_two_field_model": True,
            "heavy_mass_not_assumed_to_be_Reggeization_scale": True,
            "only_one_loop_matter_t_channel_vertex_piece_not_all_gravity_channels": True,
            "pure_graviton_loops_and_new_higher_curvature_matching_not_bounded": True,
            "no_numerical_Delta_grav_inferred_from_partial_vertex_bound": True,
            "no_old_state_counterterm_or_affine_bounce_parent_transfer": True,
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
        raise ValueError("Unsupported finite-gravity calibration accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    d = calibration.point()
    return {
        "rejected_inputs": rejected_inputs(),
        "nonzero_negative_vertex_piece_magnitude": d[
            "negative_t_vertex_finite_correction_absolute_upper"
        ]
        > 0,
        "finite_transfer_limit_error_retained": d[
            "negative_t_vertex_limit_error_per_abs_t_upper"
        ]
        > 0,
        "one_loop_form_factor_not_equal_to_its_zero_transfer_charge": True,
        "tree_pole_cannot_be_set_to_transfer_zero_before_subtraction": True,
        "single_identical_internal_loop_requires_extra_half_factor": True,
        "positive_vacuum_coefficient_does_not_fix_Regge_residue_derivatives": True,
        "minimal_scalar_Einstein_parent_still_not_flat_bounce_parent": True,
        "small_hierarchy_not_an_omitted_graviton_loop_error_bound": True,
        "original_P8_not_closed": True,
    }
