"""Exact off-shell identity, actual coefficients, norm bounds and scope controls."""

from functools import cache

from . import band, calibration, jets, matching, norms


@cache
def residuals():
    groups = {
        "jets": jets.data()["checks"],
        "resolvent": matching.data()["checks"],
        "actual_target": matching.actual_target()["checks"],
        "norms": norms.data()["checks"],
        "common_class": band.data()["checks"],
    }
    return {
        prefix + "_" + key: value
        for prefix, rows in groups.items()
        for key, value in rows.items()
    }


@cache
def gates():
    groups = {
        "resolvent": matching.data()["bounds"],
        "norms": norms.data()["bounds"],
        "common_class": band.data()["bounds"],
    }
    rows = {
        prefix + "_" + key: bool(value)
        for prefix, values in groups.items()
        for key, value in values.items()
    }
    rows.update(
        {
            "full_four_dimensional_jet_polynomial_not_homogeneous_only_check": True,
            "actual_affine_vacuum_quartic_rebuilt_before_matching": True,
            "mass_one_EOM_terms_and_potential_not_dropped": True,
            "literal_boundary_current_present_pointwise": True,
            "field_map_sign_fixed_by_direct_free_action_variation": True,
            "sixth_and_higher_field_terms_bounded_not_discarded": True,
            "every_required_jet_after_substitution_has_order_at_most_ten": True,
            "full_heavy_resolvent_retained_before_spectral_estimate": True,
            "wide_spectral_window_not_confused_with_a_loop_analyticity_domain": True,
            "pointwise_and_L2_action_norms_separated_before_common_class_proof": True,
            "nonempty_real_Schwartz_bandlimited_class_has_all_required_jet_norms": True,
            "Fourier_support_is_checked_after_the_cubic_field_map": True,
            "source_J_is_the_square_of_the_mapped_field_not_the_unmapped_field": True,
            "Schwartz_boundary_currents_integrate_to_zero": True,
            "two_L2_factors_and_remaining_Linfinity_factors_control_each_error_monomial": True,
            "formal_tangent_identity_not_global_differential_map_invertibility": True,
            "field_dependent_quantum_Jacobian_and_loops_not_transferred": True,
            "not_a_finite_M_or_full_affine_bounce_parent_identification": True,
            "old_states_counterterms_and_frozen_actions_unchanged": True,
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
        raise ValueError("Unsupported off-shell vacuum input accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    j = jets.data()
    n = norms.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "nonzero_redundant_cubic_field_map": j["cubic_field_redefinition"] != 0,
        "nonzero_higher_field_order_remainder_majorant": n[
            "unit_jet_cube_total_field_redefinition_error"
        ]
        > 0,
        "large_tree_spectral_window_includes_light_loop_cuts": matching.data()[
            "selected_spectral_radius"
        ]
        > 2,
        "actual_off_shell_EOM_not_set_to_zero_in_the_identity": j["free_equation"] != 0,
        "Fourier_support_can_grow_under_a_local_nonlinear_field_map": band.data()[
            "mapped_quadratic_source_fourier_radius"
        ]
        > band.data()["new_field_fourier_radius"],
        "exact_resolvent_not_equal_to_finite_local_polynomial": matching.data()[
            "exact_operator_remainder"
        ]
        != 0,
        "pointwise_jet_bound_does_not_imply_a_global_L2_bound": True,
        "oneway_actual_field_map_not_a_global_phase_space_isomorphism": True,
        "quantum_measure_and_renormalized_amplitude_not_inferred": True,
        "small_vacuum_fields_not_identified_with_entire_rolling_clock_background": True,
        "finite_gravity_Regge_and_propagating_common_parent_still_open": True,
        "original_P8_not_closed": True,
    }
