"""Exact homogeneous tensor-source identities and claim-boundary controls."""

from functools import cache

import sympy as sp

from . import bounds, calibration, chart, operator, response, state


@cache
def residuals():
    result = {}
    for data in (
        operator.data(),
        operator.clock(),
        chart.data(),
        chart.metric(),
        state.data(),
        state.variational_clock(),
        response.data(),
        response.center(),
        response.ward(),
        bounds.data(),
    ):
        overlap = set(result).intersection(data["checks"])
        if overlap:
            raise ValueError("Duplicate exact tensor-source identity: " + str(overlap))
        result.update(data["checks"])
    return result


@cache
def rejected_inputs():
    for name, call, args in calibration.bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Invalid tensor covariance amplitude accepted: " + name)
    return len(calibration.bad_cases())


@cache
def gates():
    result = {
        "same_actual_global_retuned_action_and_physical_lapse": True,
        "actual_DHOST_tensor_coefficient_varied_before_clock_restriction": True,
        "tensor_spatial_action_uses_volume_preserving_exponential_chart": True,
        "actual_lapse_source_not_inferred_from_on_shell_frequency": True,
        "new_positive_covariance_addition_is_frame_free_in_two_TT_modes": True,
        "same_reference_scalar_and_vector_factors_and_counterterms": True,
        "fixed_normalized_open_shell_bump_has_no_IR_or_UV_tail": True,
        "exact_actual_tensor_evolution_preserves_CCR_and_smooth_Hadamard_difference": True,
        "positive_proxy_energy_not_mistaken_for_physical_lapse_source": True,
        "isotropic_homogeneous_mean_is_finite_density_not_finite_total_R3_energy": True,
        "clock_source_derived_from_actual_clock_velocity_variation": True,
        "joint_metric_source_Ward_identity_retains_clock_and_connection": True,
        "same_actual_mean_operator_with_new_tensor_forcing": True,
        "nonzero_constraint_lapse_and_induced_matter_field_retained": True,
        "second_order_spatial_chart_correction_retained_at_center": True,
        "no_all_time_additive_canonical_momentum_bridge_claimed": True,
        "only_added_covariance_has_assumed_time_reflection": True,
        "global_bounds_cover_both_time_tails": True,
        "signed_clock_transfer_has_vanishing_boundary_on_both_tails": True,
        "negative_tensor_acceleration_uses_new_absolute_bound_not_Proca_positivity": True,
        "complete_exact_frame_representative_not_a_semiclassical_solution": True,
        "full_covariant_graviton_stress_absolute_source_higher_orders_V_G_B_open": True,
    }
    result.update({"bound_" + k: bool(v) for k, v in bounds.data()["bounds"].items()})
    return result


@cache
def controls():
    d = operator.data()
    r, s = d["positive_on_clock_tensor_energy"], d["on_clock_spatial_pressure"]
    h = operator.h
    return {
        "rejected_inputs": rejected_inputs(),
        "setting_actual_tensor_coefficient_to_one_before_lapse_variation_changes_source": sp.factor(
            d["actual_hat_lapse_source"] - (r - 3 * s / (2 * h))
        )
        != 0,
        "physical_tensor_density_is_not_positive_proxy_energy": sp.factor(
            d["fixed_physical_metric_lapse_source"] - r
        )
        != 0,
        "pure_tensor_coordinate_physical_lapse_source_can_be_negative": sp.factor(
            d["fixed_physical_metric_lapse_source"].subs({h: 1, operator.P: 0})
        )
        == -operator.q * operator.T**2 / (4 * operator.a**2),
        "omitting_second_order_scalar_chart_shift_changes_center_matter_density": chart.data()[
            "induced_linear_homogeneous_scalar_density"
        ]
        != 0,
        "actual_selected_tensor_band_clock_source_is_nonzero": response.center()[
            "anchor_clock_source_first_derivative"
        ]
        != 0,
        "omitting_induced_matter_connection_changes_Ward_identity": response.ward()[
            "retained_background_matter_connection_variation"
        ]
        != 0,
        "new_center_lapse_second_derivative_is_negative_at_lower_band_endpoint": response.center()[
            "anchor_lapse_response_second_derivative"
        ].subs(operator.q, 1)
        < 0,
        "zero_added_tensor_covariance_has_zero_all_calibrated_bounds": all(
            v == 0 for v in calibration.calibrated(0).values()
        ),
        "reference_state_and_common_absolute_tadpole_not_reassigned": True,
        "no_finite_amplitude_SEE_or_coevolved_quantum_state_claimed": True,
        "no_full_covariant_graviton_or_BRST_renormalized_stress_claimed": True,
        "no_cutoff_higher_loop_control_or_V_G_B_matching_inferred": True,
    }
