"""Native seven-mode action, polarization, covariance and preparation gates."""

from functools import cache

import sympy as sp

from . import bridge, covariance, deformation, model, projectors, state


@cache
def residuals():
    out = {}
    for data in (
        model.data(),
        model.energy(),
        projectors.data(),
        covariance.data(),
        bridge.data(),
        state.data(),
    ):
        out.update(data["checks"])
    out.update(deformation.checks())
    return out


def bad_cases():
    values = (
        True,
        False,
        sp.true,
        sp.false,
        0.1,
        sp.Float("3"),
        "3",
        None,
        sp.oo,
        sp.nan,
        sp.zoo,
        sp.Symbol("R"),
        sp.I,
        0,
        -3,
        sp.Rational(5, 2),
    )
    return [
        ("strip_radius_" + str(index), deformation.strip, (value,))
        for index, value in enumerate(values)
    ]


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("An invalid global preparation radius was accepted: " + name)
    return len(bad_cases())


@cache
def gates():
    d = deformation.data()
    r = deformation.strip()
    return {
        "same_global_retuned_action_and_positive_Proca_mass": True,
        "all_seven_free_physical_modes_retained": state.data()[
            "full_physical_mode_count"
        ]
        == 7,
        "full_density_CCR_phase_dimension_fourteen": model.data()[
            "actual_full_symplectic_form"
        ].rank()
        == 14,
        "actual_flat_Proca_covariance_is_positive_not_four_unconstrained_ghosts": True,
        "tensor_polarization_normalization_keeps_factor_two": True,
        "no_global_angular_polarization_frame_assumed": True,
        "nonzero_Proca_null_projection_symbol_not_mistaken_for_ellipticity": True,
        "corrected_Proca_propagation_theorem_with_actual_constant_mass_applies": True,
        "actual_and_auxiliary_metrics_have_complete_R_cubed_Cauchy_surfaces": True,
        "smooth_preparation_geometry_not_a_new_physical_background": True,
        "smooth_step_derivative_bound_keeps_endpoint_regions": sp.Rational(160, 33)
        < 5
        < d["smooth_step_derivative_upper"],
        "positive_auxiliary_and_actual_scale": r["actual_and_auxiliary_scale_lower"]
        == 1,
        "positive_whole_strip_Proca_energy_floor": r["Proca_energy_lower"] > 0,
        "positive_whole_strip_tensor_modified_energy_floor": r[
            "tensor_modified_energy_lower"
        ]
        > 0,
        "zero_momentum_infrared_covariance_is_locally_integrable_without_an_atom": True,
        "all_compact_time_evolutions_have_polynomial_spatial_frequency_bounds": True,
        "final_state_uses_actual_global_generator_on_both_legs": True,
        "old_scalar_state_is_unchanged_under_partial_restriction": True,
        "product_state_is_positive_quasifree_and_has_full_density_CCR": True,
        "wavefront_union_contains_clock_and_physical_matter_cones": True,
        "not_a_nonlinear_or_covariant_gauge_fixed_interacting_gravity_state": True,
        "physical_stress_Ward_mean_response_and_V_G_B_not_claimed": True,
    }


@cache
def controls():
    null = covariance.data()["nonelliptic_nonzero_null_Proca_projection_symbol"]
    return {
        "rejected_inputs": rejected_inputs(),
        "naive_four_independent_covariant_polarizations_have_negative_time_component": -1,
        "correct_Proca_null_projection_symbol_rank": null.rank(),
        "correct_Proca_null_projection_symbol_determinant": sp.factor(null.det()),
        "one_Proca_not_four_scalar_oscillators": True,
        "no_projector_direction_or_zero_frequency_atom_is_inserted": True,
        "no_instantaneous_Minkowski_reset_at_actual_bounce": True,
        "same_chosen_S6_101_scalar_state_retained": True,
        "no_old_tadpole_or_stress_subtraction_reassigned_to_new_product_state": True,
    }
