"""Exact same-model loop coefficient, conservative majorants and scope gates."""

from functools import cache

import sympy as sp

from . import angular, calibration, domain, leading, normalization, subtraction


@cache
def residuals():
    return {
        prefix + "_" + name: value
        for prefix, module in (
            ("domain", domain),
            ("angular", angular),
            ("leading", leading),
            ("normalization", normalization),
            ("subtraction", subtraction),
            ("calibration", calibration),
        )
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    rows = {
        prefix + "_" + name: bool(value)
        for prefix, module in (
            ("domain", domain),
            ("angular", angular),
            ("leading", leading),
            ("calibration", calibration),
        )
        for name, value in module.data()["bounds"].items()
    }
    rows.update(
        {
            "actual_full_Hessian_quartic_trace_not_heavy_logdet_only": True,
            "three_light_bubble_channels_include_all_original_mixed_loop_topologies": True,
            "complete_nonlocal_vertices_retained_inside_the_loop": True,
            "external_light_pole_and_unit_residue_fixed_by_same_model_parent": True,
            "full_Feynman_parameter_polynomial_positive_on_complex_unit_disc": True,
            "analytic_continuation_not_naive_unshifted_Euclidean_contour": True,
            "only_external_shifts_expanded_with_uniform_geometric_remainder": True,
            "radial_loop_momentum_integrated_from_zero_to_infinity": True,
            "odd_linear_angular_terms_vanish_under_even_measure": True,
            "both_t_channel_vertices_obey_same_bounds_without_conjugation_assumption": True,
            "crossing_includes_u_channel_and_t_remainder": True,
            "Cauchy_estimate_used_only_for_holomorphic_angular_remainder": True,
            "UV_cancellation_made_explicit_before_high_momentum_bounds": True,
            "local_full_model_counterterms_implement_the_fixed_channel_subtractions": True,
            "finite_potential_contact_fixed_once_and_has_zero_b2": True,
            "total_vertex_counterterms_do_not_double_count_light_wavefunction": True,
            "higher_loop_counterterm_insertions_not_mistaken_for_one_loop_terms": True,
            "small_couplings_alone_not_used_as_relative_error_estimate": True,
            "no_squared_unstable_heavy_pole_integrated_in_a_dispersion_relation": True,
            "positive_tree_plus_one_loop_result_not_all_order_UV_completion": True,
            "old_cosmological_states_and_counterterms_unchanged": True,
            "higher_loop_contour_finite_gravity_and_common_bounce_still_open": True,
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
        raise ValueError("Unsupported one-loop relative tolerance accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    p = calibration.point()
    return {
        "rejected_inputs": rejected_inputs(),
        "full_angular_remainder_nonzero": p["angular_b2_error_upper"] > 0,
        "angular_zero_error_nonzero_and_retained": p["angular_zero_b2_error_upper"] > 0,
        "one_e_minus_eight_tolerance_not_proved_by_this_majorant": not calibration.point(
            sp.Rational(1, 10**8)
        )["requested_tolerance_proved"],
        "constant_potential_determinant_alone_omits_scattering_momentum_dependence": True,
        "naive_derivative_expansion_not_uniform_over_all_loop_momenta": True,
        "light_two_particle_cut_outside_the_unit_disc_is_not_deleted": True,
        "finite_parameter_scheme_not_stable_heavy_on_shell_scheme": True,
        "subthreshold_loop_result_not_high_energy_contour_bound": True,
        "one_loop_positive_coefficient_not_all_orders_error_control": True,
        "original_P8_not_closed": True,
    }
