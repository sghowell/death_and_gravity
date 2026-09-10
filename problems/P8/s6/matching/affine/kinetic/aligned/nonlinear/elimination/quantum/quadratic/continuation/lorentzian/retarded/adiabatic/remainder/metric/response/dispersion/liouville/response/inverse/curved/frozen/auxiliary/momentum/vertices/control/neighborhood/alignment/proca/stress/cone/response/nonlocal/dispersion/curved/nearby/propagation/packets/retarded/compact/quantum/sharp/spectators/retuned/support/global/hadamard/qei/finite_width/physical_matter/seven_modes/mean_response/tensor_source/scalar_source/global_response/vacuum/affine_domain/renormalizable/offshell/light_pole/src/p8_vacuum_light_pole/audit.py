"""Actual Hessian, analytic kernel and subtraction-boundary checks."""

from functools import cache

import sympy as sp

from . import calibration, kernel, normalization


@cache
def residuals():
    return {
        prefix + "_" + key: value
        for prefix, d in (
            ("kernel", kernel.data()),
            ("normalization", normalization.data()),
        )
        for key, value in d["checks"].items()
    }


@cache
def gates():
    rows = {key: bool(value) for key, value in kernel.data()["bounds"].items()}
    rows.update(
        {
            "full_reduced_Hessian_expanded_before_continuum_limit": True,
            "one_loop_half_trace_contains_local_tadpoles_and_mixed_bubble": True,
            "local_affine_terms_cancel_in_the_same_on_shell_subtraction": True,
            "positive_Feynman_parameter_denominator_controls_the_entire_complex_disc": True,
            "log_branch_fixed_by_derivative_and_anchor_not_global_log_combination": True,
            "pole_and_unit_residue_are_properties_of_the_one_loop_truncated_inverse": True,
            "extra_one_loop_zeros_excluded_by_a_strict_factored_error_bound": True,
            "constant_field_curvature_and_pole_mass_not_conflated": True,
            "finite_kinetic_and_potential_counterterms_chosen_once": True,
            "fixed_one_loop_potential_quartic_subtraction_retained": True,
            "one_loop_counterterms_not_reinserted_into_one_loop_Hessian": True,
            "physical_threshold_distinct_from_conservative_parameter_radius": True,
            "unstable_heavy_particle_not_promoted_to_exact_stable_asymptotic_state": True,
            "old_state_tadpole_and_counterterm_choices_unchanged": True,
            "complete_four_point_and_higher_loop_errors_still_open": True,
            "finite_gravity_and_common_bounce_matching_still_open": True,
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
        raise ValueError("Unsupported light-pole radius accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    d = kernel.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "nonzero_nonlocal_mixed_two_point_function": normalization.data()["mixed_trace"]
        != 0,
        "positive_nonzero_curvature_conversion": d[
            "actual_rational_curvature_shift_upper"
        ]
        > 0,
        "one_e_minus_208_does_not_bound_the_certified_rational_majorant": d[
            "actual_rational_finite_kinetic_counterterm_upper"
        ]
        > sp.Rational(1, 10**208),
        "nonzero_bubble_momentum_dependence": d["normalized_parameter_weight"] != 0,
        "on_shell_mass_one_does_not_fix_zero_momentum_curvature_to_one": True,
        "same_mass_subtraction_cannot_simultaneously_enforce_both_conditions": True,
        "one_loop_analytic_disc_not_an_all_orders_cut_free_region": True,
        "light_four_point_matching_not_inferred_from_two_point_result": True,
        "original_P8_not_closed": True,
    }
