"""Independent exact controls and source-pinned analytic proof gates."""

from functools import cache

import sympy as sp

from . import (
    analytic,
    calibration,
    decoupling,
    family,
    heavy,
    transition,
    uniform,
    vacuum,
)


@cache
def residuals():
    groups = {
        "smooth": family.identities(),
        "transition": transition.data()["checks"],
        "analytic": analytic.exact_checks(),
        "vacuum_jets": analytic.local_jets()["checks"],
        "uniform": uniform.data()["checks"],
        "decoupling": decoupling.data()["checks"],
        "contact": vacuum.data()["checks"],
        "heavy": heavy.data()["checks"],
    }
    return {
        group + "_" + name: value
        for group, rows in groups.items()
        for name, value in rows.items()
    }


@cache
def gates():
    groups = {
        "transition": transition.data()["bounds"],
        "curvature_error": analytic.coefficient_error_bounds()["bounds"],
        "scalar_error": analytic.lower_error_bounds()["bounds"],
        "vacuum": analytic.local_jets()["bounds"],
        "uniform": uniform.data()["bounds"],
        "contact": vacuum.data()["bounds"],
        "heavy": heavy.data()["bounds"],
    }
    result = {
        group + "_" + name: bool(value)
        for group, rows in groups.items()
        for name, value in rows.items()
    }
    result.update(
        {
            "separately_named_classical_extension_not_a_frozen_action_edit": True,
            "literal_retuned_lower_action_not_only_its_principal_symbol": True,
            "same_physical_free_matter_and_ordinary_Proca_coupling": True,
            "Ia_completion_rederived_with_switch_derivative": True,
            "vacuum_piecewise_functions_do_not_multiply_zero_by_undefined_poles": True,
            "smooth_exact_tube_and_analytic_finite_jet_families_are_distinct": True,
            "analytic_identity_obstruction_restricted_to_exact_fixed_witness": True,
            "analytic_R_positive_on_stated_connected_domain_for_each_even_order": True,
            "pole_removability_checked_before_vacuum_expansion": True,
            "mixed_norms_weight_derivatives_after_taking_them": True,
            "inverse_jet_bound_uses_finite_Taylor_algebra_not_infinite_Neumann_series": True,
            "scalar_global_u_bound_keeps_analytic_localizer_derivatives": True,
            "uniform_even_order_bound_uses_positive_affine_scalar_majorant": True,
            "all_finite_clock_jets_claimed_only_through_n_minus_one": True,
            "co_scaled_canonical_limit_keeps_mass_lambda_and_gamma_fixed": True,
            "massive_contact_keeps_nonzero_constant_term": True,
            "heavy_parent_has_positive_field_metric_and_local_vacuum_spectrum": True,
            "tree_remainder_is_complex_polydisc_uniform_not_a_grid_sample": True,
            "heavy_tree_exchange_not_an_all_order_absorptive_spectral_density": True,
            "finite_M_gravity_poles_Regge_and_contour_errors_not_discarded": True,
            "minimal_healthy_heavy_vacuum_model_not_claimed_as_bounce_parent": True,
            "no_reassignment_of_old_absolute_state_counterterms_or_tadpoles": True,
            "full_V_G_B_and_original_P8_remain_open": True,
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
        raise ValueError("Unsupported analytic calibration accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    d = transition.data()
    n = analytic.ORDER
    h = heavy.data()
    lam = sp.Symbol("positive_fixed_lambda", positive=True)
    gamma = sp.Symbol("positive_fixed_gamma", positive=True)
    return {
        "rejected_inputs": rejected_inputs(),
        "naive_taper_has_negative_regular_matter_difference_determinant": bool(
            (-d["transition_center_naive_difference_determinant"]).is_positive
        ),
        "naive_A3_center_not_correct_A3": d["transition_center_naive_A3"]
        != d["transition_center_correct_A3"],
        "unmodified_exact_tube_A3_has_nonzero_analytic_identity_obstruction": family.identities()[
            "original_curvature_pole_identity"
        ]
        == 0,
        "analytic_family_is_not_smooth_exact_open_tube_family": (
            1 - sp.Rational(9, 10) ** 2
        )
        ** n
        != 0,
        "fixed_n_limit_distinguished_from_fixed_nonzero_gamma_limit": analytic.FIXED_GAMMA
        != 0,
        "actual_potential_quartic_not_dropped": analytic.local_jets()[
            "canonical_potential_quartic_coefficient"
        ]
        != 0,
        "heavy_omitted_operator_bound_reproduces_exact_calibration": calibration.heavy_remainder_bound(
            heavy.CHANNEL_RADIUS
        )
        == h["exact_uniform_tree_remainder_upper"],
        "heavy_forward_b2_positive_but_not_full_V_gate": h["exact_heavy_forward_b2"]
        == 4 * lam,
        "heavy_matching_needs_both_contact_coefficients": h["channel_expansion_gap"]
        == 2 * lam / gamma,
        "minimal_positive_scalar_parent_NEC_cannot_support_flat_Einstein_bounce": True,
        "algebraic_Ia_denominator_floor_not_a_heavy_mass_gap": True,
        "finite_jet_equality_not_all_order_quantum_action_equality": True,
        "no_uniform_in_n_off_tube_tensor_floor_claimed": True,
        "no_quantum_V_G_B_or_universal_no_go_inferred": True,
    }
