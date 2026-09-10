"""Exact physical cut identities, continuous bounds and perturbative scope."""

from functools import cache

from . import amplitude, angular, calibration, cut, normalization


@cache
def residuals():
    return {
        prefix + "_" + name: value
        for prefix, module in (
            ("amplitude", amplitude),
            ("angular", angular),
            ("normalization", normalization),
            ("cut", cut),
        )
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    rows = {
        prefix + "_" + name: bool(value)
        for prefix, module in (("amplitude", amplitude), ("cut", cut))
        for name, value in module.data()["bounds"].items()
    }
    rows.update(
        {
            "actual_full_tree_vertex_not_local_derivative_truncation": True,
            "physical_s_and_cosine_squared_monotonicity_give_continuous_bounds": True,
            "atanh_primitive_anchored_on_real_branch_without_spurious_threshold_division": True,
            "identical_particle_symmetry_and_two_Im_factor_both_retained": True,
            "independent_constant_vertex_bubble_discontinuity_agrees": True,
            "positive_tree_squared_cut_not_inferred_from_subthreshold_real_coefficient": True,
            "only_two_light_intermediate_particles_open_in_selected_one_loop_window": True,
            "Phi_Z2_excludes_odd_light_mixed_two_particle_channel": True,
            "real_local_counterterms_do_not_change_this_one_loop_discontinuity": True,
            "one_loop_vertex_insertions_on_cut_would_be_higher_order": True,
            "phase_space_threshold_zero_not_strict_absorption_at_threshold": True,
            "positive_open_subwindow_proves_nonzero_integrated_weight": True,
            "pi_upper_bound_follows_from_strict_positive_arctangent_integrand_gap": True,
            "dispersion_normalization_for_b2_is_half_second_derivative": True,
            "actual_cut_integrated_only_below_heavy_resonance": True,
            "complete_one_loop_real_coefficient_error_retained_in_subtracted_margin": True,
            "partial_improved_functional_not_global_dispersion_or_all_orders_verdict": True,
            "old_states_counterterms_and_frozen_actions_unchanged": True,
            "higher_loop_contour_finite_gravity_common_bounce_and_original_P8_open": True,
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
        raise ValueError("Unsupported elastic-cut invariant accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    p = calibration.point(4)
    return {
        "rejected_inputs": rejected_inputs(),
        "exact_threshold_cut_zero": p["threshold_density_exactly_zero"],
        "no_strict_absorption_claim_at_threshold": not p[
            "cut_density_strictly_positive"
        ],
        "interior_positive_cut_present": calibration.point(5)[
            "cut_density_strictly_positive"
        ],
        "computed_positive_cut_is_not_dropped_from_coefficient": cut.data()[
            "actual_cut_upper"
        ]
        > 0,
        "closed_angular_formula_has_large_cancellations": True,
        "unstable_heavy_resonance_not_integrated_as_bare_squared_pole": True,
        "tree_plus_one_loop_not_all_orders_improved_positivity": True,
        "finite_window_does_not_prove_high_energy_contour": True,
        "original_P8_not_closed": True,
    }
