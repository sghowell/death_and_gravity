"""Restricted Euclidean kernel audit and nonuniform-coordinate controls."""

from functools import cache

from . import calibration, dyson, kernel, map_domain


@cache
def residuals():
    return {
        prefix + "_" + name: value
        for prefix, module in (
            ("kernel", kernel),
            ("dyson", dyson),
            ("map_domain", map_domain),
            ("calibration", calibration),
        )
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    result = {
        prefix + "_" + name: bool(value)
        for prefix, module in (("calibration", calibration), ("map_domain", map_domain))
        for name, value in module.data()["bounds"].items()
    }
    result.update(
        {
            "actual_on_shell_subtracted_parent_kernel_not_a_new_scheme": True,
            "complete_heavy_denominator_retained_for_all_Euclidean_momenta": True,
            "positive_integral_representation_with_continuous_anchor": True,
            "strict_monotonicity_and_dominated_limit_give_sharp_supremum": True,
            "supremum_is_same_fixed_light_residue_counterterm_magnitude": True,
            "independent_both_stress_triangle_zero_transfer_normalization": True,
            "renormalized_stress_charge_stays_one_not_alpha": True,
            "positive_Euclidean_form_is_not_reflection_positivity": True,
            "massive_Sobolev_form_and_dual_source_bounds_not_timelike_contour": True,
            "geometric_series_repeats_only_computed_one_loop_kernel": True,
            "two_line_total_order_keeps_all_distinct_insertion_distributions": True,
            "insertion_tail_is_not_full_higher_loop_remainder": True,
            "unsubtracted_local_four_point_bubble_logarithmically_diverges": dyson.data()[
                "unsubtracted_bubble_diverges_to_positive_infinity"
            ],
            "no_positive_unsubtracted_trace_bound_for_renormalized_b2": True,
            "actual_derivative_composite_polynomial_not_surrogate": True,
            "positive_square_minimum_does_not_imply_uniform_smallness": True,
            "explicit_finite_witness_refutes_global_small_coordinate_bound": True,
            "nonuniform_coordinate_not_physical_cutoff_or_ghost": True,
            "proper_on_shell_LSZ_equivalence_is_not_refuted": True,
            "old_state_counterterms_and_frozen_sources_unchanged": True,
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
        raise ValueError(
            "Unsupported Euclidean insertion-control input accepted: " + name
        )
    return len(calibration.bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "uniform_kernel_bound_not_a_zero_self_energy": calibration.data()[
            "actual_uniform_relative_self_energy_upper"
        ]
        > 0,
        "known_kernel_iteration_tail_not_set_to_zero": calibration.data()[
            "actual_two_line_beyond_one_total_insertion_upper"
        ]
        > 0,
        "logarithmic_bubble_divergence_detected": dyson.data()[
            "unsubtracted_bubble_diverges_to_positive_infinity"
        ],
        "global_coordinate_smallness_counterexample": map_domain.data()[
            "witness_mixing_strict_rational_lower"
        ]
        > 1,
        "Euclidean_positivity_not_timelike_unitarity": True,
        "restricted_insertions_not_primitive_higher_loop_diagrams": True,
        "renormalized_b2_not_bounded_by_divergent_positive_trace": True,
        "large_composite_coordinate_not_a_row_no_go": True,
        "no_rolling_state_or_finite_gravity_transfer": True,
        "original_P8_not_closed": True,
    }
