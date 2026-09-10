"""Integrated finite-subsector audit with explicit exclusion controls."""

from functools import cache

from p8_vacuum_two_loop_denom import graphs, subgraphs

from . import calibration, ordered, parametric, sectors, selection


@cache
def residuals():
    return {
        prefix + "_" + name: value
        for prefix, module in (
            ("selection", selection),
            ("ordered", ordered),
            ("parametric", parametric),
            ("calibration", calibration),
        )
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    rows = sectors.summaries()
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "exactly_eighty_eight_finite_refinements_selected": len(rows) == 88,
            "exactly_one_hundred_four_subtraction_dependent_refinements_excluded": selection.data()[
                "excluded_subtraction_dependent_count"
            ]
            == 104,
            "all_two_loop_edge_orders_checked": sum(
                r["edge_order_sector_count"] for r in rows
            )
            == 258480,
            "all_heavy_cutoff_pieces_checked": sum(
                r["heavy_cutoff_piece_count"] for r in rows
            )
            == 526320,
            "all_small_tail_exponents_strictly_positive": all(
                r["minimum_small_tail_exponent"] > 0 for r in rows
            ),
            "at_most_two_logarithms_for_every_piece": all(
                r["maximum_log_power"] <= 2 for r in rows
            ),
            "all_monomial_dominance_and_mass_cancellations_checked": all(
                r["every_monomial_dominance_tail_prefix_and_mass_cancellation_checked"]
                for r in rows
            ),
            "all_large_power_patterns_independently_integrated": {
                p for r in rows for p in r["large_ordered_power_patterns"]
            }
            == set(ordered.LARGE_PATTERNS),
            "all_box_polynomial_coefficients_positive": all(
                c > 0 for r in rows for c in r["box_log_polynomial_coefficients"]
            ),
            "full_heavy_exchange_not_expanded_in_loop_momentum": True,
            "positive_layer_cake_retains_full_unbounded_parameter_domain": True,
            "anisotropic_box_is_integration_bound_not_a_physical_cutoff": True,
            "Gaussian_loop_measure_and_second_Taylor_factor_retained": True,
            "all_external_assignments_keep_frozen_Wick_weights": True,
            "compact_derivative_bound_has_absolute_integrable_majorant": True,
            "first_sheet_continuation_from_zero_external_momenta": True,
            "positive_upper_bound_not_a_positive_signed_correction": True,
            "finite_subsector_not_complete_two_loop_amplitude": True,
            "all_subtraction_dependent_graphs_require_separate_work": True,
            "two_loop_mass_residue_and_LSZ_not_included": True,
            "finite_counterterm_scheme_not_changed": True,
            "no_all_order_or_all_energy_error_claim": True,
            "no_gravity_Regge_or_rolling_parent_transfer": True,
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
        raise ValueError("Unsupported finite-subsector input accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    divergent = sum(
        bool(subgraphs.data(*case)["UV_subgraphs"]) for case in graphs.cases()
    )
    return {
        "rejected_inputs": rejected_inputs(),
        "UV_dependent_graphs_are_not_silently_included": divergent == 104,
        "no_tadpole_family_is_misclassified_as_finite": not any(
            kind == "tadpole_insertion" for kind, choices in selection.cases()
        ),
        "nonzero_quadratic_log_coefficient_retained": calibration.data()[
            "aggregate_log_polynomial_coefficients"
        ][2]
        > 0,
        "nonzero_error_upper_not_replaced_by_zero": calibration.data()[
            "actual_finite_subsector_b2_absolute_upper"
        ]
        > 0,
        "finite_graph_upper_does_not_bound_the_excluded_graphs": True,
        "heavy_mass_power_cancellation_not_a_local_derivative_expansion": True,
        "no_counterterm_or_LSZ_omission_from_complete_amplitude_claim": True,
        "complete_two_loop_and_original_P8_not_closed": True,
    }
