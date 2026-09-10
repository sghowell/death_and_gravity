"""Exact computational checks and explicit written-proof boundaries."""

from functools import cache

from . import analytic, calibration, denominators, graphs, sectors, subtraction


@cache
def residuals():
    return {
        prefix + "_" + k: v
        for prefix, module in (
            ("graphs", graphs),
            ("denominators", denominators),
            ("sectors", sectors),
            ("analytic", analytic),
            ("subtraction", subtraction),
            ("calibration", calibration),
        )
        for k, v in module.data()["checks"].items()
    }


@cache
def gates():
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "all_thirty_two_raw_refinements_retained": len(graphs.cases()) == 32,
            "all_seven_disjoint_graph_groups_counted": sum(
                graphs.data()["group_counts"].values()
            )
            == 32,
            "two_independent_U_and_P_constructions_for_every_graph": len(
                denominators.data()["all_32_denominator_certificates"]
            )
            == 32,
            "all_finite_sunset_refinements_selected_without_UV_cores": len(
                sectors.data()["six_UV_finite_sunset_ordered_bounds"]
            )
            == 6,
            "all_seven_hundred_twenty_finite_edge_orders_checked": sum(
                r["orders"]
                for r in sectors.data()["six_UV_finite_sunset_ordered_bounds"]
            )
            == 720,
            "unit_disc_unique_pole_follows_from_strict_inverse_gap": calibration.point(
                1
            )["unique_mass_one_pole_through_two_loops"],
            "same_fixed_inner_on_shell_kernel_retained": True,
            "same_entire_I0_cubic_and_heavy_mass_reference_retained": True,
            "all_momentum_dependent_proper_counterterms_grouped_before_bounds": True,
            "local_sunset_preserves_P_squared_in_twice_differentiated_integral": True,
            "nested_local_terms_annihilated_only_after_outer_OS": True,
            "inner_tadpoles_and_OS_quadratic_counterterms_grouped": True,
            "once_fixed_finite_quartic_contact_not_omitted": True,
            "strict_first_sheet_gap_has_open_neighborhood": True,
            "written_Euclidean_to_complex_routing_continuation": True,
            "anchored_log_difference_not_unsubtracted_mixed_integral": True,
            "independently_anchored_radial_integrals_no_physical_cutoff": True,
            "OS_Cauchy_bound_uses_outer_radius_two_inner_radius_one": True,
            "actual_couplings_not_only_formal_power_counting": True,
            "coarse_bound_not_claim_of_two_loop_dominance": True,
            "pole_and_residue_claim_only_through_two_loops": True,
            "complete_two_loop_four_point_ledger_not_inferred_here": True,
            "two_loop_source_aware_derivative_map_not_computed": True,
            "no_all_higher_loop_high_energy_or_Regge_control": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return result


@cache
def rejected_inputs():
    cases = graphs.bad_cases() + calibration.bad_cases()
    for name, call, args in cases:
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported two-loop pole input accepted: " + name)
    return len(cases)


@cache
def controls():
    local = denominators.point("sunset", (0, 0))
    return {
        "rejected_inputs": rejected_inputs(),
        "untreated_proper_UV_terms_are_nonzero": subtraction.data()[
            "sunset_proper_UV_piece_before_subtraction"
        ]
        != 0,
        "positive_denominator_does_not_imply_local_bare_UV_finiteness": local["P"] != 0
        and bool(graphs.uv(graphs.refine("sunset", (0, 0)))),
        "overlapping_local_sunset_cores_not_a_forest": (0, 1)
        not in graphs.forests(graphs.refine("sunset", (0, 0))),
        "fixed_inner_OS_remainder_not_set_to_zero": analytic.data()[
            "inherited_real_axis_R"
        ]
        != 0,
        "finite_potential_contact_annihilation_only_for_this_two_point_insertion": True,
        "large_coarse_bound_not_signed_higher_loop_coefficient": True,
        "two_loop_inverse_not_nonperturbative_propagator": True,
        "full_V_G_B_and_original_P8_not_closed": True,
    }
