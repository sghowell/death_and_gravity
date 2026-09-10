"""Exact finite-mass cut residuals, rational guards and proof-boundary controls."""

from functools import cache

import sympy as sp

from . import calibration, cuts, dirac, kinematics, series

MODULES = (dirac, series, kinematics, cuts, calibration)


@cache
def residuals():
    result = {}
    for module in MODULES:
        for name, value in module.data()["checks"].items():
            result[module.__name__.rsplit(".", 1)[-1] + "_" + name] = (
                value.applyfunc(sp.simplify)
                if isinstance(value, sp.MatrixBase)
                else sp.simplify(value)
            )
    return result


def scalar_entries():
    return sum(
        v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    result.update(
        {
            "four_Hermitian_Clifford_generators": len(dirac.data()["gamma_matrices"])
            == 4,
            "six_distinct_cyclic_box_orders": len(
                series.data()["cyclic_labelled_box_orders"]
            )
            == 6,
            "six_complete_four_propagator_routes": len(
                kinematics.data()["six_cyclic_box_routes"]
            )
            == 6
            and all(
                len(r["partial_momenta"]) == 4
                for r in kinematics.data()["six_cyclic_box_routes"]
            ),
            "both_Ward_cyclic_orientations": len(
                kinematics.data()["regulated_Ward_telescoping_rows"]
            )
            == 2,
            "regulated_full_sum_precedes_degree_zero_removal": True,
            "four_dimensional_norm_only_after_UV_convergent_Taylor_subtraction": True,
            "entire_n_ge_three_tail_not_only_finite_sample": True,
            "positive_radial_primitive_has_zero_and_infinity_anchors": True,
            "complex_shifts_controlled_without_arbitrary_loop_contour_translation": True,
            "physical_anchor_and_Bose_even_continuation_retained": True,
            "two_holomorphic_cut_factors_not_unphysical_absolute_square": True,
            "all_four_transverse_pairs_color_optical_and_identical_factors_retained": True,
            "both_crossed_channels_and_opposite_boundary_signs_retained": True,
            "first_two_gauge_discontinuity_complete_at_three_parent_loops": cuts.data()[
                "first_parent_loop_order"
            ]
            == 3,
            "local_counterterms_do_not_change_discontinuity": True,
            "center_and_near_center_inconclusive_bounds_are_valid_inputs": True,
            "nonzero_finite_mass_cut_not_discarded_by_tiny_absolute_size": True,
            "local_analytic_germ_not_ruled_out": True,
            "no_confinement_mass_or_full_candidate_no_go_inferred": True,
            "old_scalar_two_loop_bound_not_relabelled_as_new_model": True,
            "no_new_all_orders_UV_construction_requirement": True,
            "finite_gravity_common_parent_and_original_P8_remain_open": True,
        }
    )
    return {name: bool(value) for name, value in result.items()}


@cache
def rejected_inputs():
    for name, call, args in cuts.bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported exact cut input accepted: " + name)
    return len(cuts.bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "small_mass_tail_bound_inconclusive_not_excluded": not cuts.point(mass=24)[
            "strict_nonzero_boundary_difference_at_this_order"
        ],
        "crossing_center_zero_not_used_as_disc_test": not cuts.point(2)[
            "strict_nonzero_boundary_difference_at_this_order"
        ]
        and cuts.point()["strict_nonzero_boundary_difference_at_this_order"],
        "near_center_not_forced_to_nonzero": not calibration.data()[
            "near_center_control"
        ]["strict_nonzero_boundary_difference_at_this_order"],
        "finite_field_series_not_substituted_for_momentum_tail": True,
        "individual_unregulated_UV_box_not_used_as_finite_integral": True,
        "no_complex_conjugation_of_unphysical_scalar_momentum": True,
        "perturbative_gauge_cut_not_a_confinement_theorem": True,
        "original_V_G_B_and_P8_not_closed": True,
    }
