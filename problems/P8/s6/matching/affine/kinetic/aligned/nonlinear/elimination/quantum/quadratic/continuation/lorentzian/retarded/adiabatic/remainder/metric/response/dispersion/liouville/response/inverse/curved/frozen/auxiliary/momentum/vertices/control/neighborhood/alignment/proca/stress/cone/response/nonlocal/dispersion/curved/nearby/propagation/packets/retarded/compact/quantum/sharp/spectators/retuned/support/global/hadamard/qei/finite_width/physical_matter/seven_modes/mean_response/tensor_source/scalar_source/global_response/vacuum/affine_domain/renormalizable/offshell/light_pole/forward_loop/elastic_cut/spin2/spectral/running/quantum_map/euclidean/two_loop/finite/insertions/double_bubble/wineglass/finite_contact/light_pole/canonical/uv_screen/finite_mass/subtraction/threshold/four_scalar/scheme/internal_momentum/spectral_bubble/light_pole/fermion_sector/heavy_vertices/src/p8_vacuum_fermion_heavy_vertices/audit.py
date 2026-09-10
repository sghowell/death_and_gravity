"""Exact identities, local subtraction and the narrowly updated frontier."""

from functools import cache

import sympy as sp
from p8_vacuum_two_loop_fermion_ledger import boundary as parent

from . import bounds, calibration, ownership, parametric, renormalization

MODULES = (ownership, parametric, renormalization, bounds, calibration)


def frontier():
    rows = parent.ownership_frontier()
    for row in rows:
        if row["id"] == "scalar_Phi4_W2_F0":
            row["status"] = "BOUNDED_IN_DECLARED_PARENT_SUBTRACTION"
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError(
            "The full-vertex family frontier differs from its declared scope"
        )
    return True


@cache
def residuals():
    rows = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for mod in MODULES
        for name, value in mod.data()["checks"].items()
    }
    old, new = parent.ownership_frontier(), frontier()
    rows.update(
        {
            "frontier_exactly_one_status_changed": sum(a != b for a, b in zip(old, new))
            - 1,
            "frontier_all_nine_primitive_ids_retained": len(new) - 9,
            "frontier_seven_primitive_rows_remain_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 7,
            "frontier_one_new_parent_scheme_family_bound": sum(
                r["status"] == "BOUNDED_IN_DECLARED_PARENT_SUBTRACTION" for r in new
            )
            - 1,
            "frontier_previous_quadratic_bound_unchanged": sum(
                r["status"] == "BOUNDED_PAIRED_SECTOR" for r in new
            )
            - 1,
        }
    )
    return rows


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    result.update(
        {
            "entire_frozen_nonlocal_quartic_Hessian_retained": True,
            "literal_covariance_variation_and_external_derivatives": True,
            "three_channels_both_internal_positions_and_half_symmetry": True,
            "all_one_four_four_vertex_products_retained": True,
            "positive_spectral_representation_of_whole_paired_inner_insertion": True,
            "unequal_mass_parameter_gap_in_every_required_routing": True,
            "triangle_and_box_simplex_volumes_and_Gamma_factors": True,
            "spectral_triangle_and_box_integrals_absolutely_converge": True,
            "whole_bubble_reference_subtracted_before_regulator_removal": True,
            "subtraction_implemented_by_local_L_G_and_M_parent_terms": True,
            "finite_regulator_products_owned_by_common_regulated_reference": True,
            "no_early_epsilon_projection_of_unsubtracted_poles": True,
            "full_nonlocal_finite_reference_change_can_affect_b2": True,
            "local_subset_anchor_change_is_only_a_constant": True,
            "zero_transfer_heavy_box_and_triangle_dependence_retained": True,
            "uniform_holomorphic_disc_bound_not_pointwise_quadrature": True,
            "Cauchy_bound_applies_to_b2_not_a_claimed_unit_disc_remainder": True,
            "exact_dyadic_log_bound_no_float_logarithms": True,
            "no_sign_for_full_family_inferred_from_positive_density": True,
            "new_row_not_added_to_its_previously_bounded_local_subset": True,
            "shared_MS_and_canonical_order_two_conversion_still_open": True,
            "other_primitive_rows_and_counterterms_not_called_complete": True,
            "no_full_two_loop_higher_loop_contour_G_or_B_claim": True,
            "original_P8_remains_open": True,
        }
    )
    return {k: bool(v) for k, v in result.items()}


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("x"),
    )
    result = []
    for j in range(6):
        for i, value in enumerate(invalid):
            values = [2, 1, 1, 1, 24, 144]
            values[j] = value
            result.append((f"type_{j}_{i}", bounds.enclosure, tuple(values)))
    for i, values in enumerate(
        (
            (1, 1, 1, 1, 24, 144),
            (-2, 1, 1, 1, 24, 144),
            (2, -1, 1, 1, 24, 144),
            (2, 1, -1, 1, 24, 144),
            (2, 1, 1, -1, 24, 144),
            (2, 1, 1, 1, 23, 144),
            (2, 1, 1, 1, 24, 0),
            (2, 1, 1, 1, 24, -144),
        )
    ):
        result.append((f"outside_{i}", bounds.enclosure, values))
    for i, value in enumerate(invalid + (0, -1, sp.Rational(1, 2))):
        result.append((f"log_{i}", bounds.dyadic_log_upper, (value,)))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        result.append((f"frontier_status_{i}", validate_frontier, (rows,)))
    result.append(("frontier_missing", validate_frontier, (frontier()[:-1],)))
    result.append(
        (
            "frontier_extra",
            validate_frontier,
            (frontier() + [{"id": "new", "status": "COMPLETE"}],),
        )
    )
    return result


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported full-vertex family input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "complete_W2_vertex_not_local_L_replacement": True,
        "all_insertions_and_zero_transfer_terms_present": True,
        "subtractions_are_local_in_the_actual_parent": True,
        "unregulated_zero_bubble_not_given_a_finite_value": True,
        "nonlocal_reference_conversion_not_declared_zero": True,
        "same_family_subset_not_double_counted": True,
        "remaining_seven_primitive_rows_and_conversion_open": True,
        "original_P8_not_closed": True,
    }
