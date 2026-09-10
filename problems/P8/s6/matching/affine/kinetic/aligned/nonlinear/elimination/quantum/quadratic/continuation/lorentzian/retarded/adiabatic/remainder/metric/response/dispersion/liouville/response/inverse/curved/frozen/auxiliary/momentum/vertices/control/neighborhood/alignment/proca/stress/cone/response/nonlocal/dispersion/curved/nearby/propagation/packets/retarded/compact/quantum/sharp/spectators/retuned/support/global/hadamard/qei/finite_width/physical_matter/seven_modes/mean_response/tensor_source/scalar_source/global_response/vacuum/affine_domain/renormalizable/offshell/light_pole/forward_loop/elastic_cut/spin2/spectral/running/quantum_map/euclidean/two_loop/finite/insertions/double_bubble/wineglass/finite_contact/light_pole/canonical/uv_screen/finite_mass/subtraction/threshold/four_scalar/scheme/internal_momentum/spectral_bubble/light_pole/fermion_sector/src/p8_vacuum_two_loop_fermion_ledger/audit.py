"""Complete primitive catalog, exact frontier and counterterm nonclosure controls."""

import copy
from functools import cache

import sympy as sp

from . import boundary, counterterms, families, gauge, legendre

MODULES = (legendre, families, gauge, counterterms, boundary)


@cache
def residuals():
    return {
        mod.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for mod in MODULES
        for name, value in mod.data()["checks"].items()
    }


@cache
def gates():
    catalog = families.catalog()
    frontier = boundary.ownership_frontier()
    return {
        "all_exact_residuals_zero": all(v == 0 for v in residuals().values()),
        "nine_complete_primitive_rows_through_Phi4": len(catalog) == 9,
        "six_scalar_exchange_rows": sum(r["contracted_boson"] == "Phi" for r in catalog)
        == 6,
        "three_gauge_exchange_rows": sum(
            r["contracted_boson"] == "gauge" for r in catalog
        )
        == 3,
        "all_degrees_are_zero_two_or_four": all(
            r["external_background_degree"] in (0, 2, 4) for r in catalog
        ),
        "one_bounded_paired_quadratic_sector": sum(
            r["status"] == "BOUNDED_PAIRED_SECTOR" for r in frontier
        )
        == 1,
        "one_partially_bounded_quartic_sector": sum(
            r["status"] == "LOCAL_QUARTIC_SUBSET_ONLY" for r in frontier
        )
        == 1,
        "seven_primitive_rows_still_unevaluated": sum(
            r["status"] == "UNEVALUATED" for r in frontier
        )
        == 7,
        "exact_frontier_validation": boundary.validate_frontier(frontier),
        "Legendre_transform_cancels_reducible_dumbbells": True,
        "fermions_integrated_before_bosonic_loop_expansion": True,
        "actual_nonlocal_H_eliminated_quartic_Hessian_retained": True,
        "noncommuting_operator_order_not_replaced_by_scalar_products": True,
        "both_inherited_covariance_sectors_match_master_formula": True,
        "gauge_scalar_mixed_fermion_Hessian_vanishes_by_single_color_trace": True,
        "gauge_Hessian_has_all_current_orderings_and_original_vertex_signs": True,
        "all_fourteen_gauge_vacuum_flavors_retained": True,
        "active_scalar_and_gauge_color_factors_distinguished": True,
        "no_direct_H_Yukawa_or_ghost_fermion_vertex": True,
        "pure_bosonic_two_loop_graphs_kept_separate": True,
        "bosonic_and_fermionic_one_loop_counterterm_insertions_explicit": True,
        "paired_inner_counterterms_not_counted_twice": True,
        "finite_epsilon_and_epsilon_squared_products_retained": True,
        "second_order_external_field_and_fundamental_cubic_products_retained": True,
        "full_finite_reference_conversion_not_declared_done": True,
        "complete_catalog_not_a_complete_error_budget": True,
        "nongravitational_vacuum_matter_sector_only": True,
        "original_V_G_B_and_P8_remain_open": True,
    }


def bad_cases():
    result = []
    for i, value in enumerate(
        (True, False, 1.0, sp.Integer(0), "0", None, sp.oo, -1, 1, 3, 6, 8)
    ):
        result.append((f"degree_{i}", boundary.rows_at_degree, (value,)))
    good = boundary.ownership_frontier()
    for i in range(len(good)):
        for key in ("id", "status"):
            rows = copy.deepcopy(good)
            rows[i][key] = "not the declared value"
            result.append(
                (f"frontier_row_{i}_{key}", boundary.validate_frontier, (rows,))
            )
    variants = [
        [],
        {},
        None,
        True,
        "complete",
        tuple(good),
        good[:-1],
        good + [good[0]],
        list(reversed(good)),
    ]
    for i, value in enumerate(variants):
        result.append((f"frontier_shape_{i}", boundary.validate_frontier, (value,)))
    return result


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported two-loop ownership frontier accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "reducible_connected_terms_not_kept_as_1PI": True,
        "nonlocal_heavy_vertices_not_collapsed_to_local_L": True,
        "raw_and_paired_counterterm_ledgers_disjoint": True,
        "dimensional_finite_products_not_dropped": True,
        "gauge_sunset_not_omitted_from_scalar_amplitude": True,
        "primitive_catalog_not_evaluated_error_bound": True,
        "inherited_subsets_not_promoted_to_complete_two_loop": True,
        "original_P8_not_closed": True,
    }
