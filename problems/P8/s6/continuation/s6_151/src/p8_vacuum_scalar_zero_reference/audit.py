"""Advance the wineglass interaction-forest conversion, not the whole MS map."""

from functools import cache

import sympy as s
from p8_vacuum_finite_contact_conversion import audit as previous

from . import anchors, bounds, calibration, conversion, sector

MODULES = (sector, anchors, conversion, bounds, calibration)


def frontier():
    return previous.frontier()


def matching():
    return previous.matching() + [
        {
            "id": "complete_wineglass_interaction_forest_MS",
            "status": "BOUNDED_WITH_REGULATED_FINITE_LOCAL_ANCHOR_AND_FULL_SCALE_TERM",
        }
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Wineglass MS scope differs from the fixed ledger")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: s.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    out.update(
        {
            "all_nine_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "all_prior_matching_obligations_unchanged": sum(
                a != b for a, b in zip(matching()[:-1], previous.matching())
            ),
            "one_wineglass_conversion_obligation_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return out


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "equal_mass_derivative_in_convergence_domain_then_continuation": True,
        "six_actual_parameter_sectors_and_endpoint_subtraction": True,
        "finite_sector_constant_convergent_with_explicit_bound": True,
        "recursive_and_proper_MS_references_distinguished": True,
        "I0_epsilon_coefficient_retained_in_pole_products": True,
        "whole_sixteen_vertex_numerator_and_both_heavy_triangles": True,
        "finite_inner_scale_difference_multiplies_complete_outer_integral": True,
        "local_reference_b2_uses_exact_heavy_tree_variation": True,
        "all_radius_heavy_propagators_no_expansion_or_cutoff": True,
        "other_scalar_families_and_full_coordinate_cross_terms_open": True,
        "no_extra_LSZ_on_canonical_family": True,
        "finite_sector_bound_not_physical_higher_loop_truncation": True,
        "original_V_G_B_and_P8_remain_open": True,
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    out = []
    base = [s.Rational(1, 100), 1, 1200, 600, 144]
    for j in range(5):
        for i, v in enumerate(invalid):
            vals = list(base)
            vals[j] = v
            out.append((f"type_{j}_{i}", bounds.enclosure, tuple(vals)))
    for i, (j, v) in enumerate(
        ((0, 0), (1, -1), (2, -1), (2, 1201), (3, -1), (3, 601), (4, 0), (4, 145))
    ):
        vals = list(base)
        vals[j] = v
        out.append((f"domain_{i}", bounds.enclosure, tuple(vals)))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching())))
        obligations = matching()
        obligations[i]["status"] = "COMPLETE"
        out.append((f"matching_{i}", validate_scope, (frontier(), obligations)))
    out.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_obligation", validate_scope, (frontier(), matching()[:-1])),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported wineglass MS input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "recursive_J0_not_substituted_for_proper_MS_anchor": True,
        "finite_epsilon_times_pole_products_not_discarded": True,
        "complete_inner_scale_term_not_only_local_anchor": True,
        "both_heavy_triangles_and_all_three_channels_required": True,
        "local_tree_reference_not_counted_twice_as_coordinate_shift": True,
        "other_scalar_families_and_full_field_map_open": True,
        "sector_constant_bound_not_physical_truncation": True,
        "original_P8_not_closed": True,
    }
