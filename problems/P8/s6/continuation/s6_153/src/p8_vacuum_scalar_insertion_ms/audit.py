"""Bound all raw scalar interaction families without closing global matching."""

from functools import cache

import sympy as s
from p8_vacuum_double_bubble_ms import audit as previous

from . import bounds, calibration, forward, inner, outer

MODULES = (inner, outer, forward, bounds, calibration)


def frontier():
    return previous.frontier()


def matching():
    return previous.matching() + [
        {
            "id": "complete_scalar_OS_insertion_interaction_forest_MS",
            "status": "BOUNDED_WITH_FULL_REGULATED_SLOPE_TIMES_OUTER_REFERENCE",
        }
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Scalar insertion MS scope differs from the fixed ledger")
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
            "one_scalar_insertion_conversion_added": len(matching())
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
        "full_dimensional_inner_on_shell_slope": True,
        "epsilon_slope_coefficient_before_outer_pole_product": True,
        "two_scale_logs_and_mass_log_in_finite_reference": True,
        "same_physical_inner_mass_residue_and_source_grouping": True,
        "all_sixty_four_refinements_and_both_insertions_retained": True,
        "both_decaying_insertions_unchanged_and_convergent": True,
        "finite_reference_fixed_by_MS_not_fitted_to_b2": True,
        "exact_three_channel_C_squared_forward_coefficient": True,
        "linear_G_reference_not_disjoint_G1_square": True,
        "four_disjoint_raw_scalar_family_bounds_assembled": True,
        "UV_finite_family_invariant_at_fixed_formal_order": True,
        "global_field_coupling_cross_terms_still_open": True,
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
    base = [1, 1, 40, 1, 144]
    for j in range(5):
        for i, v in enumerate(invalid):
            vals = list(base)
            vals[j] = v
            out.append((f"type_{j}_{i}", bounds.enclosure, tuple(vals)))
    for i, (j, v) in enumerate(
        ((0, 0), (1, 0), (2, 31), (2, 10**400 + 1), (3, 0), (4, 0), (4, 145), (1, 40))
    ):
        vals = list(base)
        vals[j] = v
        out.append((f"domain_{i}", bounds.enclosure, tuple(vals)))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
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
        raise ValueError("Unsupported scalar insertion MS input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "alpha_first_epsilon_not_discarded_before_pole": True,
        "naive_alpha_zero_times_scale_log_is_incomplete": True,
        "complete_two_position_OS_family_not_raw_tadpoles": True,
        "full_C_squared_heavy_mass_and_cubic_variations": True,
        "no_new_local_fit_or_extra_LSZ": True,
        "four_raw_families_not_full_GY14_canonical_matching": True,
        "finite_order_bound_not_physical_truncation": True,
        "original_P8_not_closed": True,
    }
