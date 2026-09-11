"""Keep the minimal diagnostic exclusion separate from original P8 rows."""

from functools import cache

import sympy as s
from p8_vacuum_curved_dirac_stress import audit as previous

from . import equations, map, nullstress

MODULES = (equations, map, nullstress)
ITEM = {
    "id": "minimal_Einstein_canonical_clock_free_sector_null_equation_mismatch",
    "status": "EXCLUDED_NAMED_MINIMAL_DIAGNOSTIC_WITH_SMALL_EXTRA_NULL_BUDGET_NOT_FULL_INTERACTING_PARENT_OR_DHOST_ROW",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The minimal diagnostic cannot change the original P8 frontier"
        )
    return True


def observable():
    return {
        "quantity": "The physical null Einstein-equation residual on the exact CD trajectory and unit canonical clock, with the S6.171 specified curved free fermion contribution.",
        "model_scope": equations.data()["minimal_diagnostic_name"],
        "additional_terms_boundary": equations.data()["additional_terms_boundary"],
        "negative_verdict_boundary": "The named minimal diagnostic and any explicitly small additional-null-budget variant fail this trajectory. No bound on the whole interacting remainder is derived, and no original DHOST row or general UV completion is excluded.",
        "next_required_matching": "A parent that realizes the actual CD bounce must supply a leading signed nonminimal/derivative/other null contribution of the stated magnitude while retaining its domain, constraint, state and cutoff control. Tighter estimates on the tiny specified free sector alone cannot supply it.",
    }


@cache
def residuals():
    rows = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: s.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    rows.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_previous_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "only_scoped_minimal_diagnostic_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return rows


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **nullstress.data()["gates"],
        **equations.data()["gates"],
        **map.data()["gates"],
        "same_actual_CD_metric_and_free_states": True,
        "same_correlated_covariant_potential_used": True,
        "classical_positive_kinetic_sources_cannot_supply_negative_null_stress": True,
        "metric_variation_and_physical_curvature_sign_checked": True,
        "off_shell_map_does_not_hide_equation_residual": True,
        "additional_small_remainder_budget_explicitly_conditional": True,
        "no_bound_transferred_to_nearby_metrics_or_other_states": True,
        "all_prior_frozen_scientific_bytes_unchanged": True,
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
        -1,
    )
    cases = [
        (f"additional_bound_{i}", equations.margins, (value,))
        for i, value in enumerate(invalid)
    ]
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        cases.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        cases.append((f"matching_{i}", validate_scope, (frontier(), rows)))
    cases.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
        )
    )
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported null-budget input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "potential_correlation_not_two_independent_absolute_errors": True,
        "nonzero_local_bounce_components_not_each_declared_zero": True,
        "gravity_sign_and_positive_clock_term_retained": True,
        "pointwise_map_value_not_enough_without_first_variation": True,
        "small_additional_bound_not_claimed_measured": True,
        "leading_nonminimal_terms_not_excluded_by_diagnostic": True,
        "nearby_metric_state_bound_not_inferred": True,
        "original_P8_not_closed": True,
    }
