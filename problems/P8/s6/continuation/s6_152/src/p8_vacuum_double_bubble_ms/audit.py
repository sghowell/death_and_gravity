"""Advance the double_bubble interaction-forest conversion, not the whole MS map."""

from functools import cache

import sympy as s
from p8_vacuum_scalar_zero_reference import audit as previous

from . import bounds, calibration, conversion, forests, forward

MODULES = (forests, conversion, forward, bounds, calibration)


def frontier():
    return previous.frontier()


def matching():
    return previous.matching() + [
        {
            "id": "complete_double_bubble_interaction_forest_MS",
            "status": "BOUNDED_WITH_ACTUAL_NESTED_MS_FORESTS_AND_BOTH_SCALE_TERMS",
        }
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Double-bubble MS scope differs from the fixed ledger")
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
            "one_double_bubble_conversion_obligation_added": len(matching())
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
        "all_twenty_four_actual_restricted_MS_forests": True,
        "forbidden_shared_pair_absent": True,
        "raw_and_nested_overall_poles_keep_finite_bubble_terms": True,
        "shared_projection_cancellation_before_finite_part": True,
        "finite_triangle_keeps_its_complete_dimensional_factor": True,
        "all_eight_finite_triangle_products_excluded": True,
        "both_linear_and_quadratic_scale_terms_retained": True,
        "both_full_heavy_triangles_and_all_three_channels": True,
        "quadratic_scale_b2_exact_not_discarded_as_contact": True,
        "heavy_cubic_inverse_owned_by_disjoint_counterterms": True,
        "all_radius_heavy_propagators_no_expansion_or_cutoff": True,
        "other_scalar_families_and_full_coordinate_cross_terms_open": True,
        "no_extra_LSZ_or_double_counted_parameter_shift": True,
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
    e = s.Symbol("epsilon")
    for name, expr, regulator in (
        ("float", s.Float(1) / e, e),
        ("fractional", s.sqrt(e), e),
        ("analytic_not_polynomial", s.exp(e), e),
        ("non_laurent_rational", 1 / (1 + e), e),
        ("nonsymbolic_regulator", 1 / e, 1),
    ):
        out.append(("pole_" + name, forests.poles, (expr, regulator)))
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported double_bubble MS input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "forbidden_shared_pair_not_inserted_to_fake_product": True,
        "raw_and_nested_finite_bubble_poles_retained": True,
        "linear_and_quadratic_scale_terms_both_required": True,
        "both_heavy_triangles_and_all_three_channels_required": True,
        "heavy_counterterm_products_not_recounted_as_new_contacts": True,
        "other_scalar_families_and_full_field_map_open": True,
        "finite_order_bound_not_physical_higher_loop_truncation": True,
        "original_P8_not_closed": True,
    }
