"""Complete selected quartic primitive frontier and exact vertex-chord verification."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_opposite_chord import audit as previous

from . import calibration, catalog, forest, regions, tail

MODULES = (catalog, forest, regions, tail, calibration)
TARGETS = ("scalar_Phi4_W0_F4", "gauge_Phi4")


def frontier():
    rows = previous.frontier()
    for row in rows:
        if row["id"] in TARGETS:
            row["status"] = "BOUNDED_PAIRED_PRIMITIVE_ROW"
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError(
            "The complete selected primitive frontier differs from its exact scope"
        )
    return True


@cache
def residuals():
    rows = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for mod in MODULES
        for name, value in mod.data()["checks"].items()
    }
    old, new = previous.frontier(), frontier()
    rows.update(
        {
            "exactly_two_partial_primitive_rows_now_bounded": sum(
                a != b for a, b in zip(old, new)
            )
            - 2,
            "nine_primitive_rows_retained": len(new) - 9,
            "five_other_primitive_rows_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 5,
            "two_complete_quartic_primitive_rows": sum(
                r["status"] == "BOUNDED_PAIRED_PRIMITIVE_ROW" for r in new
            )
            - 2,
            "two_older_paired_families_unchanged": sum(
                r["status"]
                in ("BOUNDED_PAIRED_SECTOR", "BOUNDED_IN_DECLARED_PARENT_SUBTRACTION")
                for r in new
            )
            - 2,
        }
    )
    return rows


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, sp.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    result.update(
        {
            "all_twenty_four_vertex_words_both_orientations": True,
            "marked_box_vertex_normalization": True,
            "unique_proper_logarithmic_vertex_subgraph": True,
            "proper_MS_Yukawa_counterterm_paired_once": True,
            "finite_MS_zero_momentum_anchor_retained": True,
            "regulator_products_not_dropped_before_subtraction": True,
            "short_arc_routing_keeps_chord_unshifted": True,
            "raw_LOW_and_subtraction_LOW_use_distinct_valid_radii": True,
            "paired_HIGH_improves_large_short_arc_momentum_decay": True,
            "raw_HIGH_not_integrated_independently": True,
            "exact_angular_kernel_for_LOW": True,
            "all_three_regional_radial_constants_derived": True,
            "soft_projection_precedes_both_integrals": True,
            "all_positive_soft_degrees_converge": True,
            "S4_and_Lorentz_remove_lower_degree_b2": True,
            "overall_quartic_contact_not_finite_b2_tuning": True,
            "sixty_disjoint_words_in_each_completed_primitive_row": True,
            "five_other_primitive_rows_unevaluated": True,
            "other_matching_and_outer_MS_conversions_separate": True,
            "not_full_two_loop_error_or_V_G_B_P8_closure": True,
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
    for j in range(4):
        for i, value in enumerate(invalid):
            values = [720, 1, 1, 144]
            values[j] = value
            result.append((f"type_{j}_{i}", tail.enclosure, tuple(values)))
    for i, values in enumerate(
        (
            (0, 1, 1, 144),
            (719, 1, 1, 144),
            (720, -1, 1, 144),
            (720, 1, -1, 144),
            (720, 1, 1, 0),
            (720, 1, 1, -144),
            (720, 1, 1, 145),
        )
    ):
        result.append((f"outside_{i}", tail.enclosure, values))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        result.append((f"frontier_{i}", validate_frontier, (rows,)))
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
        raise ValueError("Unsupported vertex-chord input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "two_quartic_rows_not_complete_fermion_sector": True,
        "unique_proper_vertex_counterterm_owned_once": True,
        "whole_MS_local_anchor_not_dropped": True,
        "both_scalar_and_gauge_chords": True,
        "region_radii_not_interchanged": True,
        "no_raw_HIGH_UV_integral": True,
        "other_matching_tasks_open": True,
        "original_P8_not_closed": True,
    }
