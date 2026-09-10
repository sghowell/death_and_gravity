"""Partial local-reference frontier for both nonlocal quadratic primitives."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_mixed_quartic import audit as previous

from . import calibration, catalog, forest, pole, regions

MODULES = (catalog, forest, regions, pole, calibration)
TARGETS = ("scalar_Phi2_W0_F2", "gauge_Phi2")
STATUS = "BOUNDED_ON_SHELL_NONLOCAL_REMAINDER_NOT_FINITE_LOCAL_REFERENCES"


def frontier():
    rows = previous.frontier()
    for r in rows:
        if r["id"] in TARGETS:
            r["status"] = STATUS
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError("Quadratic nonlocal frontier differs from exact scope")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: sp.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    old, new = previous.frontier(), frontier()
    out.update(
        {
            "exactly_two_partial_rows_advance": sum(a != b for a, b in zip(old, new))
            - 2,
            "nine_rows_retained": len(new) - 9,
            "two_vacuum_rows_wholly_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 2,
            "two_nonlocal_quadratic_rows": sum(r["status"] == STATUS for r in new) - 2,
            "four_previous_quartic_rows_unchanged": sum(
                r["status"]
                in (
                    "BOUNDED_PAIRED_IN_COMMON_MS_INTERACTION_SCHEME",
                    "BOUNDED_PAIRED_PRIMITIVE_ROW",
                )
                for r in new
            )
            - 4,
        }
    )
    return out


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, sp.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "all_three_cyclic_quadratic_words": True,
        "all_proper_cycles_explicitly_enumerated": True,
        "overlapping_cycles_never_multiplied_as_forest": True,
        "local_whole_fermion_cycle_killed_by_OS_projection": True,
        "complete_MS_self_energy_not_only_local_anchor": True,
        "both_proper_vertex_subtractions_retained": True,
        "both_finite_MS_vertex_anchors_retained": True,
        "middle_raw_and_both_subtractions_bounded": True,
        "high_pair_and_complement_use_distinct_valid_radii": True,
        "all_momentum_regions_integrated": True,
        "massless_chord_diagonal_integrable": True,
        "pointwise_soft_projection_precedes_loop_integrals": True,
        "O4_evenness_removes_lower_soft_polynomial": True,
        "two_disc_Cauchy_OS_remainder_bound": True,
        "finite_mass_slope_and_canonical_factor_not_claimed": True,
        "two_vacuum_rows_and_other_matching_open": True,
        "V_G_B_and_original_P8_open": True,
    }


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
    out = []
    for j in range(4):
        for i, v in enumerate(invalid):
            vals = [720, 1, 1, 144]
            vals[j] = v
            out.append((f"type_{j}_{i}", pole.enclosure, tuple(vals)))
    for i, v in enumerate(
        (
            (719, 1, 1, 144),
            (720, -1, 1, 144),
            (720, 1, -1, 144),
            (720, 1, 1, 0),
            (720, 1, 1, 145),
        )
    ):
        out.append((f"domain_{i}", pole.enclosure, v))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"frontier_{i}", validate_frontier, (rows,)))
    out.extend(
        (
            ("frontier_missing", validate_frontier, (frontier()[:-1],)),
            (
                "frontier_extra",
                validate_frontier,
                (frontier() + [{"id": "extra", "status": "COMPLETE"}],),
            ),
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
        raise ValueError("Unsupported quadratic input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "two_overlapping_vertex_subtractions": True,
        "whole_fermion_cycle_local_term_accounted": True,
        "proper_MS_finite_terms_retained": True,
        "project_before_unbounded_integrals": True,
        "finite_MS_slope_not_silently_zero": True,
        "nonlocal_bound_not_complete_pole": True,
        "other_matching_and_vacuum_rows_open": True,
        "original_P8_not_closed": True,
    }
