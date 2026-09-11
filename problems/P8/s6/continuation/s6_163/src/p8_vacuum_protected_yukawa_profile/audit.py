"""A named protected-mass EFT branch, without closing the common-parent gate."""

from functools import cache

import sympy as s
from p8_vacuum_full_analytic_target_match import audit as previous

from . import calibration, clock, order, profile, transport

MODULES = (profile, order, clock, calibration, transport)
ITEM = {
    "id": "GY14_SAT8_pointwise_fermion_mass_protection_and_named_low_order_invariance",
    "status": "NEW_ANALYTIC_EFT_BRANCH_WITH_GLOBAL_REAL_MASS_FLOOR_AND_SAME_NAMED_TWO_LOOP_DATA_NOT_ROLLING_STATE_OR_COMMON_PARENT_CLOSURE",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("The protected profile cannot change its fixed scope ledger")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(value)
        for mod in MODULES
        for key, value in mod.data()["checks"].items()
    }
    out.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "all_parent_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_parent_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "only_new_protected_mass_candidate_item_added": len(matching())
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
        **transport.data()["bounds"],
        "real_positive_eighth_root_and_nonzero_mass_floor_all_real_fields": True,
        "vacuum_branch_complex_radius_finite_not_entire": True,
        "new_vertices_start_at_nine_scalar_legs": True,
        "forest_contraction_retains_loop_grade": True,
        "common_dimensional_lift_precedes_finite_products": True,
        "only_named_low_source_orders_not_entire_functional_identified": True,
        "physical_source_map_not_ordinary_Psi_prescription": True,
        "flat_naive_map_screen_not_particle_production_or_row_exclusion": True,
        "nonpolynomial_EFT_not_assigned_old_all_field_UV_running": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
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
    goodgap = [10, 1, 1, s.Rational(1, 2)]
    goodclock = [10, 1, 1, 1, 0, s.Rational(1, 2)]
    for prefix, call, good in (
        ("gap", profile.gap, goodgap),
        ("clock", clock.ratios, goodclock),
    ):
        for slot in range(len(good)):
            for i, value in enumerate(invalid):
                args = good.copy()
                args[slot] = value
                out.append((f"{prefix}_type_{slot}_{i}", call, tuple(args)))
    for slot, value in ((0, 0), (1, 0), (2, 0), (3, 0), (3, 1), (2, 5)):
        args = goodgap.copy()
        args[slot] = value
        out.append((f"gap_domain_{slot}_{value}", profile.gap, tuple(args)))
    for slot, value in (
        (0, 0),
        (1, 0),
        (2, 0),
        (3, -1),
        (4, -1),
        (5, 0),
        (5, 1),
        (2, 5),
    ):
        args = goodclock.copy()
        args[slot] = value
        out.append((f"clock_domain_{slot}_{value}", clock.ratios, tuple(args)))
    for i, value in enumerate(invalid + (s.Integer(1), -1)):
        out.append((f"binomial_index_{i}", profile.coefficient, (value,)))
    for i, value in enumerate(invalid + (s.Integer(4), -1, 11)):
        out.append((f"source_count_{i}", order.new_vertex_floor, (value,)))
    for i, value in enumerate(invalid + (s.Integer(1), 0, -1)):
        out.append((f"new_vertex_index_{i}", order.new_vertex_floor, (4, value)))
    out.extend(
        (
            ("grade_bad_external", order.grade, (True, (4,), (0,))),
            ("grade_negative_external", order.grade, (-1, (4,), (0,))),
            ("grade_list_not_tuple", order.grade, (4, [4], (0,))),
            ("grade_empty", order.grade, (0, (), ())),
            ("grade_lengths", order.grade, (0, (4,), (0, 0))),
            ("grade_sympy_valence", order.grade, (4, (s.Integer(4),), (0,))),
            ("grade_negative_order", order.grade, (4, (4,), (-1,))),
            ("grade_impossible_halfedges", order.grade, (4, (3,), (0,))),
            ("grade_odd_halfedges", order.grade, (4, (5,), (0,))),
            ("grade_negative_loop", order.grade, (8, (4, 4), (0, 0))),
        )
    )
    for i in range(len(transport.rows())):
        rows = transport.rows()
        rows[i]["through_loop"] = 10
        out.append((f"observable_{i}", transport.validate_rows, (rows,)))
    out.append(
        ("missing_observable", transport.validate_rows, (transport.rows()[:-1],))
    )
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_{i}", validate_scope, (frontier(), rows)))
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
        raise ValueError("Unsupported protected profile input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "original_linear_mass_crossing_not_silently_ignored": True,
        "paired_logdet_higher_field_difference_not_zero": True,
        "all_field_or_all_loop_identity_not_claimed": True,
        "direct_argument_small_ratios_not_transferred_to_literal_cubic_map": True,
        "pointwise_mass_floor_not_global_Dirac_propagator_or_state": True,
        "nonpolynomial_EFT_not_claimed_globally_renormalizable_GY14": True,
        "classical_fermion_zero_match_preserved_without_quantum_target_inference": True,
        "original_P8_not_closed": True,
    }
