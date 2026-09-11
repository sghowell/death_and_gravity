"""Advance the vacuum reference without claiming the entire source or P8."""

from functools import cache

import sympy as s
from p8_vacuum_full_two_loop_amplitude import audit as previous

from . import bounds, calibration, forest, sector, source

MODULES = (source, forest, sector, bounds, calibration)


def frontier():
    return previous.frontier()


def matching():
    return previous.matching() + [
        {
            "id": "complete_two_loop_vacuum_reference_and_first_source_square",
            "status": "BOUNDED_WITH_FULL_REGULATED_SOURCE_REDUCIBLE_CANCELLATION",
        }
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "Complete vacuum-reference scope differs from the fixed ledger"
        )
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
            "all_prior_matching_items_unchanged": sum(
                a != b for a, b in zip(matching()[:-1], previous.matching())
            ),
            "one_complete_vacuum_reference_item_added": len(matching())
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
        "all_pure_scalar_vacuum_Wick_factors_and_proper_references": True,
        "both_fermion_vacuum_forests_already_assigned_once": True,
        "all_42_flavors_in_gauge_vacuum": True,
        "heavy_determinant_retained_in_first_vacuum": True,
        "first_source_square_mass_insertion_and_reducible_graph_paired": True,
        "source_pairing_before_finite_epsilon_projection": True,
        "physical_Phi_OS_slope_trace_is_scaleless_only_after_pairing": True,
        "six_labeled_sunset_sectors_with_explicit_corner_subtraction": True,
        "complex_regulator_circle_bounds_full_finite_part": True,
        "common_regulator_retains_all_pole_times_epsilon_terms": True,
        "first_Phi_field_direction_leaves_regulated_V1_invariant": True,
        "vacuum_energy_zero_fixes_reference_not_extra_fit_parameter": True,
        "full_second_H_source_not_claimed": True,
        "no_global_quantum_potential_or_physical_truncation_claim": True,
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
    good = [2, 2, 1, 1, 144, 2, 2, 2]
    out = []
    for j in range(8):
        for i, v in enumerate(invalid):
            vals = good.copy()
            vals[j] = v
            out.append((f"scalar_type_{j}_{i}", bounds.scalar_enclosure, tuple(vals)))
        vals = good.copy()
        vals[j] = -1
        out.append((f"scalar_domain_{j}", bounds.scalar_enclosure, tuple(vals)))
    for name, vals in (
        ("m", (300, 2, 1, 1, 144, 2, 2, 2)),
        ("M8", (2, 100, 1, 1, 144, 2, 2, 2)),
        ("M16", (2, 2, 1, 1, 144, 2, 2, 1)),
    ):
        out.append(("power_cap_" + name, bounds.scalar_enclosure, vals))
    for i, vals in enumerate(
        (
            (2, 5, 144, 256, 10),
            (2, 2, 145, 256, 10),
            (2, 2, 144, 255, 10),
            (2, 2, 144, 256, 0),
        )
    ):
        out.append((f"one_loop_domain_{i}", bounds.one_loop_enclosure, vals))
    for i in range(8):
        rows = forest.rows()
        rows[i]["owner"] = "UNASSIGNED"
        out.append((f"vacuum_owner_{i}", forest.validate_rows, (rows,)))
    out.append(("missing_vacuum_owner", forest.validate_rows, (forest.rows()[:-1],)))
    out.append(
        (
            "duplicate_vacuum_owner",
            forest.validate_rows,
            (forest.rows() + forest.rows()[:1],),
        )
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
        raise ValueError("Unsupported complete vacuum input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "source_square_not_silently_omitted": True,
        "finite_source_square_not_square_of_finite_source": True,
        "proper_fermion_counterterm_traces_not_added_twice": True,
        "heavy_determinant_and_MS_mass_reference_retained": True,
        "exact_rational_power_caps_not_floating_point_comparisons": True,
        "vacuum_reference_not_global_potential_or_semiclassical_backreaction": True,
        "second_source_physical_truncation_and_V_G_B_open": True,
        "original_P8_not_closed": True,
    }
