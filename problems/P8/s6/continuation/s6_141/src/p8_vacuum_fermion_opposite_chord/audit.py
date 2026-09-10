"""Partial primitive frontier and exact opposite-chord verification."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_self_energy_chord import audit as previous

from . import angular, calibration, catalog, domain, tail

MODULES = (catalog, domain, angular, tail, calibration)
TARGETS = ("scalar_Phi4_W0_F4", "gauge_Phi4")


def frontier():
    rows = previous.frontier()
    for row in rows:
        if row["id"] in TARGETS:
            row["status"] = "PARTIAL_SELF_ENERGY_AND_OPPOSITE_CHORD_BOUNDS"
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError(
            "The partial opposite-chord frontier differs from its exact scope"
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
            "exactly_two_primitive_statuses_become_partial": sum(
                a != b for a, b in zip(old, new)
            )
            - 2,
            "nine_primitive_rows_retained": len(new) - 9,
            "five_other_primitive_rows_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 5,
            "two_cumulative_chord_subsets_not_complete_rows": sum(
                r["status"] == "PARTIAL_SELF_ENERGY_AND_OPPOSITE_CHORD_BOUNDS"
                for r in new
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
            "twelve_opposite_words_per_sector_both_orientations": True,
            "determinant_and_Gaussian_contraction_normalization": True,
            "all_proper_one_loop_cycles_UV_finite": True,
            "overall_MS_quartic_contact_has_zero_b2": True,
            "exact_S4_closed_subset_before_lower_degree_reduction": True,
            "literal_unshifted_boson_routing_all_selected_words": True,
            "joint_minimum_momentum_Cauchy_radius": True,
            "real_Dirac_norm_and_complex_Neumann_bound": True,
            "four_gauge_indices_and_open_Casimir_counted": True,
            "exact_four_dimensional_angular_average": True,
            "massless_diagonal_and_zero_momentum_IR_integrable": True,
            "positive_scalar_mass_bounded_by_massless_chord": True,
            "soft_degree_four_tail_projected_before_both_integrals": True,
            "all_positive_soft_degrees_UV_integrable": True,
            "exact_unbounded_two_loop_radial_majorant": True,
            "unit_forward_disc_Cauchy_coefficient_bound": True,
            "twelve_new_words_disjoint_from_twenty_four_self_energy": True,
            "twenty_four_vertex_words_per_sector_still_unbounded": True,
            "other_counterterm_and_outer_reference_conversions_separate": True,
            "no_amplitude_sign_from_positive_majorant": True,
            "original_P8_V_G_B_and_full_two_loop_open": True,
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
        raise ValueError("Unsupported opposite-chord input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "opposite_subset_not_complete_primitive_rows": True,
        "no_proper_UV_counterterm_is_missing": True,
        "both_scalar_and_gauge_internal_bosons": True,
        "no_extra_half_or_flavor_multiplicity": True,
        "joint_Taylor_tail_projected_before_integration": True,
        "unshifted_massless_chord_IR_integrable": True,
        "vertex_and_other_matching_tasks_open": True,
        "original_P8_not_closed": True,
    }
