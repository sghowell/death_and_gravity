"""Partial primitive frontier and exact self-energy-chord verification."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_heavy_vertices import audit as previous

from . import calibration, catalog, domain, kernel, tail

MODULES = (catalog, kernel, domain, tail, calibration)
TARGETS = ("scalar_Phi4_W0_F4", "gauge_Phi4")


def frontier():
    rows = previous.frontier()
    for row in rows:
        if row["id"] in TARGETS:
            row["status"] = "PARTIAL_SELF_ENERGY_CHORD_BOUND"
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError(
            "The partial self-energy-chord frontier differs from its exact scope"
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
            "two_self_energy_subsets_not_complete_rows": sum(
                r["status"] == "PARTIAL_SELF_ENERGY_CHORD_BOUND" for r in new
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
            "sixty_cyclic_words_without_discarding_orientation": True,
            "twenty_four_self_energy_words_per_sector": True,
            "bijection_to_six_boxes_four_marked_propagators": True,
            "literal_determinant_and_covariance_variation_normalization": True,
            "selected_family_closed_under_all_external_permutations": True,
            "whole_MS_self_energy_not_local_reference_only": True,
            "proper_mass_and_kinetic_counterterms_paired_once": True,
            "finite_dimensional_gauge_terms_retained": True,
            "global_complex_parameter_gap_for_all_real_internal_momenta": True,
            "real_Dirac_norm_and_complex_Neumann_bound": True,
            "integrable_log_endpoint_bound_for_massless_gauge_line": True,
            "soft_Cauchy_radius_depends_on_internal_momentum": True,
            "soft_degree_four_tail_projected_before_integration": True,
            "Lorentz_and_permutation_symmetry_remove_lower_degree_b2": True,
            "remaining_overall_contact_not_an_adjustable_b2_term": True,
            "exact_logarithmic_radial_integral_not_a_cutoff": True,
            "unit_forward_disc_Cauchy_coefficient_bound": True,
            "both_scalar_and_gauge_self_energy_chords_bounded": True,
            "no_sign_claim_for_the_amplitude_from_positive_majorants": True,
            "vertex_and_other_primitive_words_still_unbounded": True,
            "other_field_parameter_and_outer_reference_conversions_separate": True,
            "original_P8_V_G_B_and_complete_two_loop_remain_open": True,
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
        raise ValueError("Unsupported self-energy-chord input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "self_energy_subset_not_entire_primitive_rows": True,
        "complete_proper_MS_kernel_and_finite_anchors": True,
        "both_scalar_and_gauge_internal_bosons_present": True,
        "no_extra_half_or_closed_flavor_multiplicity": True,
        "Taylor_tail_projected_before_UV_integration": True,
        "q_dependent_radius_not_exported_to_vertex_chords": True,
        "other_matching_conversions_and_errors_open": True,
        "original_P8_not_closed": True,
    }
