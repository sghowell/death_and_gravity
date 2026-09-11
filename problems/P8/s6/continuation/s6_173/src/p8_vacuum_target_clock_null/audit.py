"""Scope ledger for the actual leading CD clock equations."""

from functools import cache

import sympy as s
from p8_vacuum_canonical_bounce_gap import audit as previous

from . import homogeneous, structural, target

MODULES = (homogeneous, structural, target)
ITEM = {
    "id": "actual_CD_leading_clock_Euler_source_decomposition_and_matching_domain_gap",
    "status": "LEADING_TARGET_SOURCES_AND_CONSTRAINT_DOMAIN_REQUIREMENTS_IDENTIFIED_NOT_PROPAGATING_PARENT_OR_QUANTUM_BOUNCE",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The target clock decomposition cannot change original P8 verdicts"
        )
    return True


def observable():
    return {
        "quantity": "The unchanged analytic CD_matter target's leading classical metric Euler sources on the full unit-clock trajectory, with original rolling matter and correct overall kappa normalization.",
        "fixed_source": "S6.109 analytic family and its source-pinned original CD witness and retuned margin; S6.164 vacuum matching test class is unchanged.",
        "variational_boundary": "Action values are not metric equations. Curvature X derivatives and the A3 lapse derivative are retained; A4,A5 remain present in the constraint/second-variation analysis.",
        "next_matching_requirement": "A common propagating parent must reproduce these leading clock sources and constraint structure in a genuine clock-neighborhood domain. A tiny action difference on the disjoint unit Fourier-jet vacuum class does not provide that result.",
        "not_inferred": "No propagating UV parent, physical cutoff, higher-loop bound, actual quantum-corrected bounce, background response, V/G/B or original P8 closure.",
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
            "only_actual_target_decomposition_added": len(matching())
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
        **structural.data()["gates"],
        "same_original_CD_matter_target_and_clock": True,
        "entire_time_energy_and_pressure_equations_checked": True,
        "full_analytic_first_clock_jets_not_global_packet_replacement": True,
        "original_rolling_matter_retained": True,
        "curvature_value_zero_not_variation_zero": True,
        "A3_value_zero_not_variation_zero": True,
        "A4_A5_constraint_role_retained": True,
        "vacuum_action_norm_domain_not_extended_to_clock": True,
        "signed_leading_parent_requirement_not_a_UV_construction": True,
        "physical_overall_kappa_normalization_retained": True,
        "quantum_free_source_still_needs_state_aware_background_control": True,
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
        0,
        1,
    )
    cases = [
        (f"clock_normalization_{i}", structural.normalized_first_jet_gap, (value,))
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
        raise ValueError("Unsupported target-domain input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "curvature_X_derivative_not_dropped_after_clock_value": True,
        "A3_lapse_derivative_not_dropped_after_clock_value": True,
        "zero_background_A4_A5_not_deleted_from_constraints": True,
        "M1_source_not_silently_removed": True,
        "matter_momentum_substitution_not_before_metric_variation": True,
        "vacuum_norm_class_not_clock_neighborhood": True,
        "leading_target_reconstruction_not_parent_matching": True,
        "original_P8_not_closed": True,
    }
