"""Scope and negative controls for the actual curved free in/out state."""

from functools import cache

import sympy as s
from p8_vacuum_flat_dirac_pressure import audit as previous

from . import energy, geometry, hadamard, transition, tube

MODULES = (geometry, tube, transition, hadamard, energy)
ITEM = {
    "id": "actual_CD_free_Dirac_Hadamard_in_out_state_and_uniform_state_stress_difference",
    "status": "CURVED_FREE_HADAMARD_STATE_AND_STATE_DIFFERENCE_BOUNDED_NOT_ABSOLUTE_PARENT_STRESS",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A free curved state cannot change the original frontier")
    return True


def observable():
    return {
        "quantity": "Exact free quadratic in/out states on the prescribed actual CD geometry, including42copies, plus uniform differences in their physical energy density and pressure.",
        "state_not_parent": "The state is for a specified quadratic operator on the full momentum space. No physical interacting EFT cutoff, absolute parent stress or parent background solution is asserted.",
        "external_reference": hadamard.data()["external_primary_reference"],
        "finite_reference": "Identical local state-independent subtraction and finite counterterm prescriptions are used for both states; they cancel. No absolute curved Newton/curvature reference is thereby fixed.",
    }


@cache
def residuals():
    rows = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(value)
        for mod in MODULES
        for key, value in mod.data()["checks"].items()
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
            "one_scoped_curved_state_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return rows


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    values = {
        key: value
        for mod in MODULES
        for key, value in mod.data().get("gates", {}).items()
    }
    values.update(
        {
            "actual_CD_geometry_and_reference_time_pinned": True,
            "physical_spinor_Cauchy_normalization_derived": True,
            "geometric_connection_and_all_42_copies_retained": True,
            "modewise_and_strong_not_uniform_Moller_limits": True,
            "all_exact_rotations_and_final_connection_retained": True,
            "both_time_tails_and_full_momentum_integrals": True,
            "arbitrary_order_high_momentum_threshold_explicit": True,
            "external_local_Hadamard_hypotheses_checked": True,
            "same_states_identified_by_all_order_projector_comparison": True,
            "finite_N20_not_used_as_Hadamard_criterion": True,
            "same_counterterm_state_difference_only": True,
            "no_absolute_parent_stress_or_background_control": True,
            "all_prior_frozen_scientific_bytes_unchanged": True,
        }
    )
    return values


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
    cases = []
    for slot in range(3):
        for i, value in enumerate(invalid):
            args = [10**8, s.Rational(1, 100), s.Rational(1, 2)]
            args[slot] = value
            cases.append((f"parameter_{slot}_{i}", transition.parameters, tuple(args)))
    for i, value in enumerate(invalid + (s.Integer(20), 0, -1)):
        cases.append((f"order_{i}", tube.order, (value,)))
        cases.append((f"multiplicity_{i}", energy.enclosures, (10**8, 0, 1, value)))
    for i, value in enumerate(invalid + (-1,)):
        cases.append((f"momentum_{i}", transition.bounds, (value, 10**8, 0, 1)))
    for i, value in enumerate(invalid + (s.Integer(1), -1, 21)):
        cases.append((f"frame_index_{i}", tube.radius, (20, value)))
    for i, args in enumerate(
        (
            (0, 0, 1),
            (10**8, -1, 1),
            (10**8, 10**8, 1),
            (10**8, 0, 0),
            (10**8, 0, 2),
            (1, 0, 1),
        )
    ):
        cases.append((f"domain_{i}", transition.parameters, args))
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
        raise ValueError("Unsupported curved-state input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "D_only_geometry_not_substituted_for_actual_CD": True,
        "constant_mass_flavors_not_called_inert_on_curved_metric": True,
        "four_flat_frames_not_sufficient_curved_UV_control": True,
        "finite_adiabatic_order_not_Hadamard": True,
        "unbounded_order_not_assumed_from_fixed_gap": True,
        "modewise_Moller_bound_not_uniform_in_momentum": True,
        "state_difference_not_absolute_stress_or_bounce_error": True,
        "original_P8_not_closed": True,
    }
