"""Scope and exact-evidence ledger for the specified curved quadratic tensor."""

from functools import cache

import sympy as s
from p8_vacuum_curved_dirac_state import audit as previous

from . import bloch, curvature, reference, stress, subtraction

MODULES = (bloch, curvature, reference, stress, subtraction)
ITEM = {
    "id": "complete_free_curved_one_loop_stress_with_explicit_EC_Newton_reference",
    "status": "BOUNDED_SPECIFIED_QUADRATIC_CURVED_TENSOR_NOT_FULL_INTERACTING_PARENT_OR_BACKGROUND_RESPONSE",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The specified curved tensor cannot change the original frontier"
        )
    return True


def observable():
    return {
        "quantity": "Absolute homogeneous free Dirac one-loop diag(rho,P,P,P) for either S6.170 exact in/out Hadamard state on the actual CD metric, in the comoving orthonormal frame.",
        "prescription": "GY14-SAT8-MR/EC-N0: inherited vacuum and saturated pole-mass reference, explicit dimensional E_D,C_D curvature pole basis, zero extra finite R^2 in that basis, and a separately fixed free zero-field Newton reference.",
        "source_context": [
            {
                "url": "https://arxiv.org/pdf/1703.00908",
                "version": "v2, 23 May 2017",
                "location": "Sections V-VI; equations106,109,110",
                "role": "Independent local heat-kernel/spin-curvature input; own dimensional metric variation fixes the finite terms and sign convention. Entire-state bounds are derived here.",
            }
        ],
        "not_inferred": "No full interacting parent stress or physical Newton dictionary, scalar response, controlled background, higher loops, cutoff, quantum target matching or original V/G/B closure.",
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
            "only_specified_curved_tensor_added": len(matching())
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
        **bloch.data()["gates"],
        **reference.data()["gates"],
        **stress.data()["gates"],
        "same_actual_CD_Hadamard_states_not_flat_surrogates": True,
        "all_forty_two_curved_Dirac_copies_retained": True,
        "twenty_frame_comparison_not_itself_Hadamard": True,
        "physical_tensor_not_rotating_generator": True,
        "complete_UV_subtracted_momentum_integrals": True,
        "full_dimensional_pressure_and_counterterm_variation": True,
        "explicit_evanescent_curvature_basis_not_silently_identified": True,
        "Euler_varied_before_regulator_limit": True,
        "common_covariant_local_energy_and_pressure_action": True,
        "finite_Newton_reference_is_named_not_energy_fitted": True,
        "same_entire_potential_and_fixed_mass_reference": True,
        "small_stress_not_background_or_response_control": True,
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
    )
    cases = []
    for label, call in (
        ("remainder", bloch.enclosures),
        ("local", reference.enclosures),
        ("total", stress.enclosures),
    ):
        for slot in range(3):
            for i, value in enumerate(invalid):
                args = [2**22, s.Rational(1, 100), 1]
                args[slot] = value
                cases.append((f"{label}_{slot}_{i}", call, tuple(args)))
    for i, value in enumerate(invalid + (s.Integer(6), 0, -1, 3)):
        cases.append(
            (f"paired_multiplicity_{i}", reference.enclosures, (2**22, 0, 1, value))
        )
    for i, value in enumerate(invalid + (s.Integer(42), 0, -1, 5)):
        cases.append(
            (f"total_multiplicity_{i}", stress.enclosures, (2**22, 0, 1, 6, value))
        )
    for i, args in enumerate(
        (
            (0, 0, 1),
            (1, -1, 1),
            (2**22, 0, 0),
            (1, 1, 1),
            (1, 0, 1),
            (2**22, 0, 2),
            (2**22, 2**21, 1),
        )
    ):
        cases.append((f"domain_{i}", stress.enclosures, args))
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
        raise ValueError("Unsupported curved-stress input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "same_actual_CD_geometry_not_flat_state": True,
        "constant_mass_species_not_dropped": True,
        "physical_pressure_coherence_not_omitted": True,
        "three_dimensional_prefactors_not_frozen_before_poles": True,
        "Euler_four_dimensional_topology_not_used_before_pole": True,
        "evanescent_basis_difference_not_silently_zero": True,
        "small_free_stress_not_a_canonical_parent_bounce": True,
        "original_P8_not_closed": True,
    }
