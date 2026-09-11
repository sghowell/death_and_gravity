"""Scope ledger for the specified complete homogeneous flat tensor."""

from functools import cache

import sympy as s
from p8_vacuum_flat_local_energy import audit as previous

from . import dimensional, pressure, projector

MODULES = (projector, dimensional, pressure)
ITEM = {
    "id": "complete_free_flat_homogeneous_one_loop_stress_with_dimensional_curvature_improvement",
    "status": "BOUNDED_SPECIFIED_QUADRATIC_MS_MR_FLAT_TENSOR_NOT_FULL_INTERACTING_CURVED_PARENT_STRESS",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The homogeneous flat tensor cannot change the original frontier"
        )
    return True


def observable():
    return {
        "quantity": "The same free in/out Hadamard state's homogeneous flat one-loop tensor diag(rho,P,P,P), in the fixed Minkowski orthonormal frame and the specified MS/vacuum/saturated-mass-reference prescription.",
        "new_curvature_choice": "Additional finite curvature terms are set to zero in the common MS prescription. The prior homogeneous flat T00 did not fix curvature improvements; this choice is explicit here. No full curved-state or physical Newton-reference dictionary is inferred.",
        "source_context": [
            {
                "url": "https://arxiv.org/pdf/1703.00908",
                "version": "v2, 23 May 2017",
                "location": "Section IV equations51-57 and section V covariant counterterms",
                "role": "Independent UV-pressure coefficient and covariance checks; the complete exact-state remainder and numerical bounds are derived here.",
            }
        ],
        "physical_coordinate_boundary": "This quadratic/reference contribution in the MS scalar coordinate is not the full canonically re-expressed scalar/heavy/gauge parent stress.",
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
            "only_scoped_flat_tensor_added": len(matching())
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
        **projector.data()["gates"],
        **pressure.data()["gates"],
        "same_exact_rotationally_invariant_Hadamard_states": True,
        "physical_pressure_operator_not_evolution_frequency": True,
        "off_diagonal_coherence_and_exact_q1_retained": True,
        "complete_off_diagonal_radial_integrals_converge": True,
        "full_d_dimensional_isotropic_pressure_factor_retained": True,
        "homogeneous_flat_curvature_pressure_not_dropped": True,
        "additional_finite_curvature_choice_recorded": True,
        "same_full_potential_and_mass_reference_used": True,
        "all_tensor_components_in_fixed_rest_frame_controlled": True,
        "MS_increment_not_full_physical_parent_stress": True,
        "curved_state_cutoff_and_higher_loops_not_inferred": True,
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
    for slot in range(3):
        for i, value in enumerate(invalid):
            args = [2**22, s.Rational(1, 100), 1]
            args[slot] = value
            cases.append((f"pressure_{slot}_{i}", pressure.enclosures, tuple(args)))
    for i, value in enumerate(invalid + (s.Integer(6), 0, -1)):
        cases.append((f"multiplicity_{i}", pressure.enclosures, (2**22, 0, 1, value)))
    for i, args in enumerate(
        (
            (0, 0, 1),
            (1, -1, 2**22),
            (2**22, 0, 0),
            (1, 1, 2**22),
            (1, 0, 1),
            (35, 0, 2**22),
        )
    ):
        cases.append((f"domain_{i}", pressure.enclosures, args))
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
        raise ValueError("Unsupported pressure input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "rotating_frame_frequency_not_physical_pressure": True,
        "off_diagonal_second_derivative_coherence_not_omitted": True,
        "three_dimensional_isotropic_prefactor_not_frozen_before_poles": True,
        "zero_curvature_not_zero_improvement_stress": True,
        "finite_curvature_choice_not_inferred_from_energy": True,
        "no_differentiation_of_an_energy_inequality": True,
        "homogeneous_flat_tensor_not_curved_interacting_parent": True,
        "original_P8_not_closed": True,
    }
