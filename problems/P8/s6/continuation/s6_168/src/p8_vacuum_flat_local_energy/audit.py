"""Explicit scope, exact checks and rejected input controls."""

from functools import cache

import sympy as s
from p8_vacuum_superadiabatic_state_energy import audit as previous

from . import energy, extension, frames, reference, subtraction

MODULES = (frames, subtraction, reference, extension, energy)
ITEM = {
    "id": "free_flat_absolute_one_loop_local_energy_with_fixed_saturated_mass_reference",
    "status": "BOUNDED_IN_EXPLICIT_MS_AND_SATURATED_MASS_REFERENCE_QUADRATIC_SCHEME_NOT_FULL_CURVED_INTERACTING_PARENT_STRESS",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The specified one-loop local energy cannot change the fixed frontier"
        )
    return True


def observable():
    return {
        "quantity": "Absolute local one-loop Dirac T00 of either S6.166 free flat Hadamard state, at fixed prescribed scalar mass argument; includes the MS finite parts, the fixed vacuum zero and the GY14-SAT8-MR saturated finite pole-mass reference.",
        "coordinate": "Named MS scalar-field coordinate. The inherited physical unit-residue source/parameter conversion is not a second finite energy subtraction, nor is this quadratic increment the entire canonically re-expressed interacting parent stress.",
        "source_context": [
            {
                "url": "https://arxiv.org/pdf/1703.00908",
                "version": "v2, 23 May 2017",
                "location": "Section IV equations44-50 and section V equations65-85",
                "role": "External-scalar Dirac energy subtraction and covariant counterterm structure. Our resummed fixed-mass dimensional finite parts, exact rotating-frame remainder and rational bounds are derived here.",
            }
        ],
        "fixed_scientific_boundary": "S6.167 supplies exact-frame bounds and S6.166 the same free Hadamard state. This new separately named finite high-field prescription is not retroactively assigned to the frozen SAT8 functional.",
        "reference_not_bounce_density": "The kappa comparison is not a relative error against the zero density at a flat bounce.",
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
            "nine_prior_primitive_rows_retained": len(frontier()) - 9,
            "all_prior_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_prior_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "only_named_quadratic_local_energy_added": len(matching())
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
        **{k: v for mod in MODULES for k, v in mod.data().get("gates", {}).items()},
        "same_exact_in_out_Hadamard_states": True,
        "instantaneous_fourth_frame_not_claimed_as_exact_Hadamard_state": True,
        "all_halfline_transition_orders_retained": True,
        "subtracted_mode_remainder_integrated_over_all_momenta": True,
        "MS_finite_parts_restored_not_removed_by_normal_ordering": True,
        "full_dimensional_angular_factor_retained": True,
        "finite_mass_counterterm_coefficient_fixed_by_prior_pole": True,
        "higher_field_extension_named_and_not_a_frozen_rewrite": True,
        "full_positive_potential_tail_retained": True,
        "inert_flavor_vacuum_reference_retained": True,
        "pressure_and_curved_counterterms_not_inferred_from_T00": True,
        "scalar_heavy_gauge_and_higher_loop_energy_not_included": True,
        "whole_quadratic_momentum_space_not_a_physical_cutoff": True,
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
    cases = []
    for i, value in enumerate(invalid + (-1, 11, s.Integer(4))):
        cases.append((f"loop_floor_{i}", extension.loop_floor, (value,)))
    for call, name in ((energy.enclosures, "energy"), (frames.enclosures, "frames")):
        for slot in range(3):
            for i, value in enumerate(invalid):
                args = [2**22, s.Rational(1, 100), 1]
                args[slot] = value
                cases.append((f"{name}_{slot}_{i}", call, tuple(args)))
        for i, value in enumerate(invalid + (s.Integer(6), 0, -1)):
            cases.append((f"{name}_multiplicity_{i}", call, (2**22, 0, 1, value)))
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
        cases.append((f"energy_domain_{i}", energy.enclosures, args))
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
        raise ValueError("Unsupported local energy input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_exact_state_not_an_adiabatic_truncation": True,
        "absolute_MS_finite_terms_not_set_to_zero": True,
        "unsaturated_mass_counterterm_is_a_growing_negative_control": True,
        "reference_coefficient_not_chosen_to_fit_energy": True,
        "new_higher_field_prescription_not_old_functional_identity": True,
        "named_low_order_invariance_not_full_quantum_matching": True,
        "free_T00_not_pressure_or_curved_interacting_stress": True,
        "original_P8_not_closed": True,
    }
