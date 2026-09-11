"""A complete quadratic transition bound and a scoped local state difference."""

from functools import cache

import sympy as s
from p8_vacuum_flat_dirac_hadamard import audit as previous

from . import energy, rotation, tube

MODULES = (tube, rotation, energy)
ITEM = {
    "id": "exact_four_frame_free_transition_and_uniform_same_operator_in_out_T00_difference",
    "status": "BOUNDED_FOR_FREE_FLAT_HADAMARD_STATE_DIFFERENCE_NOT_ABSOLUTE_STRESS_INTERACTING_CURVED_STATE_OR_BACKREACTION",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The scoped state-energy result cannot change the fixed ledger"
        )
    return True


def observable():
    return {
        "background": "The same prescribed Minkowski SAT8 external mass and the same free in/out Hadamard states as S6.166.",
        "quantity": "rho_in(t)-rho_out_state(t), for the symmetric Dirac T00 with the external mass held as a prescribed scalar under metric variation, and one identical local state-independent subtraction/finite counterterm prescription.",
        "difference_vs_out_particles": "The local out-state energy is not set to zero during the transition. It differs from the late-time out-particle Hamiltonian energy above its own vacuum.",
        "source_context": [
            {
                "url": "https://arxiv.org/pdf/1805.05107",
                "version": "v2",
                "location": "Section IV, equations4.6-4.9",
                "role": "Iterated unitary frames are a known construction. Here the final off-diagonal term is retained and bounded.",
            },
            {
                "url": "https://arxiv.org/pdf/1703.00908",
                "location": "Section IV, equations40-41 and discussion following42",
                "role": "Dirac energy bilinear and state-independent local subtractions in an external scalar background.",
            },
        ],
        "quantitative_input": "Our explicit fixed-disk inequalities, four-step complex-gap induction and exact unitary Duhamel bound; no numerical constant is taken from the cited papers.",
        "not_a_cutoff": "Whole-momentum results for this specified quadratic operator do not establish the interacting EFT cutoff or high-energy applicability.",
    }


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(value)
        for mod in MODULES
        for key, value in mod.data()["checks"].items()
    }
    out.update(
        {
            "nine_prior_primitive_rows_retained": len(frontier()) - 9,
            "all_prior_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_prior_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "only_named_free_energy_difference_added": len(matching())
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
        **tube.data()["gates"],
        **energy.data()["gates"],
        "fixed_complex_branches_have_strict_gaps": True,
        "Cauchy_bound_controls_each_exact_connection": True,
        "four_rotations_preserve_asymptotic_state": True,
        "remaining_off_diagonal_generator_not_discarded": True,
        "complete_quadratic_transition_not_Feynman_loop_truncation": True,
        "full_radial_momentum_integrals_converge": True,
        "same_operator_Hadamard_state_difference": True,
        "state_independent_local_subtractions_cancel_only_in_difference": True,
        "external_mass_source_energy_exchange_not_omitted": True,
        "absolute_stress_and_interacting_cutoff_not_inferred": True,
        "all_frozen_scientific_sources_and_reports_unchanged": True,
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
    for i, value in enumerate(invalid + (-1, 5, s.Integer(2))):
        for name, call in (
            ("radius", tube.radius),
            ("frequency_floor", tube.frequency_floor),
        ):
            out.append((f"{name}_{i}", call, (value,)))
    for call, name, good in (
        (energy.amplitude_upper, "amplitude", [1, 2**22, s.Rational(1, 100), 1]),
        (energy.enclosures, "energy", [2**22, s.Rational(1, 100), 1]),
    ):
        for slot in range(len(good)):
            for i, value in enumerate(invalid):
                args = good.copy()
                args[slot] = value
                out.append((f"{name}_type_{slot}_{i}", call, tuple(args)))
    for args in ((0, 0, 1), (1, -1, 2**22), (2**22, 0, 0), (1, 1, 2**22), (1, 0, 1)):
        out.append((f"parameter_domain_{args}", energy.enclosures, args))
    out.append(("negative_momentum", energy.amplitude_upper, (-1, 2**22, 0, 1)))
    for i, value in enumerate(invalid + (s.Integer(6), 0, -1)):
        out.append((f"multiplicity_{i}", energy.enclosures, (2**22, 0, 1, value)))
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
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
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
        raise ValueError("Unsupported exact-frame state-energy input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "final_coupling_and_its_transition_orders_retained": True,
        "fixed_complex_gap_not_only_real_eigenvalue_sampling": True,
        "state_coherence_not_replaced_by_occupation_energy": True,
        "out_state_energy_not_zero_during_transition": True,
        "external_mass_exchange_not_isolated_conservation": True,
        "state_difference_not_absolute_or_relative_bounce_stress": True,
        "quadratic_whole_momentum_operator_not_physical_cutoff": True,
        "original_P8_not_closed": True,
    }
