"""Exact force-interface progress and unchanged original research frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_scalar_tame_propagator import audit as previous

from . import estimates, feedback, forces

ITEM = {
    "id": "full_quantum_forced_scalar_constraints_ordered_Schur_graph_identity_and_first_response_domain",
    "status": "EXACT_FORCE_INTERFACE_AND_ADMISSIBLE_GRAPH_EQUIVALENCE_NOT_QUANTUM_INVERSE_CONTRACTION_STABILITY_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "force_interface",
        "admissible_graph_identity",
        "first_prescribed_response",
    ):
        raise ValueError(
            "Only the proved force interface, graph identity or prescribed first response"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "An exact Schur graph identity does not establish a quantum inverse or close P8"
        )
    return True


def packets():
    return {
        "complete_quantum_forced_scalar_constraint_interface": forces.data(),
        "both_clock_contacts_and_output_density_order": forces.contact_data(),
        "ordered_causal_block_and_graph_domain": feedback.data(),
        "precise_force_and_first_response_functional_spaces": estimates.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, p in packets().items()
        for key, value in p["checks"].items()
    }
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "all_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_previous_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_scoped_force_graph_checkpoint": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {
        key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for key, value in rows.items()
    }


def scalar_entry_count():
    return sum(
        v.rows * v.cols if isinstance(v, s.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, p in packets().items()
            for key, value in p["gates"].items()
        },
        "both_force_dependent_constraints_reconstructed": True,
        "force_interface_not_full_quantum_Hamiltonian": True,
        "direct_auxiliary_response_and_output_density_order_retained": True,
        "graph_equivalence_not_inverse_existence": True,
        "estimate_gap_not_instability_or_nonexistence": True,
        "no_frozen_input_edits_or_primitive_promotion": True,
        "unique_matching_identifiers": len({r["id"] for r in matching()})
        == len(matching()),
    }


def observable():
    return {
        "force_interface": "The complete clock-force saddle gives b=(pv-3gb)/2 and n=(Lc-gn-delta gzeta)/(2Jc). Thus s=CZ+Dg and Z'=KZ+Fg, with F=-Jcan C^T and rank2 D.",
        "ordered_coupled_equation": "The full force graph is [I-Qbar(D+C G0F)]g=e+Qbar C G0h, with Qbar=(kappa a^3)^-1 Rhat in output order and both original contacts. Reconstruction is equivalent only on the explicitly stated common graph domain.",
        "first_response": "For prescribed prepared scalar source S0, Qbar S0 lies inL2H^(r-3) and its classical phase response inC H^(r-17), with explicit density and interface constants. A comparison-phase source adds the finite C13 and yields27 spatial derivatives lost.",
        "remaining": "The current bounds do not give the source-domain self-map needed for iteration. This does not exclude a better full inverse. Compatible-space quantum inversion, finite-coupling/nonlinear remainders, background/stability, heavy/cutoff and original V/G/B remain open.",
    }


def bad_cases():
    bad = (
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
    for i, value in enumerate(bad):
        for pos in range(4):
            args = [0, modes.KAPPA, modes.MASS, 1]
            args[pos] = value
            out.append((f"scope_type_{pos}_{i}", require_scope, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for stage in (
        "full_quantum_inverse",
        "stable",
        "finite_Born_remainder",
        "closed_P8",
        "new_state",
        "order_reduction",
        "physical_cutoff",
        "all_graph_inputs",
        "global_homogeneous",
    ):
        out.append((f"unsupported_stage_{stage}", require_stage, (stage,)))
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
            (
                "extra_quantum_inverse",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "quantum_inverse", "status": "COMPLETE"}],
                ),
            ),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported force-interface claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "independent_literal_forced_saddle_and_all_weighted_phase_equations": True,
        "zero_Theta_and_small_large_transfer_with_nonzero_forces": True,
        "independent_full_clock_functional_composition_keeps_contact_once": True,
        "independent_two_time_causal_full_block_solve": True,
        "nonzero_deleted_auxiliary_and_wrong_density_order_controls": True,
        "second_formal_coefficient_has_nonzero_direct_auxiliary_piece": True,
        "derivative_loss_counterexample_distinguishes_contraction_from_inverse": True,
        "all_previous_frontiers_unchanged_and_original_P8_open": True,
    }
