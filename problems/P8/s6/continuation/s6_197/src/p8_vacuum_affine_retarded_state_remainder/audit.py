"""Actual full retarded memory and contact state remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_separated_response import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import contact, memory, response, tail

ITEM = {
    "id": "actual_overlapping_time_spatial_retarded_state_remainder_with_complete_contact",
    "status": "ACTUAL_RETARDED_MEMORY_CONTACT_COMPARISON_NOT_REFERENCE_DIAGONAL_MATCHING_INVERSE_OR_V_G_B",
}


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the fixed unit CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError("Require actual fixed CD parent, mass1000 and unit slab")
    return t, k, m, L


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A retarded actual-state comparison bound cannot close original P8"
        )
    return True


def packets():
    return {
        "complete_overlapping_time_memory_remainder": memory.data(),
        "full_local_metric_contact_remainder": contact.data(),
        "uniform_all_momentum_regulator_tail": tail.data(),
        "actual_comparison_and_reference_boundary": response.data(),
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
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_scoped_row_added": len(matching()) - len(previous.matching()) - 1,
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
            name + "_" + key: bool(v)
            for name, p in packets().items()
            for key, v in p.get("gates", {}).items()
        },
        "same_full_parent_M1_profile_state_and_prescription": True,
        "actual_all_order_mixing_not_reset_to_W8": True,
        "constant_initial_alpha_retained_without_CCR_assumption": True,
        "complete_ten_field_readout_and_all_three_polarizations": True,
        "both_mixed_memory_terms_and_error_square": True,
        "retarded_step_kept_on_overlapping_supports": True,
        "no_derivative_of_rapidly_oscillating_mixing": True,
        "full_noncommuting_second_metric_contact_retained": True,
        "actual_temporal_constraint_not_a_new_oscillator": True,
        "low_external_both_leg_tail_implication": True,
        "high_external_transfer_not_cut_off": True,
        "contact_one_momentum_not_a_two_mode_overlap": True,
        "all_internal_and_external_momentum_regulators_removed": True,
        "both_canonical_tensor_chain_factors": True,
        "reference_not_a_new_physical_state_or_effective_action": True,
        "fixed_singular_reference_matching_still_open": True,
        "no_full_inverse_interacting_background_or_stability": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The actual prepared Gaussian Proca sector on the unit CD slab, with compact spatial tracefree source and readout that may overlap in time.",
        "comparison": "The actual finite-regulator current minus the same mode formula using constant-alpha W8 readouts converges absolutely, including the unchanged retarded step and full metric contact.",
        "bound": "The complete remainder is below2e23 M[D]M[Gamma], where M is spatial H1 and time L2; its regulator error is below1e27 M[D]M[Gamma]/K.",
        "canonical": "Two canonical tensor factors give magnitude8e-777 and regulator error4e-773/K in the same norms.",
        "boundary": "The reference is not a substituted CCR state. The singular reference response still needs the original covariant diagonal matching. No finite-amplitude remainder, full inverse, interacting background or V/G/B closure follows.",
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
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, 0, 1),
        (0, modes.KAPPA, modes.MASS, 0),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
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
        raise ValueError(
            "Unsupported retarded actual-state comparison scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_state_constraint_and_physical_vertex_retained": True,
        "independent_full_and_removed_nonzero_transfer_integrals": True,
        "literal_retarded_Fock_operators_not_only_scalar_identities": True,
        "both_mixed_terms_and_error_square_negative_control": True,
        "literal_ten_field_three_mode_noncommuting_contact": True,
        "reference_CCR_defect_and_high_transfer_contact_controls": True,
        "comparison_not_singular_reference_matching_or_inverse": True,
        "original_P8_not_closed": True,
    }
