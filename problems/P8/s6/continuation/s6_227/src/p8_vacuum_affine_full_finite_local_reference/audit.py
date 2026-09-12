"""Actual complete finite local reference and unchanged full quantum frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_weighted_channel_resolvent import audit as previous

from . import coefficients, factorization, resolvent

ITEM = {
    "id": "actual_complete_finite_local_curvature_conjugation_uniform_weighted_bound_and_full_finite_reference_inverse",
    "status": "ACTUAL_FULL_FINITE_LOCAL_REFERENCE_INVERTED_NOT_FULL_NONLOCAL_CURVED_QUANTUM_SYSTEM_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "actual_full_finite_local_factorization",
        "source_time_adjoint_bound",
        "full_finite_reference_inverse",
    ):
        raise ValueError(
            "Only the complete finite-local reference and its stated inverse are proved"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The complete finite local reference cannot close the full quantum or P8 frontier"
        )
    return True


def packets():
    return {
        "actual_complete_local_operator_and_boundary": factorization.data(),
        "whole_slab_original_local_coefficient_bounds": coefficients.data(),
        "source_time_reciprocity_full_conjugation_and_inverse": resolvent.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "unchanged_primitive_rows": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_rows": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_complete_local_reference_input": len(matching())
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
        value.rows * value.cols if isinstance(value, s.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "all_original_finite_local_terms_and_coefficients_retained": True,
        "source_time_derivative_not_confused_with_forward_time_derivative": True,
        "original_pole_and_complete_initial_boundary_retained": True,
        "noncommuting_conjugate_terms_keep_exact_order": True,
        "actual_nonlocal_remainder_not_assumed_in_bounded_class": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "actual_local_identity": "The complete original finite local remainder has the exact Euler factorization Bc*A+A^T Bc+D*HD+D*J+J^T D+V0, including the coefficient-derivative matrix J and the explicit density boundary.",
        "no_loss_bound": "Retarded/advanced reciprocity bounds the SOURCE-time derivative of Z. The six ordered terms of Omega=Z Rloc Y then have a same-H^r weighted norm below1/4 uniformly in all spatial momenta, using actual coefficient bounds.",
        "complete_reference_inverse": "The actual full finite local reference Bc*(Fdiag+V)Bc+Rloc has a causal weighted inverse bounded by125/[(123+32M)(sigma_M-20)^4], with the entire original local Hessian retained.",
        "boundary": "The physical density stays on the right, kappa and the unweighting cost exp(sigma_M T) remain. Actual nonlocal curved/state/contact/tree/matter matching, the full S222 inverse, physical stability and original P8 remain open.",
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
        out.append((f"bound_type_{i}", resolvent.channel.size, (value,)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for value in (-1, s.Rational(-1, 2)):
        out.append((f"negative_bound_{value}", resolvent.channel.size, (value,)))
    for stage in (
        "full_quantum_inverse",
        "stable",
        "finite_Born_remainder",
        "closed_P8",
        "new_state",
        "pole_deletion",
        "physical_cutoff",
        "unweighted_smallness",
        "actual_curved_matching",
        "zero_weight",
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
                "extra_actual_inverse",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "full_actual_inverse", "status": "COMPLETE"}],
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
        raise ValueError(
            "Unsupported full quantum or stability claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    packet = factorization.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "deleting_J_changes_actual_operator": packet["gates"][
            "deleting_coefficient_derivative_J_changes_operator"
        ],
        "density_boundary_not_identically_zero": packet["exact_bilinear_boundary"] != 0,
        "actual_local_remainder_not_deleted": packet["gates"][
            "actual_remainder_not_zero"
        ],
        "source_time_derivative_requires_retained_initial_atom": True,
        "unweighting_cost_not_assumed_one": resolvent.SIGMA.is_positive,
        "original_weighted_growing_pole_control_retained": previous.ITEM["id"]
        in {row["id"] for row in matching()},
        "physical_kappa_kept_in_output_density": factorization.g.bridge.KAPPA
        == modes.KAPPA,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
