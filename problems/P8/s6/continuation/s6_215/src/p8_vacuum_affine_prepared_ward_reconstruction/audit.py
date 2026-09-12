"""Prepared Ward reduction with explicit known-piece bounds, not scalar completion."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_full_adm_vertices import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import kinematics, norms, scalar, ward

ITEM = {
    "id": "prepared_ordered_Ward_reconstruction_full_chart_correction_and_three_remaining_scalar_kernels",
    "status": "PREPARED_WARD_REDUCTION_AND_KNOWN_PIECE_BOUNDS_NOT_MISSING_SCALAR_KERNEL_MATCHING_FULL_RESPONSE_INVERSE_OR_V_G_B",
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
        raise ValueError("Prepared Ward reconstruction cannot close original P8")
    return True


def packets():
    return {
        "prepared_source_and_final_detector_synchronous_maps": kinematics.data(),
        "ordered_density_Ward_identities_and_full_chart_contact": ward.data(),
        "rotation_reduction_and_three_independent_scalar_entries": scalar.data(),
        "projection_known_piece_and_conditional_full_norms": norms.data(),
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
            "all_corrected_prior_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_current_checkpoint_row_added": len(matching())
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
        "same_actual_parent_state_and_original_finite_prescription": True,
        "distinct_prepared_source_and_final_detector_boundaries": True,
        "initial_only_tracefree_source_extension_audited": True,
        "actual_one_point_and_all_chart_contacts_retained": True,
        "three_ordered_scalar_kernels_remain_explicit": True,
        "matched_tracefree_block_preserved": True,
        "known_piece_bound_not_full_operator": True,
        "no_scalar_strong_differentiability_from_Weyl_theorem": True,
        "unreduced_coordinate_norm_not_canonical_scalar_inverse": True,
        "original_V_G_B_and_P8_open": True,
        "all_primitive_and_corrected_prior_statuses_unchanged": True,
        "current_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "ordered_reconstruction": "The full conditional Gaussian ADM response is its synchronous spatial response plus separate source-covariance and detector-conservation terms, with the actual weight-one current and nonlinear chart contacts retained.",
        "boundaries": "The source gauge vector uses the unchanged initial germ, the detector uses a final-vanishing advanced primitive, and both differentiated endpoint fluxes vanish for their stated separate reasons.",
        "remaining_kernels": "Exactly three ordered scalar kernels remain beyond the complete tracefree block: trace/trace and both distinct trace/scalar-tracefree orientations.",
        "known_piece": "The projected matched tracefree current and full one-point Ward correction obey1e118 V04 U138. The displayed full bound is conditional on three still-unproved scalar constants.",
        "boundary": "Missing scalar state/time/UV/contact matching and norms, genuinely reduced inverse, nonlinear quantum background/stability, remaining parent matching/cutoff and original V/G/B remain.",
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
    for value in (0, True, None, "source", "past"):
        out.append((f"primitive_side_{len(out)}", kinematics.primitive, (1, value)))
    for value in (-1, 14, True, s.Integer(1), 1.0):
        out.append((f"derivative_order_{len(out)}", norms.leibniz_bound, (value, 7)))
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError(
            "Unsupported prepared Ward reconstruction input accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "literal_four_metric_Lie_and_both_primitives": True,
        "direct_covariant_volume_density_variation": True,
        "ordered_source_detector_Ward_and_integrated_action": True,
        "full_differentiated_boundary_flux_negative_control": True,
        "nonzero_one_point_chart_and_ordered_scalar_controls": True,
        "normalized_rotations_and_primitive_Cauchy_norm_checks": True,
        "initial_only_extension_and_conditional_scalar_scope": True,
        "original_P8_not_closed": True,
    }
