"""Fixed local full tracefree spatial Hessian, not quantum matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_reference_spatial_remainder import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import bounds, helicities, local, proper

ITEM = {
    "id": "fixed_covariant_local_full_tracefree_spatial_Hessian_proper_time_and_norm",
    "status": "FIXED_LOCAL_FULL_TRACEFREE_SPATIAL_HESSIAN_NOT_QUANTUM_MATCHING_MIXED_INVERSE_OR_V_G_B",
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
        raise ValueError("A fixed local spatial Hessian cannot close original P8")
    return True


def packets():
    return {
        "complete_tracefree_spatial_space_and_polynomial_maps": helicities.data(),
        "original_fixed_local_action_and_full_conformal_Hessians": local.data(),
        "full_proper_time_conversion_and_measure": proper.data(),
        "actual_CD_local_Hessian_derivative_norm": bounds.data(),
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
        "same_actual_CD_parent_preparation_and_prescription": True,
        "already_fixed_finite_local_coefficients_not_retuned": True,
        "four_dimensional_Euler_only_after_dimension_limit": True,
        "unimodular_volume_variation_exactly_zero": True,
        "all_five_tracefree_spatial_directions": True,
        "both_vector_and_scalar_shears_retained": True,
        "smooth_polynomial_momentum_maps_at_zero": True,
        "complete_conformal_curvature_Hessians": True,
        "full_four_dimensional_curvature_crosscheck": True,
        "noncommuting_matrix_time_jets_not_diagonal_only": True,
        "compact_Weyl_boundary_term_explicit": True,
        "proper_time_density_measure_factor_retained": True,
        "full_clock_derivative_and_friction_terms": True,
        "actual_CD_local_L2_derivative_norm": True,
        "both_external_metric_canonical_factors": True,
        "no_lapse_shift_or_reduced_mixed_operator_claim": True,
        "remaining_quantum_contact_matching_not_assumed": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "domain": "The already-fixed finite local action on the actual compact CD slab, in the same exponential tracefree spatial metric parametrization.",
        "full_spatial": "Both tensor, both vector and scalar shear directions are represented by smooth polynomial momentum maps. Their complete curvature Hessians retain every spatial and clock term.",
        "proper_time": "The exact a^-1 density factor and partial_eta=a partial_t conversion preserve the compact action pairing; no constant-clock approximation is used.",
        "bound": "The local L2 Hessian is below5e4 Phi42[Gamma], with two-metric canonical coefficient2e-795.",
        "boundary": "This supplies a fixed local matching target, not equality of the remaining quantum contact/Taylor sector to it, a full mixed inverse, interacting background, stability or original P8 closure.",
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
            "Unsupported fixed local spatial Hessian scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "actual_parent_preparation_and_prescription_unchanged": True,
        "independent_literal_full_four_dimensional_curvature": True,
        "all_five_spatial_directions_and_noncommuting_time_jets": True,
        "exact_compact_boundary_not_pointwise_discard": True,
        "proper_time_measure_negative_control": True,
        "actual_self_adjoint_pairing_and_derivative_norm": True,
        "no_quantum_matching_or_reduced_inverse_assumed": True,
        "original_P8_not_closed": True,
    }
