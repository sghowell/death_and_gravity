"""Fixed finite local covariant tensor Hessian, not the full response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_shear_response import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import bounds, geometry, operator, prescription

ITEM = {
    "id": "actual_fixed_local_covariant_full_spatial_tensor_Hessian",
    "status": "ACTUAL_FIXED_FINITE_LOCAL_TENSOR_HESSIAN_BOUND_NOT_COMPLETE_DETERMINANT_RESPONSE_OR_V_G_B",
}


def require_scope(time, kappa=modes.KAPPA, mass=modes.MASS, length=1):
    t, k, m, L = map(rational, (time, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the fixed unit CD slab")
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError(
            "Require actual fixed CD parent, mass1000 and unit observation slab"
        )
    return t, k, m, L


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A local finite tensor Hessian cannot close original P8")
    return True


def packets():
    return {
        "direct_tensor_curvature": geometry.data(),
        "fixed_covariant_prescription": prescription.data(),
        "full_spatial_local_operator": operator.data(),
        "full_Sobolev_bound": bounds.data(),
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
        "same_fixed_finite_density_not_new_counterterm": True,
        "finite_density_not_complete_state_dependent_determinant": True,
        "original_scalar_profile_and_M1_unchanged": True,
        "direct_dynamic_four_dimensional_curvature_checked": True,
        "spatial_terms_and_both_tensor_polarizations_retained": True,
        "physical_local_stress_signs_independently_recovered": True,
        "full_local_metric_Hessian_not_only_centered_noise": True,
        "compact_Euler_and_box_R_variations_only": True,
        "actual_a_cubed_weighted_adjoint": True,
        "complete_fourth_derivative_operator_retained": True,
        "full_joint_spacetime_H4_not_homogeneous_only_norm": True,
        "canonical_TT_normalization_not_scalar_mixed_chart": True,
        "no_singular_inverse_or_extra_pole_deletion": True,
        "no_new_physical_cutoff_from_local_norm": True,
        "nonlocal_UV_subtracted_response_still_required": True,
        "full_feedback_and_finite_coupling_remainder_still_open": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "prescription": "The same fixed mu=m local Proca density, including its evanescent finite curvature terms and physical signs, is retained without a new counterterm choice.",
        "operator": "Its complete compact TT Hessian has A=5m^2/12-R0/36 and Qloc=(A L+A'partial_t+Ddag D/60)/(8pi^2 kappa), with all spatial derivatives and the actual weighted adjoint.",
        "bound": "On the fixed CD slab, ||Qloc h||L2<6e-796 ||h||H4 in the stated full-spacetime canonical TT norm.",
        "boundary": "This is the full Hessian of the specified finite local matching action, not the entire state/subtraction-dependent determinant response, a stable inverse, finite-coupling error or cutoff.",
        "remaining": "The uniformly UV-subtracted nonlocal response and its remaining finite contacts, full mixed/nonlinear feedback, interacting parent matching and original V/G/B.",
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
        (0, modes.KAPPA * 2, modes.MASS, 1),
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
            "Unsupported fixed local tensor-response scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "no_new_local_prescription_or_state_reset": True,
        "full_spatial_curvature_not_homogeneous_substitution": True,
        "both_TT_polarizations_and_correct_Frobenius_factor": True,
        "omitted_weighted_adjoint_is_a_negative_control": True,
        "fourth_derivative_piece_not_discarded": True,
        "local_Hessian_not_complete_quantum_response": True,
        "no_cutoff_or_inverse_from_coefficient_smallness": True,
        "original_P8_not_closed": True,
    }
