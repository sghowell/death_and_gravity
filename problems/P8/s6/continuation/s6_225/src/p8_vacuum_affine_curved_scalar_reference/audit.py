"""Curved geometry, ordered reference inversion, and unpromoted frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_flat_scalar_quotient_inverse import audit as previous

from . import geometry, local, resolvent

ITEM = {
    "id": "actual_curved_scalar_curvature_coordinates_complete_finite_local_Hessian_and_ordered_reference_inverse",
    "status": "CURVATURE_ADAPTED_REFERENCE_UNIFORM_INVERSE_NOT_ACTUAL_FULL_CURVED_QUANTUM_INVERSE_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "actual_curvature_coordinates",
        "complete_finite_local_Hessian",
        "ordered_curved_reference_inverse",
    ):
        raise ValueError(
            "Only the actual local geometry and specified reference inverse are proved"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A curvature-adapted reference cannot close the actual quantum or P8 frontier"
        )
    return True


def packets():
    return {
        "actual_curvature_coordinates": geometry.data(),
        "complete_original_finite_local_Hessian": local.data(),
        "ordered_curvature_adapted_reference_inverse": resolvent.data(),
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
            "one_new_curvature_reference_input": len(matching())
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
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "original_finite_coefficients_and_all_background_local_terms_retained": True,
        "actual_formal_adjoint_coefficient_derivative_and_factor_order_retained": True,
        "fixed_flat_scalar_factors_not_asserted_actual_curved_mass_or_state": True,
        "uniform_reference_bound_keeps_pole_and_complete_initial_boundary": True,
        "no_full_quantum_inverse_kappa_smallness_or_primitive_promotion": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "actual_geometry": "In fixed conformal time the actual exponential spatial scalar chart gives S=(D^2+3hD+2q/3)w-(D^2+3hD)c/3 and W=-qw-D^2c. Their matrix is Bc=B0+C D with ||C||<8 and ||C'||<35 on the original slab.",
        "complete_local_Hessian": "The unchanged finite R_old^2, Weyl, Einstein and volume Hessians are retained completely. The curvature-square factor is -4SD SG-(4/45)WD WG and the explicit remainder has weighted total differential order at most2, not an asserted small compatible-space norm. Its homogeneous unit-Frobenius TF restriction agrees with S189.",
        "ordered_reference_inverse": "The causal inverses Y=Bc^-1 and Z=(Bc*)^-1 retain Bc*=B0^T-D(C^T), including Cprime. The curvature-adapted reference inverse is the ordered Y Kdiag Z, represented by Y Jdiag partial Z. Its norm is at most375T^4 exp(108T)/64 on C_tH^r and L2_tH^r, all momenta and every realr.",
        "boundary": "The physical conformal-density reference inverse retains right multiplication by64pi^2 kappa a^4 and costs at most2250pi^2 kappa T^4 exp(108T). The scalar factors are still the specified flat reference factors. Actual curved mass/state/contact/tree/matter differences and the full quantum-force graph remain open.",
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
        "pole_deletion",
        "physical_cutoff",
        "half_line_L1",
        "commuting_adjoint",
        "small_local_remainder",
        "actual_curved_scalar_factors",
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
        raise ValueError("Unsupported actual curved quantum claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    geom = geometry
    detector = s.Matrix([geom.wd, geom.cd])
    return {
        "rejected_inputs": rejected_inputs(),
        "deleting_coefficient_derivative_changes_actual_adjoint": geom.C.diff(
            geom.eta
        ).T
        * detector
        != s.zeros(2, 1),
        "flat_coordinate_operator_misses_actual_Hubble_terms": geom.C != s.zeros(2),
        "deleting_lower_order_local_remainder_changes_Hessian": local.data()[
            "complete_lower_order_local_remainder"
        ]
        != 0,
        "finite_order_two_not_a_zero_local_remainder": local.data()[
            "weighted_derivative_orders_in_remainder"
        ]
        == (0, 1, 2),
        "physical_conformal_density_not_constant": s.diff(
            geom.profile()["a"] ** 4, geom.profile()["t"]
        )
        != 0,
        "retained_flat_quotient_and_shear_pole_inputs": previous.ITEM["id"]
        in {row["id"] for row in matching()},
        "physical_inverse_bound_not_kappa_small": 2250 * s.pi**2 * modes.KAPPA > 1,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
