"""All-momentum isolated quotient inverse with the full research frontier intact."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_isolated_shear_resolvent import audit as previous

from . import estimates, geometry, resolvent

ITEM = {
    "id": "full_flat_scalar_gauge_quotient_two_channel_causal_reference_inverse_uniform_all_momentum_bound",
    "status": "FLAT_REFERENCE_QUOTIENT_INVERSE_NO_SPATIAL_LOSS_NOT_FULL_CURVED_COUPLED_INVERSE_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "flat_gauge_quotient",
        "two_channel_reference_inverse",
        "uniform_spatial_bound",
    ):
        raise ValueError(
            "Only the flat quotient reference and its uniform inverse bound are proved"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A flat quotient reference inverse cannot close the curved or original P8 frontier"
        )
    return True


def packets():
    return {
        "complete_covariant_scalar_quotient": geometry.data(),
        "two_channel_causal_reference_inverse": resolvent.data(),
        "uniform_time_spatial_bounds_and_units": estimates.data(),
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
            "unchanged_primitive_rows": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_rows": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_isolated_quotient_input": len(matching())
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
        "original_finite_and_physical_channel_normalizations_retained": True,
        "both_shear_pole_and_continuum_retained": True,
        "uniform_bound_not_uniform_half_line_scalar_L1": True,
        "flat_quotient_not_full_three_source_or_curved_inverse": True,
        "no_kappa_smallness_or_primitive_promotion": True,
        "unique_matching_identifiers": len({r["id"] for r in matching()})
        == len(matching()),
    }


def observable():
    return {
        "full_quotient": "The full spatial covariant projectors and original finite Hessian give B^T diag(Ftrace(D^2+q),(8/3)F2(D^2+q)) B. Both ordered full-metric gauge legs reduce to this quotient; the three-source reference still has a gauge kernel.",
        "actual_inverse": "The causal coordinate kernel Rq has no inverse-q singularity atq0 and uniformly bounded first derivative. The shifted positive spectral measures give a diagonal primitive bounded by45/2. Their ordinary convolution Eq=Rq'*Jq*Rq^T is the complete quotient inverse.",
        "uniform_result": "||Eq(t)||<=375t^3/16, hence the normalized reference inverse is bounded by375T^4/64<6T^4 onC_tH^r andL2_tH^r for every realr, with no spatial derivative loss. The first time derivative is also controlled.",
        "normalization_and_boundary": "The physical Hessian inverse retains64pi^2; the S222 force-normalized reference also retainskappa, giving375pi^2 kappa T^4. This is not a small full feedback bound. The complete curved/state/contact/tree/matter remainder, common graph bridge and original V/G/B remain open.",
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
        out.append((f"transfer_type_{i}", resolvent.transfer, (value,)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for value in (-1, s.Rational(-1, 2)):
        out.append((f"negative_transfer_{value}", resolvent.transfer, (value,)))
    for stage in (
        "full_quantum_inverse",
        "stable",
        "finite_Born_remainder",
        "closed_P8",
        "new_state",
        "pole_deletion",
        "physical_cutoff",
        "half_line_L1",
        "full_three_source_inverse",
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
        raise ValueError("Unsupported full quotient claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    M = geometry.M
    wrong = M * s.Matrix([-geometry.lam, 0, geometry.q])
    metric = geometry.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "wrong_source_time_gauge_sign_detected": wrong != s.zeros(2, 1),
        "full_three_source_gauge_block_not_invertible": metric[
            "full_three_source_reference"
        ]
        .det()
        .expand()
        == 0,
        "dropping_shear_channel_loses_quotient_rank": s.factor(
            metric["two_channel_quotient_reference"].det().subs(geometry.f2, 0)
        )
        == 0,
        "full_ordered_three_source_not_naively_symmetric": s.factor(
            metric["full_three_source_reference"][0, 2]
            - metric["full_three_source_reference"][2, 0]
        )
        != 0,
        "q0_coordinate_extension_not_division_by_q": resolvent.coordinate_kernel(0)[
            0, 0
        ]
        == resolvent.t,
        "retained_isolated_shear_pole": previous.ITEM["id"]
        in {r["id"] for r in matching()},
        "physical_force_inverse_bound_not_kappa_small": 375 * s.pi**2 * modes.KAPPA > 1,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
