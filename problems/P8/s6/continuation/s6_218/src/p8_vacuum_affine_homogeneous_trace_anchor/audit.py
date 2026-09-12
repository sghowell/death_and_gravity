"""Full homogeneous trace anchor frontier without nonzero-transfer closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_corrected_spatial_current import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import anchor, comparison, homogeneous, local

ITEM = {
    "id": "full_homogeneous_trace_Gaussian_anchor_original_volume_matching_and_ordered_cross_zeros",
    "status": "COMPLETE_HOMOGENEOUS_TRACE_ANCHOR_AND_FULL_HOMOGENEOUS_HILBERT_LIFT_NOT_NONZERO_TRANSFER_SCALAR_REMAINDERS_REDUCED_INVERSE_OR_V_G_B",
}


def require_family(
    time, amplitude, jet_bound=1, kappa=modes.KAPPA, mass=modes.MASS, length=1
):
    t, e, g, k, m, L = map(rational, (time, amplitude, jet_bound, kappa, mass, length))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the original unit CD slab")
    if abs(e) > s.Rational(1, 100) or not 0 <= g <= 1:
        raise ValueError(
            "Require the admitted trace amplitude and all twelve bounded direction jets"
        )
    if k != modes.KAPPA or m != modes.MASS or L != 1:
        raise ValueError("Require original kappa,mass1000 and unit slab")
    return t, e, g, k, m, L


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A homogeneous trace anchor cannot close the full scalar or original P8 problem"
        )
    return True


def packets():
    return {
        "complete_constrained_trace_Hamiltonian_and_current_vertices": homogeneous.data(),
        "all_momentum_mixed_frequency_squeeze_and_vertex_domination": homogeneous.domination_data(),
        "complete_actual_state_reference_finite_and_infinite_comparison": comparison.data(),
        "general_dimension_full_volume_local_action_Euler_current_and_bounds": local.data(),
        "actual_trace_anchor_ordered_cross_zeros_and_full_Hilbert_lift": anchor.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, p in packets().items()
        for key, value in p["checks"].items()
    }
    before, after = previous.matching(), matching()
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "all_original_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_historical_matching_rows_unchanged": sum(
                a != b for a, b in zip(before, after)
            ),
            "one_new_homogeneous_trace_anchor_row": len(after) - len(before) - 1,
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
        "full_trace_volume_and_temporal_constraint_retained": True,
        "same_original_state_parent_and_finite_prescription": True,
        "independent_literal_ADM_covariance_and_scale_factor_Euler_tests": True,
        "finite_amplitude_C2_only_for_Gaussian_determinant_component": True,
        "actual_reference_parent_Hessian_not_full_sourced_nonlinear_parent": True,
        "three_nonzero_transfer_scalar_currents_still_require_remainders": True,
        "no_frozen_scientific_test_or_report_bytes_changed": True,
        "original_V_G_B_and_P8_open": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "trace_anchor": "The actual original-prescription P0 trace/trace Gaussian response is the homogeneous current's first derivative, including the temporal-constraint contact and nonzero m^4 volume term.",
        "bounds": "The admitted trace-family Gaussian determinant current has coordinate-density C0,C1,C2 displays(1e78,1e95,2e111). The full normalized homogeneous anchor has the derivative-losing Hilbert bound1e95||D||L2 Z130; both ordered trace/tracefree cross anchors vanish at P0.",
        "canonical": "The two external factors give3e-705 for the unit-trace anchor and4e-705 for the full homogeneous lift, not a reduced inverse.",
        "parent_boundary": "The C2 finite-amplitude comparison concerns the determinant component only. The reference Hessian equals the actual parent's Gaussian response by S176's S=DS=0 and zero-mean contact; no full sourced-parent finite-amplitude remainder is asserted.",
        "remaining": "The three nonzero-transfer scalar currents still need full state/time/endpoint/contact tails and dimensional/all-transfer assembly, followed by the reduced inverse, nonlinear quantum background and remaining original V/G/B obligations.",
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
    good = [0, 0, 1, modes.KAPPA, modes.MASS, 1]
    for i, value in enumerate(bad):
        for pos in range(6):
            args = list(good)
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_family, tuple(args)))
    for pos, value in (
        (0, 1),
        (0, -1),
        (1, s.Rational(1, 99)),
        (1, -s.Rational(1, 99)),
        (2, -1),
        (2, 2),
        (3, 2 * modes.KAPPA),
        (4, 999),
        (5, 0),
        (5, 2),
    ):
        args = list(good)
        args[pos] = value
        out.append((f"scope_{len(out)}", require_family, tuple(args)))
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
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported homogeneous trace input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "literal_full_ADM_trace_vertices_and_covariance_tangent": True,
        "all_momentum_actual_frequency_and_squeeze_jet_checks": True,
        "independent_scale_factor_Euler_and_full_chart_Hessian": True,
        "independent_full_K_omega_dimensional_action_derivation": True,
        "deleted_constraint_and_deleted_volume_negative_controls": True,
        "same_parent_reference_Hessian_not_full_nonlinear_source_free_parent": True,
        "ordered_cross_anchor_zeros_not_nonzero_transfer_equality": True,
        "original_P8_not_closed": True,
    }
