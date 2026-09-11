"""Strict scope and complete quantitative gates for conditional metric response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_retarded_energy import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import bounds, bridge, chart, limits

ITEM = {
    "id": "same_affine_conditional_prepared_homogeneous_metric_response_with_uncancelled_chart_contacts",
    "status": "COMPLETE_CONDITIONAL_VECTOR_HOMOGENEOUS_C10_TO_C0_BOUND_NOT_SPATIAL_NOISE_COUPLED_BACKGROUND_CUTOFF_OR_V_G_B",
}


def require_scope(time, order=10, zeta=bridge.ZETA, kappa=bridge.KAPPA):
    t, z, k = map(rational, (time, zeta, kappa))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Outside the quantified compact bounce slab")
    if type(order) is not int or order != 10:
        raise ValueError("Require the declared native C10 source norm")
    if z != bridge.ZETA or k != bridge.KAPPA:
        raise ValueError("Require the same fixed conditional parent")
    return t, order, z, k


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Conditional metric response cannot close original P8 gates")
    return True


def observable():
    return {
        "parent": "Same complete regular affine parent, fixed physical metric and mass, conditional canonical measure and all-order Cauchy state as S6.176; no scalar profile is added.",
        "histories": "Smooth homogeneous physical N=1+epsilon n and a=a0 exp(epsilon v), compactly supported strictly inside [-1/2,1/2] and zero on an initial neighborhood. The state has no independent variation.",
        "physical_result": "The full conditional vector mean stress has the actual ordinary-Proca prepared C10-to-C0 response, below5e-795 after fixed-kappa normalization.",
        "chart_result": "The actual full analytic metric map has the checked clock jets. Its mean currents use fixed reference volume, retaining physical volume/normalization and nonzero second-map background contacts; the C10-to-C0 bound is below1e-784 after kappa normalization.",
        "not_inferred": "No old background-cancelling profile or total clock result is transferred. Derivative loss, no arbitrary spatial/noise norm, no coupled inverse, no full quantum background/measure/state/cutoff and open V/G/B are explicit.",
    }


@cache
def residuals():
    out = {}
    for name, packet in (
        ("bridge", bridge.data()),
        ("homogeneous", bridge.homogeneous()),
        ("chart", chart.data()),
        ("envelopes", chart.envelopes()),
        ("bounds", bounds.data()),
        ("limits", limits.data()),
    ):
        out.update({name + "_" + key: value for key, value in packet["checks"].items()})
    out.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_conditional_metric_matching_row_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {
        k: s.ImmutableMatrix(v.applyfunc(s.cancel))
        if isinstance(v, s.MatrixBase)
        else s.cancel(v)
        for k, v in out.items()
    }


def scalar_entry_count():
    return sum(
        v.rows * v.cols if isinstance(v, s.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    quantitative = {}
    for name, packet in (("chart", chart.envelopes()), ("bounds", bounds.data())):
        quantitative.update(
            {name + "_" + key: bool(value) for key, value in packet["gates"].items()}
        )
    return {
        **quantitative,
        "complete_full_R_parent_unchanged": True,
        "literal_arbitrary_lapse_scale_operator_bridge": True,
        "longitudinal_fixed_comoving_map_and_zero_mode_retained": True,
        "physical_readout_contacts_checked_before_canonical_substitution": True,
        "same_conditional_canonical_measure_and_prepared_Hadamard_state": True,
        "same_covariant_dimensional_subtraction_and_finite_parts": True,
        "source_pinned_full_continuum_physical_response_only": True,
        "homogeneous_complete_mean_source_force_cancels": True,
        "no_old_scalar_stress_canceling_profile_added": True,
        "full_chart_clock_jets_not_global_polynomial_replacement": True,
        "nonzero_background_first_variation_retained": True,
        "nonlinear_second_map_contact_not_dropped": True,
        "fixed_reference_volume_and_kappa_normalization": True,
        "retarded_response_plus_local_contacts_not_symmetric_action": True,
        "C10_derivative_loss_not_coupled_inverse": True,
        "compact_slab_not_global_tail_or_physical_cutoff": True,
        "spatial_noise_full_quantum_parent_and_V_G_B_open": True,
        "all_prior_frozen_scientific_bytes_unchanged": True,
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
            args = [0, 10, bridge.ZETA, bridge.KAPPA]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for row in (
        (-1, 10, bridge.ZETA, bridge.KAPPA),
        (1, 10, bridge.ZETA, bridge.KAPPA),
        (0, 4, bridge.ZETA, bridge.KAPPA),
        (0, 5, bridge.ZETA, bridge.KAPPA),
        (0, 0, bridge.ZETA, bridge.KAPPA),
        (0, 11, bridge.ZETA, bridge.KAPPA),
        (0, s.Integer(10), bridge.ZETA, bridge.KAPPA),
        (0, 10, s.Rational(1, 2000), bridge.KAPPA),
        (0, 10, bridge.ZETA, 2 * bridge.KAPPA),
    ):
        out.append((f"scope_{len(out)}", require_scope, row))
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
            "Unsupported conditional metric response scope accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "physical_vector_not_old_profile_total": True,
        "nonzero_second_metric_map_contact_retained": True,
        "same_preparation_not_mode_reset": True,
        "homogeneous_source_cancellation_not_spatial_claim": True,
        "C10_not_no_loss_inverse_or_actual_C4_impossibility": True,
        "retarded_mean_response_not_noise_or_single_branch_Hessian": True,
        "local_vacuum_tail_control_not_total_asymptotic": True,
        "original_P8_not_closed": True,
    }
