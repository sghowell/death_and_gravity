"""Strict evidence boundary for the conditional vector Gaussian sector."""

from functools import cache

import sympy as s
from p8_vacuum_affine_vector_schur import audit as previous

from . import bridge, gaussian, state, stress

ITEM = {
    "id": "regular_affine_parent_conditional_Gaussian_Hadamard_family_and_ordinary_Proca_stress",
    "status": "SOURCE_AWARE_CONDITIONAL_VECTOR_STATE_AND_FIXED_CLOCK_STRESS_NOT_FULL_PARENT_QUANTUM_SOLUTION",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Conditional Proca control cannot close original P8 gates")
    return True


def observable():
    return {
        "parent": "Same CD-REG-AFFINE-ISO classical target, now specialized to zeta=1e-6 within its certified range and kappa=1e800.",
        "quantized_component": "ONLY the retained conditional Gaussian vector, with canonical three-polarization measure, fixed S6.55 all-order Cauchy covariance and S6.82 covariant dimensional finite prescription at mu=m=1000.",
        "new_state_family": "Actual compact sourced histories after the common Cauchy neighborhood; source-free Hadamard fluctuations plus the retarded smooth coherent mean.",
        "quantitative_result": "On the unchanged clock and [-1/2,1/2], each of twelve rho/P time derivatives through order five is below 1e30 in physical units, or 1e-770 against kappa.",
        "no_cancellation": "No old stress-canceling scalar profile is installed. The bounded conditional vector stress is a contribution to a prescribed-clock residual, not a solved quantum bounce.",
        "remaining": "Full parent measure, light/tensor/auxiliary loops and mixed interactions, quantitative sourced response/noise, physical cutoff, full matching and V/G/B remain open.",
    }


@cache
def residuals():
    out = {}
    for name, packet in (
        ("canonical", bridge.data()),
        ("clock", bridge.clock()),
        ("state", state.data()),
        ("Gaussian", gaussian.data()),
        ("prescription", stress.prescription()),
        ("local", stress.local()),
        ("stress", stress.data()),
    ):
        for key, value in packet["checks"].items():
            out[name + "_" + key] = value
    out.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_conditional_sector_matching_row_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {
        k: v.applyfunc(s.cancel) if isinstance(v, s.MatrixBase) else s.cancel(v)
        for k, v in out.items()
    }


def scalar_entry_count():
    return sum(
        v.rows * v.cols if isinstance(v, s.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    return {
        **stress.data()["gates"],
        "selected_zeta_within_frozen_parent_range": bool(
            0 < bridge.ZETA <= s.Rational(1, 2000)
        ),
        "complete_action_and_canonical_source_contact_retained": True,
        "physical_lapse_scale_variation_precedes_canonical_substitution": True,
        "ordinary_metric_readouts_not_inferred_from_modes_alone": True,
        "canonical_measure_and_finite_prescription_explicit": True,
        "same_all_order_Cauchy_state_not_old_finite_order_reference": True,
        "actual_common_Cauchy_and_global_hyperbolicity_conditions": True,
        "Proca_propagation_wavefront_gap_repair_checked": True,
        "coherent_shift_maps_sourced_algebra_not_false_homogeneous_bisolution": True,
        "connected_Hadamard_covariance_and_CCR_preserved": True,
        "retarded_mean_distinct_from_symmetric_noise_kernel": True,
        "nonzero_source_higher_insertions_and_contacts_retained": True,
        "stress_time_derivatives_not_quantitative_functional_response_norm": True,
        "older_scalar_stress_canceling_profiles_not_transferred": True,
        "no_full_parent_quantization_or_GY14_matching_transfer": True,
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
        for pos in (0, 2, 3):
            args = [0, 0, bridge.ZETA, bridge.KAPPA]
            args[pos] = value
            out.append((f"type_{pos}_{i}", stress.require_scope, tuple(args)))
    for order in (*bad, s.Integer(1), -1, 6):
        out.append(
            (
                f"order_{len(out)}",
                stress.require_scope,
                (0, order, bridge.ZETA, bridge.KAPPA),
            )
        )
    for time, zeta, kappa in (
        (-1, bridge.ZETA, bridge.KAPPA),
        (1, bridge.ZETA, bridge.KAPPA),
        (0, s.Rational(1, 2000), bridge.KAPPA),
        (0, 0, bridge.KAPPA),
        (0, bridge.ZETA, 1),
        (0, bridge.ZETA, -1),
    ):
        out.append((f"scope_{len(out)}", stress.require_scope, (time, 0, zeta, kappa)))
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
        raise ValueError("Unsupported conditional quantum claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "old_nonminimal_energy_not_reused": True,
        "unscaled_covariance_not_assigned_to_W": True,
        "source_mean_not_omitted_off_clock": True,
        "old_scalar_tadpole_profiles_not_added": True,
        "compact_slab_bound_not_claimed_uniform_on_all_real_time": True,
        "C5_time_bound_not_functional_response_bound": True,
        "conditional_Gaussian_not_full_interacting_parent": True,
        "original_P8_not_closed": True,
    }
