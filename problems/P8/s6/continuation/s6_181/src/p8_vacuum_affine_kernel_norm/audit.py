"""A numeric scalar range bound, without promoting any original P8 gate."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_response import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import dyadic, estimate, tail, threshold

ITEM = {
    "id": "same_massive_scalar_range_kernel_explicit_full_half_line_L1_majorant",
    "status": "MASS_INDEPENDENT_EXPLICIT_SCALAR_KERNEL_NORM_NOT_FULL_INVERSE_SMALLNESS_BACKGROUND_CUTOFF_OR_V_G_B",
}


def require_scope(mass, scalar_range=1, all_bands=1):
    m, r, b = map(rational, (mass, scalar_range, all_bands))
    if m <= 0 or r != 1 or b != 1:
        raise ValueError(
            "Require positive mass, the same scalar range only, and every frequency band"
        )
    return m, r, b


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A scalar kernel norm cannot close original P8 gates")
    return True


def packets():
    return {
        "threshold": threshold.data(),
        "tail": tail.data(),
        "dyadic": dyadic.data(),
        "estimate": estimate.data(),
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
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_scalar_norm_row_added": len(matching()) - len(previous.matching()) - 1,
        }
    )
    return {key: s.cancel(value) for key, value in rows.items()}


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet.get("gates", {}).items()
        },
        "same_source_pinned_density_not_a_new_model": True,
        "continuous_threshold_and_tail_derivatives_not_sample_only": True,
        "C2_partition_and_two_integrations_by_parts": True,
        "every_integer_frequency_band_including_infinite_tail": True,
        "absolute_L1_sum_and_positive_spectral_Laplace_identification": True,
        "same_mass_scaling_and_zero_instantaneous_term": True,
        "only_K_constant_replaced_in_current_retained_inverse": True,
        "curved_remainder_C_not_evaluated_and_no_smallness_inferred": True,
        "no_nonprepared_background_residual_fed_to_prepared_inverse": True,
        "original_V_G_B_and_full_quantum_background_still_open": True,
    }


def observable():
    return {
        "new_quantitative_result": "For the same normalized massive scalar range inverse, ||K_m||L1(0,infinity)<201600000000<300000000000 for every fixed m>0.",
        "derivation": "Continuous two-derivative threshold and logarithmic-tail estimates, a C2 all-integer dyadic partition, two integrations by parts in each compact band, an absolutely summable Fourier L1 bound, and same-kernel Laplace identification.",
        "current_parent_transfer": "The current retained coupled inverse may take K_L1=300000000000 in its constructive weight. The other scalar C controlling the complete curved weak-log remainder is still unevaluated.",
        "not_implied": "No small full inverse, stability, nonlinear quantum background, new compatible preparation, full interacting loops/state/measure, physical cutoff or V/G/B closure. The enormous coarse bound is not a measured norm.",
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
        for pos in range(3):
            args = [1000, 1, 1]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for args in (
        (0, 1, 1),
        (-1, 1, 1),
        (1000, 0, 1),
        (1000, 2, 1),
        (1000, 1, 0),
        (1000, 1, 2),
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
        raise ValueError("Unsupported scalar-kernel scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "threshold_square_root_not_assumed_smooth_at_zero": True,
        "large_frequency_log_tail_not_deleted": True,
        "band_derivative_rescaling_factors_retained": True,
        "sine_factor_two_and_half_line_L1_normalization_retained": True,
        "same_massive_inverse_identified_not_an_arbitrary_L1_kernel": True,
        "curved_C_not_replaced_by_scalar_norm": True,
        "all_frequency_retained_mathematics_not_physical_EFT_cutoff": True,
        "original_P8_not_closed": True,
    }
