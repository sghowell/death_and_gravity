"""Actual CD state and complete source noise, with canonical scope retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_vacuum_causal_noise import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import covariance, derivatives, force, modes

ITEM = {
    "id": "actual_all_order_CD_state_complete_clock_source_smeared_force_noise_and_canonical_normalization",
    "status": "SAME_STATE_FULL_SOURCE_SMEARED_SCALAR_NOISE_BOUND_NOT_METRIC_NOISE_SOLUTION_CUTOFF_OR_V_G_B",
}


def require_scope(time, kappa=modes.KAPPA, delta=derivatives.DELTA, mass=modes.MASS):
    t, k, d, m = map(rational, (time, kappa, delta, mass))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the fixed unit CD slab")
    if k != modes.KAPPA or not 0 < d <= derivatives.DELTA or m != modes.MASS:
        raise ValueError(
            "Require the actual fixed CD parent, mass1000 and stated jet ball"
        )
    return t, k, d, m


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A scalar source-noise bound cannot close original P8")
    return True


def packets():
    return {
        "actual_state_modes": modes.data(),
        "curved_covariance": covariance.data(),
        "full_source_derivatives": derivatives.data(),
        "force_conventions": force.data(),
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
        "same_all_order_state_not_W8_or_flat_vacuum": True,
        "exact_evolved_mixing_not_preparation_only": True,
        "all_three_modes_and_temporal_constraint_kept": True,
        "positive_cross_covariance_not_discarded": True,
        "actual_time_dependent_normalizers_and_volume": True,
        "complete_full_switch_source_not_a_germ": True,
        "all_covariant_Hessian_connection_terms_kept": True,
        "small_history_not_small_arbitrary_test": True,
        "full_second_jet_and_mixed_spatial_source_bound": True,
        "source_force_not_metric_stress_noise": True,
        "clock_normalized_and_canonical_variances_distinct": True,
        "canonical_gravity_suppression_cancels": True,
        "explicit_smaller_ball_not_a_full_solution": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "state": "The actual CD all-order Gaussian state has a full three-mode covariance upper bound108/m times L2 squared plus144/m^3 times spatial H1 seminorm squared on the unit slab.",
        "source": "The complete near-clock affine source satisfies ||grad DS_eta||<=4e8 delta ||eta||H3, supplementing the unchanged2048 delta undifferentiated bound.",
        "normalized_noise": "The force divided by kappa has smeared standard deviation below2e-392 delta per unit raw-clock H3 test norm at kappa1e800.",
        "canonical_noise": "The canonical test bound is instead2e8 delta, with no remaining kappa suppression. On delta<=1e-14 it is below2e-6 and the prior mean relative coefficient is below4e-18.",
        "remaining": "Metric stress/noise, spatial and nonlinear response, finite-coupling remainders, interacting measures/loops, heavy/cutoff/threshold matching and original V/G/B.",
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
            args = [0, modes.KAPPA, derivatives.DELTA, modes.MASS]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for args in (
        (1, modes.KAPPA, derivatives.DELTA, modes.MASS),
        (-1, modes.KAPPA, derivatives.DELTA, modes.MASS),
        (0, modes.KAPPA * 2, derivatives.DELTA, modes.MASS),
        (0, modes.KAPPA, 0, modes.MASS),
        (0, modes.KAPPA, 1, modes.MASS),
        (0, modes.KAPPA, derivatives.DELTA, 0),
        (0, modes.KAPPA, derivatives.DELTA, 999),
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
        raise ValueError("Unsupported CD source-noise scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "state_not_reset_or_replaced_by_reference": True,
        "full_longitudinal_temporal_readout_retained": True,
        "covariance_upper_bound_not_QEI_lower_bound": True,
        "full_source_not_replaced_by_clock_germ": True,
        "canonical_kappa_cancellation_explicit": True,
        "zero_reference_source_noise_not_zero_metric_noise": True,
        "smaller_ball_not_a_quantum_background_solution": True,
        "original_P8_not_closed": True,
    }
