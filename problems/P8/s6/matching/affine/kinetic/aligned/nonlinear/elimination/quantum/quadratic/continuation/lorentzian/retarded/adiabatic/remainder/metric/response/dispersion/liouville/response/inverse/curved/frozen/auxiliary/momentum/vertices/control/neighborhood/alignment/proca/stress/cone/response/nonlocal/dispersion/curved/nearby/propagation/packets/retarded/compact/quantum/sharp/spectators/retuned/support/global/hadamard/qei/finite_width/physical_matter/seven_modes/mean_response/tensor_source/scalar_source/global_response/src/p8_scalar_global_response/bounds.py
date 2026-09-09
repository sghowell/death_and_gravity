"""Global scalar covariance, compensated mean and physical observable bounds."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import tails
from p8_scalar_mean_source import geometry, model, observable

from . import phase, source


@cache
def data():
    ph = phase.data()
    src = source.envelopes()
    u = model.u
    t = 1 + u * u
    Q = sp.Integer(3) ** 59
    Ccov = (40 * Q) ** 2
    If = (
        src["mean_scale_forcing"]["half_entry_L1_upper"]
        + src["weighted_mean_trace_forcing"]["half_entry_L1_upper"]
    )
    CY = 9 * If * Ccov
    nstate = (
        sum(phase.uniform_upper(row) for row in src["lapse_state_term"]["entries"]) / 2
    )
    CN = sp.Rational(9, 40) * CY + nstate * Ccov
    CB = CY + CN / 2
    CF = (CN / 10 + 3 * CY / 10) * tails.half_line_integral_upper(6) + src[
        "matter_field_state_term"
    ]["half_entry_L1_upper"] * Ccov
    B = geometry.canonical()["old_to_natural_phase"]
    Ti = ph["old_to_weighted_phase"].inv()
    physical = {}
    for key in ("intrinsic_density_Hessian", "intrinsic_pressure_Hessian"):
        matrix = (t**4 * Ti.T * B.T * observable.data()[key] * B * Ti).applyfunc(
            sp.factor
        )
        entries = []
        for i in range(4):
            for j in range(4):
                entries.append(
                    phase.envelope(
                        matrix[i, j],
                        int(i >= 2) + int(j >= 2),
                        require_integrable=False,
                    )
                )
        intrinsic = sum(phase.uniform_upper(row) for row in entries) / 2
        physical[key] = {
            "scaled_kernel_entry_envelopes": entries,
            "absolute_coefficient_per_t_minus_four": intrinsic * Ccov
            + sp.Rational(3, 100) * CB,
        }
    return {
        "phase_growth_integer_upper": Q,
        "improved_phase_covariance_entry_upper_per_eta": Ccov,
        "weighted_mean_forcing_L1_per_covariance_bound": If,
        "weighted_mean_phase_upper_per_eta": CY,
        "lapse_state_kernel_uniform_half_entry_sum": nstate,
        "global_lapse_upper_per_eta": CN,
        "global_linear_physical_scale_upper_per_eta": CB,
        "global_actual_frame_scale_upper_per_eta": CY + CN,
        "global_matter_field_upper_per_eta": CF,
        "global_physical_observable_bounds": physical,
        "global_original_trace_bound_per_eta": CY / t**3 + 36 * sp.Abs(u) * Q**2 / t,
        "lapse_state_all_kernel_entries_decay": all(
            2 * (row["power"] + sp.Rational(row["extra_half_order"], 2)) - row["parity"]
            > 0
            for row in src["lapse_state_term"]["entries"]
            if row["constant"] != 0
        ),
        "checks": {
            "global_weighted_phase_integral_below_59": ph["integrable_entry_sum_upper"]
            < 59,
            "natural_anchor_initial_norm_is_one": ph[
                "natural_anchor_initial_infinity_norm"
            ]
            == 1,
            "improved_momentum_bound_uses_retained_damping": all(
                1 + v < 40 for v in ph["momentum_remainder_row_bounds_per_one_over_t"]
            ),
        },
    }
