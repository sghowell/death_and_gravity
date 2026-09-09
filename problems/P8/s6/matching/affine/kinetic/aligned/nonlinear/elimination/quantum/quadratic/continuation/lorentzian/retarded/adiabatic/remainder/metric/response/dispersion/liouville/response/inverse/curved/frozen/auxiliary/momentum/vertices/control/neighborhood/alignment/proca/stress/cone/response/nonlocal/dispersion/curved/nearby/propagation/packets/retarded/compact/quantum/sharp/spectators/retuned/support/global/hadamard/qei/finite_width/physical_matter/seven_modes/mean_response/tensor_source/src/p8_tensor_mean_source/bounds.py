"""Independent exact global response and integrated clock-exchange bounds."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import mean as prior_mean
from p8_proca_mean_response import tails as prior_tails

from . import operator, response


@cache
def data():
    u = prior_mean.u
    t = 1 + u * u
    inherited = prior_tails.data()
    C = sp.Integer(6400)
    L = sp.Integer(1500)
    r0 = sp.Integer(4)
    B = sp.Integer(7200)
    lapse = sp.Rational(9, 40) * C + sp.Rational(25, 4) * r0
    field = (L / 10 + 3 * C / 10) * prior_tails.half_line_integral_upper(6)
    c = response.center()
    q = operator.q
    accel = c["anchor_proper_Hubble_derivative_response"]
    ndd = c["anchor_lapse_response_second_derivative"]
    xi = c["anchor_hat_scale_second_derivative"]
    n = c["anchor_lapse_response"]
    coefficient_bound = lambda expr: sum(
        abs(co) * 4 ** p[0] for p, co in sp.Poly(expr, q).terms()
    )
    weighted_clock = 2 * r0 * (8 * prior_tails.half_line_integral_upper(3) + 30)
    mu = sp.Symbol("fixed_normalized_band_mean_squared_momentum", real=True)
    return {
        "inherited_actual_mean_generator_integrable_upper": inherited[
            "weighted_generator_integrable_entry_sum_upper"
        ],
        "source_specific_proxy_density_anchor_upper": r0,
        "source_specific_global_proxy_density_upper_per_eta": r0 / t**4,
        "actual_hat_lapse_forcing_absolute_upper_per_proxy_density": sp.Rational(3, 2),
        "global_weighted_mean_phase_upper_per_eta": C,
        "global_lapse_upper_per_eta": L,
        "global_linear_physical_log_scale_upper_per_eta": B,
        "global_actual_frame_log_scale_upper_per_eta": sp.Integer(8000),
        "global_matter_field_upper_per_eta": sp.Integer(1000),
        "global_scalar_density_upper_per_eta": sp.Integer(216),
        "global_relative_scalar_density_upper_per_eta": sp.Integer(43200),
        "global_lapse_tail_bound_per_eta": sp.Rational(3, 10) * C * sp.Abs(u) / t**3
        + sp.Rational(3, 40) * C / t**11
        + sp.Rational(25, 4) * r0 / t**3,
        "global_hat_scale_limit_error_per_eta_for_u_ge_one": sp.Rational(3, 220)
        * C
        / t**11
        + sp.Rational(241, 1000) * C / u**5
        + sp.Rational(25, 6) * r0 / t**3,
        "absolute_volume_weighted_clock_source_integral_upper_per_eta": weighted_clock,
        "rounded_absolute_volume_weighted_clock_source_integral_upper_per_eta": sp.Integer(
            280
        ),
        "future_volume_weighted_clock_source_integral_per_eta": 4 - mu,
        "whole_line_volume_weighted_clock_source_integral": sp.Integer(0),
        "band_mean_squared_momentum_domain": "1<mu<4 for the fixed normalized open-shell bump",
        "center_proper_acceleration_lower_per_eta": sp.Integer(-15),
        "center_proper_acceleration_upper_per_eta": sp.Integer(-9),
        "center_lapse_second_derivative_absolute_coefficient_bound": coefficient_bound(
            ndd
        ),
        "complete_first_order_frame_representative_not_a_SEE_solution": True,
        "checks": {
            "actual_clock_volume_coefficient_is_three_H_over_two": sp.factor(
                3 * prior_mean.data()["H"]
                - sp.diff(t**3, u) / t**3
                - sp.Rational(3, 2) * prior_mean.data()["H"]
            ),
            "clock_transfer_radial_time_integral": sp.integrate(u / t**2, (u, 0, sp.oo))
            - sp.Rational(1, 2),
            "global_mean_scalar_density_bound": sp.Rational(3, 100) * B - 216,
            "global_relative_matter_density_bound": 6 * B - 43200,
        },
        "bounds": {
            "inherited_weighted_generator_integral_still_below_two": inherited[
                "weighted_generator_integrable_entry_sum_upper"
            ]
            < 2,
            "new_tensor_global_weighted_phase_bound": 9 * 170 * r0 < C,
            "new_tensor_global_lapse_bound": lapse < L,
            "new_tensor_linear_physical_scale_bound": C + L / 2 < B,
            "new_tensor_actual_frame_scale_bound": C + L < 8000,
            "new_tensor_matter_field_bound": field < 1000,
            "global_clock_source_is_absolutely_integrable": weighted_clock < 280,
            "actual_center_acceleration_is_strictly_increasing_on_band": sp.diff(
                accel, q
            ).subs(q, 1)
            > 0
            and sp.diff(accel, q, 2) > 0,
            "actual_center_acceleration_lower": accel.subs(q, 1) > -15,
            "actual_center_acceleration_upper": accel.subs(q, 4) < -9,
            "actual_center_lapse_below_one_and_positive": n.subs(q, 1) < 1
            and n.subs(q, 4) > 0,
            "actual_center_hat_scale_second_derivative_positive_and_below_two": xi.subs(
                q, 1
            )
            < 2
            and xi.subs(q, 4) > 0,
            "actual_center_lapse_second_derivative_absolute_below_thirty": coefficient_bound(
                ndd
            )
            < 30,
            "calibrated_actual_frame_positive_everywhere": sp.Rational(8000, 10**20)
            < sp.Rational(1, 10**10),
            "calibrated_exact_frame_center_acceleration_numerator_above_three": 4
            - sp.Rational(35, 10**20)
            - sp.Rational(3, 2 * 10**40)
            > 3,
        },
    }
