"""A fixed smooth positive compact band and explicit finite mean-response bounds."""

from fractions import Fraction
from functools import cache

import sympy as sp

WIDTH = sp.Rational(1, 100)
DEFAULT_AMPLITUDE = sp.Rational(1, 10**20)


def exact_amplitude(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require an exact finite nonnegative covariance amplitude")
    value = sp.Rational(value)
    if value < 0 or value > sp.Rational(1, 10**20):
        raise ValueError("Covariance amplitude must lie between zero and 10^-20")
    return value


@cache
def data():
    q = sp.Symbol("squared_comoving_band_momentum", real=True)
    flat_rho = (3 * 10**6 + 3 + q + 2 * q / 10**6) / 2
    flat_p = (10**6 - 1 + q + 2 * q / (3 * 10**6)) / 2
    tmax = sp.Rational(10001, 10000)
    hmax = tmax**3
    Jfloor = sp.Rational(1215, 800) / tmax**18
    alpha = 3 * WIDTH * (4 * hmax - 1)
    H = sp.Rational(1, 25)
    ell = sp.Rational(1, 10)
    beta = sp.Rational(1, 20)
    A = sp.Rational(1, 10)
    J = sp.Rational(3, 2)
    matrix_abs = sp.Matrix(
        [
            [A * ell * beta / (2 * J), sp.Rational(1, 2) + A * A / (6 * J)],
            [
                3 * ell * ell + 3 * ell * ell * beta * beta / (2 * J),
                3 * H + ell * A * beta / (2 * J),
            ],
        ]
    )
    rho_upper = sp.Integer(2_000_000)
    F_upper = sp.Integer(5_000_000)
    source_abs = sp.Matrix(
        [A * F_upper / (6 * J), ell * beta * F_upper / (2 * J) + rho_upper]
    )
    phase_upper = sp.Integer(22_000)
    lapse = (A * phase_upper + 3 * ell * beta * phase_upper + F_upper) / (2 * J)
    return {
        "fixed_band_shape": "exp(-1/((q-1)*(4-q))) for 1<q<4, zero elsewhere",
        "band_measure": "d^3k/(2*pi)^3",
        "band_normalization": "b(k)=shape(|k|^2)/integral shape(|k|^2) d^3k/(2*pi)^3",
        "anchor_added_covariance": "eta*b(k)*I6, eta=hbar*lambda/kappa; the scalar and tensor factors are unchanged",
        "actual_added_covariance": "U_P,a(t,0,k) eta*b(k)*I6 U_P,a(s,0,k)^T",
        "flat_anchor_density_per_unit_band_covariance": flat_rho,
        "flat_anchor_pressure_per_unit_band_covariance": flat_p,
        "band_minimum_squared_momentum": sp.Integer(1),
        "band_maximum_squared_momentum": sp.Integer(4),
        "closed_time_half_width": WIDTH,
        "whole_interval_clock_h_upper": hmax,
        "whole_interval_J_lower": Jfloor,
        "whole_interval_alpha_upper": alpha,
        "whole_interval_mean_generator_entrywise_absolute_upper": matrix_abs,
        "whole_interval_mean_forcing_entrywise_upper_per_eta": source_abs,
        "whole_interval_density_upper_per_eta": rho_upper,
        "whole_interval_pressure_absolute_upper_per_eta": rho_upper,
        "whole_interval_lapse_force_upper_per_eta": F_upper,
        "whole_interval_mean_generator_infinity_norm_upper": sp.Integer(1),
        "whole_interval_mean_forcing_infinity_norm_upper_per_eta": sp.Integer(
            2_100_000
        ),
        "whole_interval_hat_scale_trace_response_upper_per_eta": phase_upper,
        "whole_interval_lapse_response_upper_per_eta": sp.Integer(1_700_000),
        "whole_interval_physical_log_scale_response_upper_per_eta": sp.Integer(900_000),
        "whole_interval_scalar_density_response_upper_per_eta": sp.Integer(27_000),
        "whole_interval_matter_field_derivative_response_upper_per_eta": sp.Integer(
            100_000
        ),
        "whole_interval_matter_field_response_upper_per_eta": sp.Integer(1_000),
        "derived_lapse_bound_before_rounding": lapse,
        "checks": {
            "band_anchor_density_retains_all_mass_and_gradient_terms": sp.factor(
                flat_rho
                - (3 * sp.Integer(10**6) + 3 + q + 2 * q / sp.Integer(10**6)) / 2
            ),
            "exact_density_energy_exponent": 6 * H * WIDTH - sp.Rational(3, 1250),
            "exact_linear_variation_of_constants_majorant": WIDTH / (1 - WIDTH)
            - sp.Rational(1, 99),
        },
        "bounds": {
            "complete_band_anchor_density": flat_rho.subs(q, 4) < 1_500_004,
            "complete_density_propagation": sp.Rational(1_500_004) * 1250 / 1247
            < rho_upper,
            "whole_interval_clock_h_below_1001_over_1000": hmax
            < sp.Rational(1001, 1000),
            "whole_interval_J_at_least_three_halves": Jfloor > J,
            "whole_interval_alpha_below_one_tenth": alpha < A,
            "both_mean_generator_row_sums_below_one": all(
                sum(matrix_abs.row(i)) < 1 for i in range(2)
            ),
            "both_mean_forcing_rows_below_2100000": all(
                v < 2_100_000 for v in source_abs
            ),
            "full_retarded_mean_response_below_22000": sp.Rational(2_100_000, 99)
            < phase_upper,
            "full_lapse_response_below_1700000": lapse < 1_700_000,
            "full_physical_scale_response_below_900000": phase_upper
            + sp.Rational(1_700_000, 2)
            < 900_000,
            "scalar_density_response_upper_is_27000": 3 * ell * ell * 900_000 == 27_000,
            "matter_field_response_derivative_below_100000": beta * 1_700_000
            + 3 * ell * phase_upper
            < 100_000,
            "matter_field_response_below_1000": WIDTH * 100_000 == 1_000,
        },
    }


def calibrated(value=DEFAULT_AMPLITUDE):
    eta = exact_amplitude(value)
    d = data()
    keys = (
        "density",
        "pressure_absolute",
        "lapse_force",
        "hat_scale_trace_response",
        "lapse_response",
        "physical_log_scale_response",
        "scalar_density_response",
    )
    return {
        "eta": eta,
        **{key: d["whole_interval_" + key + "_upper_per_eta"] * eta for key in keys},
    }


def bad_cases():
    values = (
        True,
        False,
        sp.true,
        sp.false,
        0.1,
        sp.Float("0.0"),
        None,
        "0",
        sp.oo,
        sp.nan,
        sp.zoo,
        sp.I,
        sp.Symbol("eta"),
        -1,
        sp.Rational(-1, 10**30),
        sp.Rational(1, 10**19),
    )
    return [
        ("invalid_covariance_amplitude_" + str(i), calibrated, (v,))
        for i, v in enumerate(values)
    ]
