"""Exact rational whole-strip bounds for the actual scalar generator and sources."""

from functools import cache

import sympy as sp

from . import model, observable, response, source

WIDTH = sp.Rational(1, 100)


def exact_rational_bound(value):
    if isinstance(value, bool) or not isinstance(value, (int, sp.Expr)):
        raise TypeError(
            "Require an exact rational expression in actual time and momentum"
        )
    value = sp.sympify(value)
    if value.has(
        sp.Float, sp.oo, -sp.oo, sp.nan, sp.zoo, sp.I
    ) or not value.free_symbols <= {model.u, model.k}:
        raise ValueError(
            "Require real exact coefficients and only the declared variables"
        )
    numerator, denominator = sp.fraction(sp.cancel(value))
    try:
        den = sp.Poly(denominator, model.u, domain=sp.QQ)
        num = sp.Poly(numerator, model.u, model.k, domain=sp.QQ)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed) as exc:
        raise ValueError(
            "Require a rational numerator and a time-only polynomial denominator"
        ) from exc
    if any(co < 0 or power[0] % 2 for power, co in den.terms()):
        raise ValueError("Require a positive even denominator coefficient proof")
    if den.nth(0) <= 0:
        raise ValueError("Require a positive denominator constant")
    upper = sum(
        abs(co) * WIDTH ** power[0] * 2 ** power[1] for power, co in num.terms()
    )
    return upper / den.nth(0)


@cache
def data():
    d = source.data()
    matrix = d["natural_phase_generator"].applyfunc(exact_rational_bound)
    rows = [sum(matrix.row(j)) for j in range(4)]
    kernels = {
        key: value.applyfunc(exact_rational_bound)
        for key, value in d["source_Hessians"].items()
    }
    sums = {key: sum(value) / 2 for key, value in kernels.items()}
    return {
        "generator_entry_bounds": matrix,
        "generator_row_sums": rows,
        "source_half_entry_sums": sums,
        "generator_max_row": max(rows),
    }


@cache
def response_bounds():
    d = data()
    eta = sp.Rational(1, 10**20)
    intrinsic = {
        key: sum(observable.data()[key].applyfunc(exact_rational_bound)) / 2
        for key in ("intrinsic_density_Hessian", "intrinsic_pressure_Hessian")
    }
    c = response.center()
    polynomial_bound = lambda val: sum(
        abs(co) * 2 ** powers[0] for powers, co in sp.Poly(val, model.k).terms()
    )
    center_bounds = {
        key: polynomial_bound(c[key])
        for key in (
            "anchor_induced_lapse",
            "anchor_hat_scale_second_derivative",
            "anchor_lapse_second_derivative",
            "anchor_proper_Hubble_derivative_response",
        )
    }
    forcing_xi = sp.Rational(350, 90) + 10
    forcing_p = sp.Rational(350, 600) + 400
    lapse = (sp.Rational(1, 10) * 5 + sp.Rational(3, 200) * 5 + 350) / 3
    field = sp.Rational(1, 20) * 120 + sp.Rational(3, 10) * 5 + 60
    return {
        "closed_time_half_width": WIDTH,
        "closed_momentum_root_interval": [sp.Integer(1), sp.Integer(2)],
        "exact_generator_infinity_norm_upper": sp.Integer(90),
        "actual_fundamental_matrix_infinity_norm_upper": sp.Integer(3),
        "added_covariance_each_entry_absolute_upper_per_eta": sp.Integer(9),
        "four_scalar_source_absolute_uppers_per_eta": {
            "lapse": sp.Integer(350),
            "scale": sp.Integer(10),
            "trace": sp.Integer(400),
            "matter_field": sp.Integer(60),
        },
        "complete_two_component_mean_forcing_upper_per_eta": sp.Integer(410),
        "mean_scale_and_trace_upper_per_eta": sp.Integer(5),
        "mean_lapse_upper_per_eta": sp.Integer(120),
        "mean_linear_physical_log_scale_upper_per_eta": sp.Integer(66),
        "mean_exact_frame_log_scale_upper_per_eta": sp.Integer(130),
        "mean_matter_field_derivative_upper_per_eta": sp.Integer(70),
        "mean_matter_field_upper_per_eta": sp.Integer(1),
        "intrinsic_observable_half_entry_sum_bounds": intrinsic,
        "intrinsic_physical_density_and_pressure_absolute_upper_per_eta": sp.Integer(
            36
        ),
        "additional_mean_density_and_pressure_absolute_upper_per_eta": sp.Integer(2),
        "complete_physical_density_and_pressure_absolute_upper_per_eta": sp.Integer(40),
        "exact_center_coefficient_bounds": center_bounds,
        "center_lapse_absolute_upper_per_eta": sp.Integer(10),
        "center_hat_scale_second_derivative_absolute_upper_per_eta": sp.Integer(40),
        "center_lapse_second_derivative_absolute_upper_per_eta": sp.Integer(200),
        "center_proper_Hubble_derivative_absolute_upper_per_eta": sp.Integer(200),
        "checks": {
            "whole_strip_induced_matter_field_integration": 70 * WIDTH
            - sp.Rational(7, 10),
            "whole_strip_additional_mean_physical_density": sp.Rational(3, 100) * 66
            - sp.Rational(99, 50),
        },
        "bounds": {
            "actual_generator_whole_strip_bound": d["generator_max_row"] < 90,
            "fundamental_matrix_exponent_below_one": 90 * WIDTH < 1,
            "lapse_covariance_source_bound": 9
            * d["source_half_entry_sums"]["actual_mean_lapse_source"]
            < 350,
            "scale_covariance_source_bound": 9
            * d["source_half_entry_sums"]["mean_hat_scale_direct_source"]
            < 10,
            "trace_covariance_source_bound": 9
            * d["source_half_entry_sums"]["mean_trace_direct_source"]
            < 400,
            "matter_field_covariance_source_bound": 9
            * d["source_half_entry_sums"]["mean_matter_field_direct_source"]
            < 60,
            "mean_forcing_bound": max(forcing_xi, forcing_p) < 410,
            "mean_retarded_phase_bound": sp.Rational(410, 99) < 5,
            "mean_constraint_lapse_bound": lapse < 120,
            "mean_physical_linear_scale_bound": 5 + sp.Rational(120, 2) < 66,
            "mean_exact_frame_scale_bound": 5 + 120 < 130,
            "mean_field_derivative_bound": field < 70,
            "mean_field_bound": 70 * WIDTH < 1,
            "intrinsic_physical_density_bound": 9
            * intrinsic["intrinsic_density_Hessian"]
            < 36,
            "intrinsic_physical_pressure_bound": 9
            * intrinsic["intrinsic_pressure_Hessian"]
            < 36,
            "additional_mean_density_bound": sp.Rational(99, 50) < 2,
            "complete_physical_density_pressure_bound": 36 + 2 < 40,
            "center_lapse_absolute_bound": center_bounds["anchor_induced_lapse"] < 10,
            "center_hat_scale_second_derivative_bound": center_bounds[
                "anchor_hat_scale_second_derivative"
            ]
            < 40,
            "center_lapse_second_derivative_bound": center_bounds[
                "anchor_lapse_second_derivative"
            ]
            < 200,
            "center_acceleration_absolute_bound": center_bounds[
                "anchor_proper_Hubble_derivative_response"
            ]
            < 200,
            "calibrated_frame_remains_positive_on_strip": 130 * eta
            < sp.Rational(1, 10**10),
            "exact_frame_center_bounce_numerator_above_three": 4
            - 300 * eta
            - 150 * eta * eta
            > 3,
        },
    }
