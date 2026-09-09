"""Exact compact-strip calibration and rejection controls."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import band as parent

from . import bounds, model


def calibrated(value=parent.DEFAULT_AMPLITUDE):
    eta = parent.exact_amplitude(value)
    d = bounds.response_bounds()
    return {
        "eta": eta,
        **{
            key.removesuffix("_per_eta"): d[key] * eta
            for key in (
                "mean_scale_and_trace_upper_per_eta",
                "mean_lapse_upper_per_eta",
                "mean_linear_physical_log_scale_upper_per_eta",
                "mean_exact_frame_log_scale_upper_per_eta",
                "mean_matter_field_upper_per_eta",
                "complete_physical_density_and_pressure_absolute_upper_per_eta",
                "center_proper_Hubble_derivative_absolute_upper_per_eta",
            )
        },
    }


@cache
def bad_cases():
    amplitude = [(name, calibrated, args) for name, _, args in parent.bad_cases()]
    u, k = model.u, model.k
    invalid = (
        True,
        False,
        0.1,
        sp.Float("0.0"),
        sp.oo,
        sp.nan,
        sp.I,
        "0",
        None,
        sp.Symbol("undeclared_coefficient"),
        sp.exp(u),
        1 / (1 - u * u),
        1 / (u * u),
        1 / (1 + u),
        1 / (1 + k * k),
        sp.sqrt(2),
    )
    return amplitude + [
        (
            "invalid_rational_strip_bound_" + str(i),
            bounds.exact_rational_bound,
            (value,),
        )
        for i, value in enumerate(invalid)
    ]
