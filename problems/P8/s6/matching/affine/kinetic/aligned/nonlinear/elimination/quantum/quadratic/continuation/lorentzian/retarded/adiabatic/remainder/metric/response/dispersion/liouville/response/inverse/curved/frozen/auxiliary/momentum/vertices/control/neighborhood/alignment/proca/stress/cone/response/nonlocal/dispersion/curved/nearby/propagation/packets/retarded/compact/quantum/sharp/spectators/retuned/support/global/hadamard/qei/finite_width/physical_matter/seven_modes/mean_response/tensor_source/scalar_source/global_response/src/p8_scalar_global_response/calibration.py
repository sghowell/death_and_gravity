"""Exact amplitudes for the same scalar state family, on a smaller global domain."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import band as parent

from . import bounds, phase

DEFAULT_AMPLITUDE = sp.Rational(1, 10**80)


def exact_amplitude(value):
    eta = parent.exact_amplitude(value)
    if eta > DEFAULT_AMPLITUDE:
        raise ValueError("Global scalar amplitude must be at most 10^-80")
    return eta


def calibrated(value=DEFAULT_AMPLITUDE):
    eta = exact_amplitude(value)
    d = bounds.data()
    keys = (
        "weighted_mean_phase",
        "global_lapse",
        "global_linear_physical_scale",
        "global_actual_frame_scale",
        "global_matter_field",
    )
    return {
        "eta": eta,
        **{key: d[key + "_upper_per_eta"] * eta for key in keys},
        "physical_density_and_pressure_scaled_by_t_four": 10**64 * eta,
    }


@cache
def bad_cases():
    amplitude = [(name, exact_amplitude, args) for name, _, args in parent.bad_cases()]
    amplitude += [
        (
            "scalar_amplitude_in_local_but_not_global_domain",
            exact_amplitude,
            (sp.Rational(1, 10**79),),
        )
    ]
    u, k = phase.u, phase.k
    bad = (
        True,
        False,
        0.1,
        sp.Float("0"),
        sp.oo,
        sp.nan,
        sp.zoo,
        sp.I,
        "0",
        None,
        sp.Symbol("undeclared"),
        sp.exp(u),
        1 / (1 - u * u),
        1 / u**2,
        1 / (1 + u),
        1 / (1 + k * k),
        sp.sqrt(2),
        (1 + u) / (1 + u * u) ** 2,
        1,
        u / (1 + u * u),
    )
    return (
        amplitude
        + [
            ("invalid_global_envelope_" + str(i), phase.envelope, (value,))
            for i, value in enumerate(bad)
        ]
        + [
            ("invalid_integral_exponent_" + str(i), phase.half_line_integral, (value,))
            for i, value in enumerate(
                (
                    True,
                    0.5,
                    sp.Float("1.5"),
                    sp.oo,
                    sp.nan,
                    "2",
                    0,
                    sp.Rational(1, 2),
                    sp.Rational(5, 3),
                )
            )
        ]
        + [
            (
                "invalid_momentum_half_order_" + str(i),
                phase.envelope,
                (1 / (1 + u * u), value),
            )
            for i, value in enumerate((True, -1, 3, 0.0, sp.Rational(1, 2)))
        ]
    )
