"""Exact-radius bounds for the same on-shell-subtracted one-loop kernel."""

from fractions import Fraction

import sympy as sp

from . import kernel


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational")
    return sp.Rational(value)


def point(radius=1):
    r = rational(radius)
    d = kernel.data()
    M = d["actual_heavy_mass_squared"]
    g = d["actual_cubic_squared"]
    if not 0 < r < M:
        raise ValueError(
            "Require a positive radius strictly below the proven parameter radius"
        )
    quadratic = g / (864 * M * M * (1 - r / M))
    return {
        "complex_disc_radius_about_mass_one": r,
        "quadratic_self_energy_upper_coefficient": quadratic,
        "factored_inverse_propagator_error_upper": r * quadratic,
        "no_extra_one_loop_pole_proved_in_disc": bool(r * quadratic < 1),
        "finite_kinetic_counterterm_upper": g / (288 * M),
        "positive_curvature_conversion_upper": g / (864 * M * M),
        "scope": "A radius below M proves parameter analyticity, not automatically absence of additional zeros. That conclusion additionally requires the displayed factored error to be below one. No full higher-loop cut or pole statement.",
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        sp.I,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.nan,
        None,
        sp.Symbol("r"),
    )
    rows = [
        ("inexact_or_nonreal_radius_" + str(i), point, (v,))
        for i, v in enumerate(invalid)
    ]
    M = kernel.data()["actual_heavy_mass_squared"]
    rows += [
        ("outside_parameter_disc_" + str(i), point, (v,))
        for i, v in enumerate((0, -1, M, M + 1))
    ]
    return rows
