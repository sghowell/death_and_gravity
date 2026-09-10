"""Exact physical invariant inputs and positive cut-density enclosures."""

from fractions import Fraction

import sympy as sp

from . import amplitude


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational physical invariant")
    return sp.Rational(value)


def point(invariant=5):
    s = rational(invariant)
    if not 4 <= s <= 6:
        raise ValueError(
            "Require the physical invariant in the certified interval [4,6]"
        )
    lam = amplitude.data()["actual_fixed_lambda"]
    beta = sp.sqrt(1 - 4 / s)
    return {
        "physical_invariant": s,
        "two_body_phase_beta": beta,
        "cut_density_lower": beta * (23 * lam) ** 2 / 128,
        "cut_density_upper": beta * (73 * lam) ** 2 / 96,
        "cut_density_strictly_positive": bool(s > 4),
        "threshold_density_exactly_zero": bool(s == 4),
        "scope": "Exact invariant-dependent enclosure of the same one-loop elastic density only, with the continuous threshold value handled directly",
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "5",
        sp.I,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.nan,
        None,
        sp.Symbol("s"),
    )
    rows = [
        ("inexact_or_nonreal_invariant_" + str(i), point, (v,))
        for i, v in enumerate(invalid)
    ]
    rows += [
        ("outside_physical_cut_window_" + str(i), point, (v,))
        for i, v in enumerate((0, 3, 7, -1))
    ]
    return rows
