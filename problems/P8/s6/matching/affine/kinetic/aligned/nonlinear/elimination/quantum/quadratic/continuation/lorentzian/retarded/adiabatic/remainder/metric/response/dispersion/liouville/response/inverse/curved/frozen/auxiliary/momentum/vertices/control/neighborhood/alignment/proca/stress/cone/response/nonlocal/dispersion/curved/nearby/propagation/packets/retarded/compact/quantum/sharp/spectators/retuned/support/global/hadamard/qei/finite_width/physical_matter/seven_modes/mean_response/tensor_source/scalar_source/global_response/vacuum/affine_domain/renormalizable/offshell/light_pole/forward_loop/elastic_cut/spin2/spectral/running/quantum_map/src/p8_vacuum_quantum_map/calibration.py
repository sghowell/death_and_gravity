"""Actual composite overlap and local one-loop propagator-disc bounds."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_offshell_vacuum import matching
from p8_polynomial_vacuum import model
from p8_vacuum_light_pole import kernel

from . import wick


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require an exact finite rational mass-squared invariant")
    return sp.Rational(value)


def point(invariant=1):
    s = rational(invariant)
    if not 0 <= s <= 2:
        raise ValueError("The declared real slice is invariant zero through two")
    d = data()
    B = sp.Symbol("free_Box_eigenvalue", real=True)
    return {
        "mass_squared_invariant": s,
        "finite_linear_composite_mixing": d["actual_finite_composite_polynomial"].subs(
            B, -s
        ),
        "whole_complex_disc_mixing_upper": d["unit_disc_composite_mixing_upper"],
        "scope": "This evaluates the one-loop composite polynomial, not an all-order propagator or a new resummed higher-derivative theory.",
    }


@cache
def data():
    p = model.data()["actual_parameters"]
    lam, gamma = p["lambda"], p["gamma"]
    m = matching.data()
    c = m["redundant_c"].subs({m["lambda"]: lam, m["gamma"]: gamma})
    jl, jg, jc = sp.symbols("quartic_lambda quartic_gamma redundant_c", real=True)
    B = sp.Symbol("free_Box_eigenvalue", real=True)
    W = sp.expand(
        wick.data()["MSbar_finite_linear_composite_mixing"].subs(
            {jl: lam, jg: gamma, jc: c}
        )
    )
    w = sp.factor(W.subs(B, -1))
    numerator = c / 2 - 2 * lam + 9 * gamma / 4
    upper = (c / 2 + 5 * lam + 5 * gamma) / 144
    on_shell_upper = numerator / 144
    old_epsilon = kernel.data()["actual_unit_disc_quadratic_self_energy_upper"]
    epsilon = 2 * upper + old_epsilon
    s, hbar = sp.symbols("Minkowski_invariant loop_order")
    old_sigma = sp.Symbol("parent_on_shell_self_energy")
    wp = sp.Function("finite_composite_polynomial")
    G_new = 1 / (s - 1) - hbar * old_sigma / (s - 1) ** 2 - 2 * hbar * wp(s) / (s - 1)
    D_new = s - 1 + hbar * (old_sigma + 2 * (s - 1) * wp(s))
    map_sigma = 2 * (s - 1) * wp(s)
    return {
        "actual_lambda": lam,
        "actual_gamma": gamma,
        "actual_redundant_c": c,
        "actual_finite_composite_polynomial": W,
        "actual_on_shell_composite_overlap": w,
        "actual_overlap_positive_numerator": numerator,
        "actual_on_shell_overlap_rational_upper": on_shell_upper,
        "unit_disc_composite_mixing_upper": upper,
        "parent_unit_disc_self_energy_factor_upper": old_epsilon,
        "mapped_unit_disc_inverse_factor_error_upper": epsilon,
        "mapped_one_loop_residue": 1 - 2 * w,
        "mapped_residue_rational_lower": 1 - 2 * on_shell_upper,
        "formal_mapped_connected_propagator": G_new,
        "formal_mapped_inverse_propagator": D_new,
        "scope": "Mass-one light pole, explicit MSbar composite prescription and |s-1|<=1 through one loop only. Positive residue and no extra zeros refer to the displayed truncated inverse in this disc. No global differential inverse, resummed higher-derivative spectrum, higher-loop pole theorem or rolling-state transfer.",
        "checks": {
            "same_actual_S711_coefficient": sp.factor(c - 4 * lam**2 / gamma),
            "actual_on_shell_finite_mixing": sp.factor(w - numerator / (16 * sp.pi**2)),
            "actual_disc_polynomial_majorant": sp.expand(
                144 * upper - c / 2 - 5 * lam - 5 * gamma
            ),
            "full_one_loop_inverse_multiplication": sp.cancel(
                sp.expand(G_new * D_new).coeff(hbar, 1)
            ),
            "field_map_self_energy_vanishes_at_light_pole": map_sigma.subs(s, 1),
            "field_map_self_energy_derivative_retains_overlap": sp.diff(
                map_sigma, s
            ).subs(s, 1)
            - 2 * wp(1),
            "inverse_residue_matches_composite_source_factor": sp.series(
                1 / (1 + 2 * hbar * wp(1)), hbar, 0, 2
            ).removeO()
            - 1
            + 2 * hbar * wp(1),
            "mapped_disc_factor_error_keeps_both_sources": epsilon
            - 2 * upper
            - old_epsilon,
        },
        "bounds": {
            "actual_composite_overlap_strictly_positive": numerator > 0,
            "actual_overlap_upper_below_one_e_minus_404": 0
            < on_shell_upper
            < sp.Rational(1, 10**404),
            "whole_complex_disc_mixing_below_one_e_minus_404": 0
            < upper
            < sp.Rational(1, 10**404),
            "actual_residue_rational_lower_positive": 1 - 2 * on_shell_upper > 0,
            "whole_disc_inverse_factor_error_below_three_e_minus_404": 0
            < epsilon
            < sp.Rational(3, 10**404),
            "whole_disc_inverse_factor_cannot_vanish": epsilon < 1,
            "parent_on_shell_disc_bound_remains_below_one_e_minus_405": old_epsilon
            < sp.Rational(1, 10**405),
        },
    }


def bad_cases():
    bad = (
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
        sp.Symbol("s"),
    )
    rows = [("inexact_invariant_" + str(i), point, (v,)) for i, v in enumerate(bad)]
    rows += [
        ("outside_invariant_" + str(i), point, (v,))
        for i, v in enumerate((-1, 3, Fraction(-1, 2), Fraction(5, 2)))
    ]
    return rows
