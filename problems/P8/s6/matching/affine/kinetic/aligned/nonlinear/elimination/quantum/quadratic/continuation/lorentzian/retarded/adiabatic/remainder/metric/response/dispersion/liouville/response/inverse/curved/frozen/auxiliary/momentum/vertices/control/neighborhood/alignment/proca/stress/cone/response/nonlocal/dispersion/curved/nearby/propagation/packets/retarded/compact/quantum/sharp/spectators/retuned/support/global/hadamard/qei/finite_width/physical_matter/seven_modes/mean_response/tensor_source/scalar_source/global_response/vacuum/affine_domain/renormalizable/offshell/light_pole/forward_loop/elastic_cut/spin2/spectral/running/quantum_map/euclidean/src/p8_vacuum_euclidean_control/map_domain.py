"""Actual derivative-coordinate multiplier is small locally, not uniformly."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_quantum_map import calibration as previous


@cache
def data():
    p = model.data()["actual_parameters"]
    lam, gamma, c = sp.symbols(
        "positive_lambda positive_gamma redundant_c", positive=True
    )
    y = sp.Symbol("Euclidean_momentum_squared", nonnegative=True)
    N = c / 2 - 3 * lam + gamma / 2 - (lam + 5 * gamma / 4) * y + gamma * y * y / 2
    minimum = 3 * lam * lam / (2 * gamma) - 17 * lam / 4 - 9 * gamma / 32
    center = lam / gamma + sp.Rational(5, 4)
    square = gamma * (y - center) ** 2 / 2 + minimum
    actual = {
        lam: p["lambda"],
        gamma: p["gamma"],
        c: previous.data()["actual_redundant_c"],
    }
    actualN = sp.expand(N.subs(actual))
    M = p["heavy_mass_squared"]
    window = (
        c / 2 + 3 * lam + gamma / 2 + (lam + 5 * gamma / 4) * M + gamma * M * M / 2
    ).subs(actual) / 144
    witness = sp.Integer(10) ** 400
    witnessN = sp.factor(actualN.subs(y, witness))
    old = previous.data()["actual_finite_composite_polynomial"].subs(
        sp.Symbol("free_Box_eigenvalue", real=True), y
    )
    return {
        "y": y,
        "general_composite_numerator": N,
        "completed_square": square,
        "minimum_numerator": minimum,
        "actual_composite_numerator": actualN,
        "actual_minimum_numerator": minimum.subs(actual),
        "actual_heavy_window_endpoint": M,
        "actual_heavy_window_mixing_upper": window,
        "explicit_nonuniformity_witness_invariant": witness,
        "witness_numerator": witnessN,
        "witness_mixing_strict_rational_lower": witnessN / 256,
        "selected_Planck_squared": sp.Integer(10) ** 800,
        "scope": "The explicit composite coordinate fails uniform perturbative smallness at sufficiently large Euclidean momentum. This is not a physical ghost, breakdown of on-shell equivalence, a quantum cutoff determination or a no-go for a P8 row.",
        "checks": {
            "actual_parent_composite_polynomial": sp.expand(
                old - actualN / (16 * sp.pi**2)
            ),
            "general_completed_square_with_actual_c_relation": sp.expand(
                N.subs(c, 4 * lam**2 / gamma) - square
            ),
            "minimum_at_square_center": sp.factor(
                N.subs(c, 4 * lam**2 / gamma).subs(y, center) - minimum
            ),
            "positive_large_momentum_leading_coefficient": sp.limit(N / y**2, y, sp.oo)
            - gamma / 2,
            "same_actual_massive_reference_scale": M - p["heavy_mass_squared"],
        },
        "bounds": {
            "actual_completed_square_minimum_positive": minimum.subs(actual) > 0,
            "heavy_window_field_map_mixing_below_one_e_minus_404": 0
            < window
            < sp.Rational(1, 10**404),
            "explicit_witness_above_heavy_window": witness > M,
            "explicit_witness_below_selected_Planck_squared": witness
            < sp.Integer(10) ** 800,
            "witness_numerator_strictly_above_256": witnessN > 256,
            "witness_mixing_strictly_above_one": witnessN / 256 > 1,
        },
    }
