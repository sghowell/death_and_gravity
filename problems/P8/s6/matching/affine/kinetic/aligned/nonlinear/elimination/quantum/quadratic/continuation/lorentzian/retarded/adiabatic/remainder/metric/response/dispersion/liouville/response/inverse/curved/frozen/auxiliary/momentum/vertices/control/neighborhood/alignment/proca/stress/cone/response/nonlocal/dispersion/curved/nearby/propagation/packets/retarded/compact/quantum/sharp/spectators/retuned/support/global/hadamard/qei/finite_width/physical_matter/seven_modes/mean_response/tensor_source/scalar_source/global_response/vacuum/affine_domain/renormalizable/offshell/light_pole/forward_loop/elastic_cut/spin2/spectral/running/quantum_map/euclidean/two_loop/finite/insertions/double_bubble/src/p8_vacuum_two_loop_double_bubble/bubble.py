"""The inherited zero-momentum light bubble and a full complex-disc bound."""

from functools import cache

import sympy as sp


@cache
def data():
    s = sp.Symbol("complex_channel_invariant")
    x, t = sp.symbols("unit_Feynman_parameter unit_auxiliary_parameter", real=True)
    a = x * (1 - x)
    logarithm = -sp.log(1 - a * s)
    primitive = -sp.log(1 - t * a * s)
    return {
        "s": s,
        "x": x,
        "finite_light_bubble": sp.Integral(logarithm, (x, 0, 1)) / (16 * sp.pi**2),
        "finite_logarithmic_parameter_integrand": logarithm,
        "complete_complex_logarithm_upper": sp.Integer(2),
        "finite_light_bubble_modulus_upper": 2 / (16 * sp.pi**2),
        "positive_exp_two_partial_sum": 1 + sp.Integer(2) + sp.Rational(2**2, 2),
        "scope": "For |s|<=3 the anchored principal logarithm is on the first sheet and its modulus is at most ln(4)<2. This is I(s)-I(0) with the entire same regulated reference, not a value assigned to I(0).",
        "checks": {
            "anchored_logarithm_derivative": sp.factor(
                sp.diff(primitive, t) - a * s / (1 - t * a * s)
            ),
            "anchored_logarithm_zero": primitive.subs(t, 0),
            "anchored_logarithm_endpoint": primitive.subs(t, 1) - logarithm,
            "unit_parameter_weight_maximum": sp.expand(
                sp.Rational(1, 4) - a - (x - sp.Rational(1, 2)) ** 2
            ),
            "disc_denominator_gap": 1 - sp.Rational(3, 4) - sp.Rational(1, 4),
            "zero_transfer_finite_bubble_vanishes": logarithm.subs(s, 0),
            "bubble_measure_uses_no_quartic_symmetry_factor": 2 / (16 * sp.pi**2)
            - 1 / (8 * sp.pi**2),
            "positive_exp_two_partial_sum": 1 + 2 + sp.Rational(2**2, 2) - 5,
        },
    }
