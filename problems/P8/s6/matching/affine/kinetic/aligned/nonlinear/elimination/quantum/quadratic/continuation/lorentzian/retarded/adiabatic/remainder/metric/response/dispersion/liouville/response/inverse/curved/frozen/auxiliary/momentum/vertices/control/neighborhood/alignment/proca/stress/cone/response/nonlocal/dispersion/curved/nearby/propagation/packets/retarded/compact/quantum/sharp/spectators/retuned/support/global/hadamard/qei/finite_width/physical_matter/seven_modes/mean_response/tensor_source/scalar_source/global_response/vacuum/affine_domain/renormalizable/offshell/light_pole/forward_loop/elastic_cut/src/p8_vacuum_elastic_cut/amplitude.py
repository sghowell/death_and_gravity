"""Exact physical tree amplitude and strict continuous-window envelopes."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_uv import vacuum


@cache
def data():
    s = sp.Symbol("physical_s", positive=True)
    z = sp.Symbol("physical_cosine", real=True)
    M, g, lam = sp.symbols(
        "heavy_mass_squared cubic_squared polynomial_quartic", positive=True
    )
    r = sp.Symbol("cosine_squared", nonnegative=True)
    k = (s - 4) / 2
    B = M + k
    t = -k * (1 - z)
    u = -k * (1 + z)
    C = -lam + g / (M - s)
    A = -lam + g * (1 / (M - s) + 1 / (M - t) + 1 / (M - u))
    compact = C + 2 * g * B / (B * B - k * k * z * z)
    even = C + 2 * g * B / (B * B - k * k * r)
    derivative = (
        g / (M - s) ** 2
        - g * ((1 - z) / (M + k * (1 - z)) ** 2 + (1 + z) / (M + k * (1 + z)) ** 2) / 2
    )
    parent = model.data()
    actual_map = {
        M: parent["heavy_mass_squared"],
        g: parent["cubic_coupling_squared"],
        lam: parent["polynomial_quartic"],
    }
    actual_tree = parent["exact_tree_amplitude"].subs(
        {vacuum.s: s, vacuum.transfer: t, vacuum.w: u}, simultaneous=True
    )
    p = parent["actual_parameters"]
    fixed = {
        M: p["heavy_mass_squared"],
        g: p["cubic_coupling_squared"],
        lam: p["bare_polynomial_quartic"],
    }
    threshold = sp.factor(A.subs({s: 4, z: 0}).subs(fixed))
    positive_window = sp.factor(A.subs({s: 5, z: 0}).subs(fixed))
    maximum = sp.factor(A.subs({s: 6, z: 1}).subs(fixed))
    ell = p["lambda"]
    return {
        "s": s,
        "cosine": z,
        "M": M,
        "g": g,
        "lambda4": lam,
        "t": t,
        "u": u,
        "physical_tree_amplitude": A,
        "compact_physical_tree_amplitude": compact,
        "selected_invariant_window": [sp.Integer(4), sp.Integer(6)],
        "strict_positive_subwindow": [sp.Integer(5), sp.Integer(6)],
        "actual_threshold_minimum": threshold,
        "actual_positive_subwindow_minimum": positive_window,
        "actual_full_window_maximum": maximum,
        "actual_fixed_lambda": ell,
        "actual_heavy_mass_squared": p["heavy_mass_squared"],
        "scope": "Exact physical on-shell tree vertices retained on the cut; no derivative truncation or squared unstable heavy resonance is integrated.",
        "checks": {
            "physical_mass_one_Mandelstam_sum": sp.expand(s + t + u - 4),
            "actual_parent_tree_vertex_rebuilt": sp.factor(
                actual_tree - A.subs(actual_map)
            ),
            "exact_even_angle_form": sp.factor(A - compact),
            "cosine_squared_form_not_linear_cosine_monotonicity": sp.factor(
                compact - even.subs(r, z * z)
            ),
            "physical_s_derivative": sp.factor(sp.diff(A, s) - derivative),
            "positive_cosine_squared_derivative": sp.factor(
                sp.diff(even, r) - 2 * g * B * k * k / (B * B - k * k * r) ** 2
            ),
            "positive_s_derivative_conservative_gap": sp.factor(
                1 / (M - s) ** 2 - 1 / M**2 - s * (2 * M - s) / (M * M * (M - s) ** 2)
            ),
        },
        "bounds": {
            "window_strictly_below_heavy_resonance": p["heavy_mass_squared"] > 6,
            "window_below_two_heavy_threshold": 4 * p["heavy_mass_squared"] > 6,
            "actual_threshold_amplitude_above_twenty_three_lambda": threshold
            > 23 * ell,
            "actual_positive_subwindow_amplitude_above_forty_two_lambda": positive_window
            > 42 * ell,
            "actual_all_angle_window_amplitude_below_seventy_three_lambda": maximum
            < 73 * ell,
        },
    }
