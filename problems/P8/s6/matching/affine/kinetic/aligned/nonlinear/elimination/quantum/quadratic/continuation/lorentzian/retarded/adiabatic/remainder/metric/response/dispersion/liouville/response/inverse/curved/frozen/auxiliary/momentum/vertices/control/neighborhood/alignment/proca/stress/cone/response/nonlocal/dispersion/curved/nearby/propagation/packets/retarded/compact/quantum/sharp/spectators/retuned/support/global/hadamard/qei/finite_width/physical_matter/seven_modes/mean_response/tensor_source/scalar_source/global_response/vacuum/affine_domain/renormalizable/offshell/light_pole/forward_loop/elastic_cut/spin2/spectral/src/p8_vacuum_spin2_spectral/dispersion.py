"""Positive pushforward measure and exact once-subtracted vertex dispersion."""

from functools import cache

import sympy as sp
from p8_vacuum_spin2 import triangles


@cache
def data():
    d = triangles.data()
    x, z, t, M, g = d["x"], d["z"], d["t"], d["M"], d["g"]
    d1 = d["on_shell_two_point_denominator"]
    hL = (1 - x) ** 2 * z * (1 - z)
    hH = x * x * z * (1 - z)
    wL = x * x * (1 - x)
    wH = x * (1 - x) ** 2
    tau, weight, h, den = sp.symbols(
        "spectral_invariant positive_vertex_weight transfer_weight positive_denominator",
        positive=True,
    )
    generic = weight / (den - h * t) - weight / den
    measure_weight = weight / h
    kernel = t * measure_weight / (tau * (tau - t))
    pref = g / (16 * sp.pi**2)
    return {
        "light_spectral_location": d1 / hL,
        "heavy_spectral_location": d1 / hH,
        "light_pushforward_weight": pref * wL / hL,
        "heavy_pushforward_weight": pref * wH / hH,
        "once_subtracted_spectral_kernel": kernel,
        "slope_moment_kernel": measure_weight / tau**2,
        "support": "Light stress cut [4,infinity); heavy stress cut [4M,infinity). Endpoints are support limits, not finite delta atoms.",
        "measure_boundary": "The positive pushforward measure need not have finite total mass. Its inverse-first and inverse-second moments are finite. Tonelli and the actual parameter majorants justify the subtracted identity and slope.",
        "scope": "Exact dispersion representation for this one-matter-loop stress form factor only. Not a scattering-amplitude contour, UV completion or gravitational positivity sum rule.",
        "checks": {
            "literal_parameter_to_subtracted_Stieltjes_kernel": sp.factor(
                generic - kernel.subs(tau, den / h)
            ),
            "positive_slope_inverse_second_moment": sp.factor(
                sp.diff(kernel, t).subs(t, 0) - measure_weight / tau**2
            ),
            "light_threshold_positive_gap": sp.factor(
                d1 - 4 * hL - (x * M + (1 - x) ** 2 * (2 * z - 1) ** 2)
            ),
            "heavy_threshold_positive_gap": sp.factor(
                d1
                - 4 * M * hH
                - (M * x * (1 - x) + (1 - x) ** 2 + M * x * x * (2 * z - 1) ** 2)
            ),
            "actual_parent_subtracted_vertex_integrand": sp.factor(
                d["once_Ward_subtracted_unit_square_integrand"]
                - pref * (wL / (d1 - hL * t) + wH / (d1 - hH * t) - (wL + wH) / d1)
            ),
            "once_subtracted_kernel_zero_anchor": kernel.subs(t, 0),
        },
    }
