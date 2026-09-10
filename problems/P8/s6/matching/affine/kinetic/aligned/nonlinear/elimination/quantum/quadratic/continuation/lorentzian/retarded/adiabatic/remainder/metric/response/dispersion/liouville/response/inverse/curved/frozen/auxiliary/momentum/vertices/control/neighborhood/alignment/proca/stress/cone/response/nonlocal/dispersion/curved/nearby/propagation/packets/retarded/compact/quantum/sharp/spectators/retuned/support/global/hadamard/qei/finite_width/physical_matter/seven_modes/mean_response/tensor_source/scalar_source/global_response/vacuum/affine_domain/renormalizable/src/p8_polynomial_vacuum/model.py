"""Stable polynomial two-scalar realization of the exact massive tree exchange."""

from functools import cache

import sympy as sp
from p8_exceptional_vacuum import analytic
from p8_exceptional_vacuum import heavy as previous
from p8_uv import vacuum


@cache
def data():
    lam, gamma = sp.symbols("positive_fixed_lambda positive_fixed_gamma", positive=True)
    D = 2 * lam / gamma
    mu2 = D + 2
    g2 = 8 * lam * D
    G2 = g2 * D**2 / 4
    quartic = g2 * (3 * D - 2) / 4
    s, t, w = vacuum.s, vacuum.transfer, vacuum.w
    amplitude = -quartic + G2 * sum(1 / (mu2 - z) for z in (s, t, w))
    old = previous.data()["exact_crossing_symmetric_tree_exchange"]
    stable = sp.factor(quartic - 3 * G2 / mu2)
    Phi, H = sp.symbols("canonical_light_field canonical_heavy_field", real=True)
    G = sp.Symbol("positive_cubic_coupling", positive=True)
    V = Phi**2 / 2 + mu2 * H**2 / 2 + G * H * Phi**2 / 2 + quartic * Phi**4 / 24
    square = (
        mu2 * (H + G * Phi**2 / (2 * mu2)) ** 2 / 2
        + Phi**2 / 2
        + (quartic - 3 * G * G / mu2) * Phi**4 / 24
    )
    sub = {lam: analytic.VACUUM_LAMBDA_BAR, gamma: analytic.FIXED_GAMMA}
    actual = {
        name: sp.factor(value.subs(sub))
        for name, value in {
            "lambda": lam,
            "gamma": gamma,
            "D": D,
            "heavy_mass_squared": mu2,
            "cubic_coupling_squared": G2,
            "bare_polynomial_quartic": quartic,
            "positive_completed_square_quartic_margin": stable,
            "relative_heavy_tree_width_upper": G2 / (96 * mu2),
        }.items()
    }
    return {
        "light_mass_squared": sp.Integer(1),
        "D": D,
        "heavy_mass_squared": mu2,
        "cubic_coupling_squared": G2,
        "polynomial_quartic": quartic,
        "canonical_potential": V,
        "completed_square_potential": square,
        "positive_quartic_margin": stable,
        "exact_tree_amplitude": amplitude,
        "actual_parameters": actual,
        "coupling_relation": "G=+sqrt(G_squared); all fields canonical",
        "width_convention": "Gamma(H->Phi Phi)=G_squared sqrt(1-4/mu_squared)/(32 pi mu); bound uses pi>3",
        "scope": "Power-counting renormalizable massive flat-vacuum model; exact on-shell tree four-point match, not all higher-point equivalence, a nonperturbative all-energy UV theorem, gravity or a bounce parent",
        "checks": {
            "full_on_shell_tree_amplitude_matches_previous_exchange": sp.factor(
                (amplitude - old).subs(w, 4 - s - t)
            ),
            "completed_square_is_literal_potential": sp.factor(V - square),
            "positive_quartic_margin_factorization": sp.factor(
                stable - 4 * G2 * (D - 1) / (D**2 * (D + 2))
            ),
            "exact_channel_polynomial_division": sp.factor(
                (s - 2) ** 2 / (D - (s - 2)) - D**2 / (D - (s - 2)) + D + s - 2
            ),
            "canonical_vacuum_mass_matrix": sp.hessian(V, (Phi, H)).subs({Phi: 0, H: 0})
            - sp.diag(1, mu2),
            "tree_origin_is_stationary": sp.Matrix(
                [sp.diff(V, Phi), sp.diff(V, H)]
            ).subs({Phi: 0, H: 0}),
            "eliminated_stationarity_has_unique_real_zero_for_positive_margin": sp.factor(
                sp.diff(V, Phi).subs(H, -G * Phi**2 / (2 * mu2)).subs(G, sp.sqrt(G2))
                - Phi * (1 + stable * Phi**2 / 6)
            ),
            "exact_forward_coefficient": sp.factor(
                vacuum.forward_coefficient(amplitude).subs(vacuum.mass, 1) - 4 * lam
            ),
            "stable_margin_uses_same_cubic": sp.factor(
                (quartic - 3 * G * G / mu2).subs(G, sp.sqrt(G2)) - stable
            ),
        },
        "bounds": {
            "actual_D_above_one": actual["D"] > 1,
            "global_completed_square_margin_positive": actual[
                "positive_completed_square_quartic_margin"
            ]
            > 0,
            "canonical_heavy_mass_above_two_light_threshold": actual[
                "heavy_mass_squared"
            ]
            > 4,
            "cubic_squared_in_light_mass_units_below_one_e_minus_seven": actual[
                "cubic_coupling_squared"
            ]
            < sp.Rational(1, 10**7),
            "quartic_coupling_below_one_e_minus_204": actual["bare_polynomial_quartic"]
            < sp.Rational(1, 10**204),
            "relative_tree_width_upper_below_one_e_minus_207": actual[
                "relative_heavy_tree_width_upper"
            ]
            < sp.Rational(1, 10**207),
        },
    }
