"""Validated rational-parameter enclosure of the complete paired mixed quartic row."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_heavy_vertices import calibration as earlier
from p8_vacuum_fermion_self_energy_chord.tail import rational


def dyadic(n):
    n = rational(n)
    if n < 1:
        raise ValueError("Need dyadic argument >=1")
    j = max(0, int(n.p).bit_length() - int(n.q).bit_length())
    while sp.Integer(2) ** j < n:
        j += 1
    while j > 0 and sp.Integer(2) ** (j - 1) >= n:
        j -= 1
    return j


def enclosure(m, Y, L, g, M, Q):
    m, Y, L, g, M, Q = map(rational, (m, Y, L, g, M, Q))
    if m < 720 or min(Y, L, g) < 0 or M < 10000 or not 0 < Q <= 144:
        raise ValueError("Need m>=720, nonnegative Y,L,g, M>=10000, 0<Q<=144")
    T = 4 * m * m
    nm, nM, nT = dyadic(m * m), dyadic(2 * M), dyadic(T)
    Vc = L + g / (M - 3)
    V = Vc + 4 * g / M
    pref = 6 * Y * Y / Q**2
    poly = 24 * nM * nM + 184 * nM + 184
    delta = 2123366400000 * pref * V / sp.sqrt(m)
    bubble = 3 * 3432 * pref * Vc
    triangles = 3 * 64 * pref * g * poly / M
    F = 2 * pref * (46 + 32 * nm + (300 + 100 * nT) / T)
    tree = 2 * g / (M - 2) ** 3
    conversion = tree * F
    return {
        "mF": m,
        "N": 6,
        "Y_upper": Y,
        "L": L,
        "g": g,
        "M": M,
        "Q_lower": Q,
        "spectral_threshold": T,
        "dyadic_m_squared": nm,
        "dyadic_twice_M": nM,
        "dyadic_T": nT,
        "full_vertex_upper": V,
        "external_vertex_upper": Vc,
        "zero_soft_triangle_log_polynomial": poly,
        "nonzero_soft_remainder_upper": delta,
        "subtracted_zero_soft_bubble_upper": bubble,
        "zero_soft_heavy_triangles_upper": triangles,
        "finite_MS_reference_absolute_upper": F,
        "finite_MS_reference_conversion_b2_upper": conversion,
        "tree_b2": tree,
        "complete_paired_mixed_b2_upper": delta + bubble + triangles + conversion,
        "scope": "Only scalar_Phi4_W1_F2 including its proper quartic MS forest and finite outer MS conversion; not other order-two matching or canonical terms.",
    }


@cache
def data():
    p = earlier.data()["actual_reference_parameters"]
    r = enclosure(p["mF"], p["Y_upper"], p["L"], p["g"], p["M"], 144)
    E = r["complete_paired_mixed_b2_upper"]
    return {
        "actual_reference_parameters": p,
        "enclosure": r,
        "relative_upper": E / (4 * p["lambda"]),
        "checks": {
            "same_fermion_mass": r["mF"] - 10**200,
            "same_actual_tree": r["tree_b2"] - 4 * p["lambda"],
            "all_disjoint_pieces_added_once": E
            - sum(
                r[k]
                for k in (
                    "nonzero_soft_remainder_upper",
                    "subtracted_zero_soft_bubble_upper",
                    "zero_soft_heavy_triangles_upper",
                    "finite_MS_reference_conversion_b2_upper",
                )
            ),
            "dimensionless_MS_leading_conversion_factor": r[
                "finite_MS_reference_absolute_upper"
            ]
            - 12
            * p["Y_upper"] ** 2
            / 144**2
            * (
                46
                + 32 * r["dyadic_m_squared"]
                + (300 + 100 * r["dyadic_T"]) / r["spectral_threshold"]
            ),
            "actual_threshold_dyadic": r["dyadic_T"] - 1331,
            "light_bubble_radial_bound": 33
            * (136 * sp.Rational(1, 2) + 48 * sp.Rational(3, 4))
            - 3432,
            "triangle_polynomial_quadratic": 136 + 48 - 184,
            "triangle_polynomial_constant": 136 + 48 - 184,
            "triangle_polynomial_log_squared": 48 / sp.Integer(2) - 24,
        },
        "bounds": {
            "actual_mass_in_domain": bool(p["mF"] >= 720),
            "actual_heavy_mass_squared_in_domain": bool(p["M"] >= 10000),
            "actual_all_terms_strictly_positive": all(
                r[k] > 0
                for k in (
                    "nonzero_soft_remainder_upper",
                    "subtracted_zero_soft_bubble_upper",
                    "zero_soft_heavy_triangles_upper",
                    "finite_MS_reference_conversion_b2_upper",
                )
            ),
            "actual_momentum_remainder_below_one_e_minus_700": bool(
                r["nonzero_soft_remainder_upper"] < sp.Rational(1, 10**700)
            ),
            "actual_MS_reference_below_one_e_minus_400": bool(
                r["finite_MS_reference_absolute_upper"] < sp.Rational(1, 10**400)
            ),
            "actual_total_below_one_e_minus_609": bool(E < sp.Rational(1, 10**609)),
            "actual_relative_below_one_e_minus_10": bool(
                E / (4 * p["lambda"]) < sp.Rational(1, 10**10)
            ),
            "dyadic_bounds_cover_all_log_arguments": bool(
                p["mF"] ** 2 <= 2 ** r["dyadic_m_squared"]
                and 2 * p["M"] <= 2 ** r["dyadic_twice_M"]
                and r["spectral_threshold"] <= 2 ** r["dyadic_T"]
            ),
        },
    }
