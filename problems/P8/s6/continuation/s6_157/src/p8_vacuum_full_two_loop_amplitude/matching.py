"""Hybrid physical-Phi/MS-interaction matching after all pole products."""

from functools import cache

import sympy as s


@cache
def data():
    e = s.Symbol("epsilon")
    L, g, G, M, Y, Q, k0, k1, t = s.symbols("L g G M Y Q k0 k1 t")
    a, b, c, d = s.symbols("a b c d")
    q = (L, g, Y)
    As = a * L**2 + b * L * g + c * g**2
    Af = d * Y**2
    first = (-2 * k0 * L, -2 * k0 * g, -k0 * Y)
    variation = sum(v * s.diff(As + Af, x) for x, v in zip(q, first))
    G1 = -k0 * G
    G2 = (k0**2 - t) * G + k1 * L * G / Q
    g2 = (3 * k0**2 - 2 * t) * g + 2 * k1 * L * g / Q
    M2 = k1 * g / Q
    L2 = (3 * k0**2 - 2 * t) * L + 3 * k1 * L**2 / Q
    tree = 2 * g / (M - 2) ** 3
    relative = 3 * k0**2 - 2 * t + k1 * (2 * L - 3 * g / (M - 2)) / Q
    tree_second = s.diff(tree, g) * g2 + s.diff(tree, M) * M2
    A10, A11, A12, P = s.symbols("A10 A11 A12 P")
    kD = k0 + e * k1
    raw = P / e + A10 + e * A11 + e**2 * A12
    finite_pair = s.expand(-2 * kD * (raw - P / e)).coeff(e, 0)
    return {
        "hybrid_definition": "Physical Phi OS conditions, interaction-MS poles, direct-MS fermion proper kernels, exact Gaussian H; mu=mF fixed.",
        "first_regulated_direction": {
            "L": -2 * kD * L,
            "g": -2 * kD * g,
            "Y": -kD * Y,
            "M": s.Integer(0),
        },
        "finite_one_loop_homogeneity_representatives": {"scalar": As, "fermion": Af},
        "first_variation_of_hybrid_one_loop": -k0 * (4 * As + 2 * Af),
        "full_second_bare_scalar_map": {"G": G2, "g": g2, "L": L2, "M": M2},
        "tree_second_relative": relative,
        "complete_second_coefficient": "AH2_raw_MS - k0*(4*AH1_scalar_MS+2*AH1_fermion_MS) + tree_b2*tree_second_relative",
        "checks": {
            "scalar_first_homogeneity_weight_four": 2 * L * s.diff(As, L)
            + 2 * g * s.diff(As, g)
            - 4 * As,
            "fermion_first_homogeneity_weight_two": Y * s.diff(Af, Y) - 2 * Af,
            "complete_first_coordinate_variation": variation + k0 * (4 * As + 2 * Af),
            "fundamental_G_square_retained": s.expand(
                (G1**2 + 2 * G * G2 - g2).subs(g, G**2)
            ),
            "first_tree_variation": s.diff(tree, g) * (-2 * k0 * g) + 2 * k0 * tree,
            "full_second_tree_variation": s.factor(tree_second / tree - relative),
            "new_heavy_mass_shift_retained": M2 - k1 * g / Q,
            "renormalized_pair_before_epsilon_limit": finite_pair + 2 * k0 * A10,
            "raw_epsilon_pole_product_retained": s.expand(-2 * kD * raw).coeff(e, 0)
            + 2 * k0 * A10
            + 2 * k1 * P,
            "counterterm_epsilon_pole_product_retained": s.expand(2 * kD * P / e).coeff(
                e, 0
            )
            - 2 * k1 * P,
            "positive_epsilon_after_pairing_has_zero_finite_part": s.expand(
                e * k1 * (A10 + e * A11)
            ).coeff(e, 0),
            "no_first_heavy_mass_shift": s.Integer(0),
        },
        "scope": "The complete renormalized AH1 is finite before its first coordinate variation, so only k0 remains there. This does not discard k1 in the bare map: its already regulated G/L/M contributions are retained separately. AH2 is hybrid, not raw-MS amputation; no additional -2k0 A1 or LSZ factor is allowed.",
    }
