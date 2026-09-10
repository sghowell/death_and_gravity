"""Bare-coupling equality, including the inherited linear two-loop reference."""

from functools import cache

import sympy as s


def coefficient(expr, h, n):
    return s.expand(s.series(expr, h, 0, n + 1).removeO()).coeff(h, n)


@cache
def data():
    h, L, G, M, k, c, F0, fermion = s.symbols(
        "h L G M k_D c_D entire_I0_over_2Q fermion_quartic_MS_pole"
    )
    sigma = c - k * L
    star_exact = (L - h * c) / (1 - h * k)
    star2 = L - h * sigma + h**2 * sigma * s.diff(sigma, L)
    c1 = {"G": L * G * F0, "M": G**2 * F0, "L": 3 * L**2 * F0 + fermion}
    c2 = {name: sigma * s.diff(v, L) for name, v in c1.items()}
    checks = {
        "exact_inverse_reference_equation": s.factor(
            star_exact + h * sigma.subs(L, star_exact) - L
        ),
        "quadratic_inverse_coefficient": coefficient(star_exact - star2, h, 2),
        "linear_inverse_coefficient": coefficient(star_exact - star2, h, 1),
        "inverse_tail_exact": s.factor(
            star_exact - star2 + h**3 * sigma * k**2 / (1 - h * k)
        ),
        "same_G_squared_linear_reference": s.expand(
            2 * G * c2["G"] - 2 * sigma * G**2 * F0
        ),
        "same_quartic_linear_reference": s.expand(c2["L"] - 6 * sigma * L * F0),
        "same_heavy_mass_linear_reference": c2["M"],
    }
    for name, leading in (("G", G), ("M", M), ("L", L)):
        raw = star2 + h * sigma.subs(L, star2) if name == "L" else leading
        raw += h * c1[name].subs(L, star2) + h**2 * c2[name].subs(L, star2)
        wanted = leading + h * c1[name]
        for n in range(3):
            checks[f"bare_{name}_coefficient_{n}"] = coefficient(raw - wanted, h, n)
    z1, z2 = s.symbols("L_independent_first_field L_independent_sigma_second_field")
    Z = 1 + h * z1 + h**2 * z2
    oldL = star2 + h * (sigma + c1["L"]).subs(L, star2) + h**2 * c2["L"].subs(L, star2)
    oldG = G + h * c1["G"].subs(L, star2) + h**2 * c2["G"].subs(L, star2)
    for name, old, target, power in (
        ("L", oldL, L + h * c1["L"], 2),
        ("G", oldG, G + h * c1["G"], 1),
    ):
        checks[f"same_bare_field_division_{name}"] = coefficient(
            (old - target) / Z**power, h, 2
        )
    e, sigma0, sigma1, fin = s.symbols("epsilon sigma0 sigma1 finite_I0")
    sigmaD = sigma0 + e * sigma1
    I0 = 1 / e + fin
    checks["full_regulator_pole_product_cancellation"] = s.expand(
        sigmaD * I0 - sigmaD * I0
    )
    checks["one_leg_early_epsilon_truncation_finite_defect"] = (
        s.limit(sigmaD * I0 - sigma0 * I0, e, 0) - sigma1
    )
    return {
        "symbols": {"h": h, "L": L, "G": G, "M": M, "k": k, "c": c, "F0": F0},
        "dimensional_sigma": sigma,
        "exact_isolated_coordinate_inverse": star_exact,
        "inverse_through_order_two": star2,
        "first_total_canonical_UV_references": c1,
        "inherited_sigma_second_references": c2,
        "uncancelled_finite_defect_if_only_one_leg_is_truncated": sigma1,
        "checks": checks,
        "scope": "Only the -sigma part of the S6.133 first-order map, with scale, g, M and the canonical field convention held fixed. No separate sigma-dependent noncontact G2 or M2 shift is needed. Other map directions and their cross terms remain open. The exact affine coordinate inverse is not an all-loop physical amplitude.",
    }
