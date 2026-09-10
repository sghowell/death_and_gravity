"""Both unequal-mass triangle cuts and an anchored universal Q2 primitive."""

from functools import cache

import sympy as sp


@cache
def data():
    T, M, g = sp.symbols(
        "physical_transfer heavy_mass_squared cubic_squared", positive=True
    )
    chi, mA, mB = sp.symbols(
        "single_parameter stress_mass_squared other_mass_squared",
        real=True,
    )
    d = (1 - chi) * mA + chi * mB - chi * (1 - chi)
    discriminant = T * (1 - chi) ** 2 - 4 * d
    general = (T - 4) * chi**2 - 2 * (T - 2 * mA + 2 * mB - 2) * chi + T - 4 * mA
    Z, w, y = sp.symbols(
        "Legendre_argument shifted_parameter rescaled_parameter", positive=True
    )
    rad = sp.sqrt(w * w - Z * Z + 1)
    coefficient = (3 * Z * Z - 1) / 2
    primitive = (
        coefficient * sp.log((w + rad) / sp.sqrt(Z * Z - 1)) + (w / 2 - 2 * Z) * rad
    )
    Q2 = (3 * Z * Z - 1) * sp.log((Z + 1) / (Z - 1)) / 4 - 3 * Z / 2
    A = T - 4
    CH = T - 4 * M
    ZL = 1 + 2 * M / A
    ZH = (T - 2 * M) / sp.sqrt(A * CH)
    rhoL = g * Q2.subs(Z, ZL) / (8 * sp.pi * sp.sqrt(T * A))
    rhoH = g * CH * Q2.subs(Z, ZH) / (8 * sp.pi * sp.sqrt(T) * A ** sp.Rational(3, 2))
    C, B, Ag = sp.symbols(
        "positive_cut_gap positive_middle_half positive_external_gap", positive=True
    )
    beta = sp.Symbol("positive_pair_discriminant", positive=True)
    roots = ((1 - beta) / 2, (1 + beta) / 2)
    return {
        "T": T,
        "M": M,
        "g": g,
        "Z": Z,
        "general_cut_discriminant": general,
        "universal_anchored_cut_primitive": primitive,
        "universal_positive_Q2": Q2,
        "light_pair_Q2_argument": ZL,
        "heavy_pair_Q2_argument": ZH,
        "light_pair_imaginary_form_factor_above_threshold": rhoL,
        "heavy_pair_imaginary_form_factor_above_threshold": rhoH,
        "universal_cut_moment": "integral_0^(Z-sqrt(Z^2-1)) y^2/sqrt(y^2-2Zy+1) dy = Q2(Z), Z>1",
        "threshold_prescription": "Each displayed one-loop density is used only strictly above its own threshold and equals zero below it. At threshold its continuous value is zero. Heavy lines have their zeroth-order mass: the heavy particle decays in the interacting theory, so this is not an exact stable-heavy spectral threshold or a uniform higher-loop threshold approximation.",
        "checks": {
            "literal_general_cut_discriminant": sp.expand(discriminant - general),
            "physical_cut_first_root_inside_parameter_interval": sp.expand(
                general.subs(chi, 1) + 4 * mB
            ),
            "light_cut_polynomial_assignment": sp.factor(
                general.subs({mA: 1, mB: M}) - A * (chi * chi - 2 * ZL * chi + 1)
            ),
            "heavy_cut_polynomial_assignment": sp.expand(
                general.subs({mA: M, mB: 1})
                - (A * chi * chi - 2 * (T - 2 * M) * chi + CH)
            ),
            "universal_cut_rescaling": sp.simplify(
                (Ag * chi * chi - 2 * B * chi + C).subs(chi, sp.sqrt(C / Ag) * y)
                - C * (y * y - 2 * B / sp.sqrt(Ag * C) * y + 1)
            ),
            "anchored_Q2_primitive_derivative": sp.simplify(
                sp.diff(primitive, w) - (Z - w) ** 2 / rad
            ),
            "primitive_zero_at_lower_root": sp.simplify(
                primitive.subs(w, sp.sqrt(Z * Z - 1))
            ),
            "primitive_upper_endpoint_before_positive_log_identity": sp.simplify(
                primitive.subs(w, Z)
                - coefficient * sp.log((Z + 1) / sp.sqrt(Z * Z - 1))
                + 3 * Z / 2
            ),
            "positive_log_identity_squared_argument": sp.factor(
                ((Z + 1) / sp.sqrt(Z * Z - 1)) ** 2 - (Z + 1) / (Z - 1)
            ),
            "heavy_cut_Legendre_argument_above_one_gap": sp.factor(
                (T - 2 * M) ** 2 - A * CH - 4 * (T + M * M - 4 * M)
            ),
            "pair_delta_root_first": sp.factor(
                1 - 4 * roots[0] * (1 - roots[0]) - beta * beta
            ),
            "pair_delta_root_second": sp.factor(
                1 - 4 * roots[1] * (1 - roots[1]) - beta * beta
            ),
            "two_pair_roots_delta_Jacobian_sum": 2 / ((1 - chi) ** 2 * T * beta)
            - sum(1 / ((1 - chi) ** 2 * T * beta) for root in roots),
        },
    }
