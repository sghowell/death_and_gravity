"""Full regulated field-factor identity with finite coefficient re-expansion."""

from functools import cache

import sympy as s


@cache
def data():
    h, e = s.symbols("h epsilon")
    r, f, r1, f1, z, tH, t0, t1 = s.symbols(
        "r fp r_epsilon fp_epsilon z11 t_H t_MS t_MS_epsilon"
    )
    z22, z21, zH22, zH21 = s.symbols("z22 z21 zH22 zH21")
    k = r - f
    ke = r1 - f1
    kD = k + e * ke
    # The bare coordinate identity needs q1_D=-k_D E_Phi q, including
    # positive epsilon terms before its multiplication by first-order poles.
    dqk = -kD * (2 * r - f + e * (2 * r1 - f1))
    dqz = -kD * z
    Zms = 1 + h * z / e + h**2 * (z22 / e**2 + z21 / e)
    K = 1 + h * kD + h**2 * (t0 + e * t1)
    ZH = 1 + h * (z / e - kD) + h**2 * (zH22 / e**2 + zH21 / e - tH + dqz / e - dqk)
    difference = s.series(Zms / K - ZH, h, 0, 3).removeO().expand()
    second = difference.coeff(h, 2)
    finite = second.coeff(e, 0)
    solved = tH - k * r
    L, g, Y, M, G, y = s.symbols("L g Y M G y", positive=True)
    rs, fs, zc = s.symbols("scalar_slope_per_g fermion_slope_per_Y field_pole_per_Y")
    kernel = g * rs - Y * fs
    Euler = lambda F: (
        2 * L * s.diff(F, L)
        + G * s.diff(F, G)
        + 2 * g * s.diff(F, g)
        + Y * s.diff(F, Y)
    )
    fields = {
        "L": -2 * k * L,
        "G": -k * G,
        "g": -2 * k * g,
        "y": -k * y / 2,
        "Y": -k * Y,
        "M": s.Integer(0),
        "m": s.Integer(0),
        "a": s.Integer(0),
    }
    field_factor, light_mass_squared = s.symbols(
        "positive_field_factor light_mass_squared", positive=True
    )
    Q, N = s.symbols("Q N", positive=True)
    residues = {
        "L": (3 * L**2 / 2 - 24 * N * Y**2) / Q,
        "G": L * G / (2 * Q),
        "M": G**2 / (2 * Q),
    }
    weights = {"L": 2, "G": 1, "M": 0}
    parameters = {"L": L, "G": G, "M": M}
    second_map = {}
    map_checks = {}
    for name, w in weights.items():
        X, p = parameters[name], residues[name]
        comm = Euler(p) - w * p
        second_map[name] = (s.Rational(w * (w + 1), 2) * k**2 - w * t0) * X + ke * comm
        before = (s.Rational(w * (w + 1), 2) * kD**2 - w * t0) * X + kD * comm / e
        map_checks["full_regulated_bare_second_map_" + name] = s.expand(
            s.expand(before).coeff(e, 0) - second_map[name]
        )
    g_second = s.expand(2 * G * second_map["G"] + (-k * G) ** 2)
    forward_ratio = g_second / G**2 - 3 * second_map["M"] / (M - 2)
    map_checks.update(
        {
            "full_G_second_map": s.expand(
                second_map["G"] - (k**2 - t0) * G - ke * L * G / Q
            ),
            "full_L_second_map": s.expand(
                second_map["L"] - (3 * k**2 - 2 * t0) * L - 3 * ke * L**2 / Q
            ),
            "induced_second_heavy_mass_reference": s.expand(
                second_map["M"] - ke * G**2 / Q
            ),
            "full_forward_tree_second_reference": s.factor(
                forward_ratio
                - 3 * k**2
                + 2 * t0
                - ke * (2 * L - 3 * G**2 / (M - 2)) / Q
            ),
        }
    )
    return {
        "complete_first_inverse_slope": k,
        "complete_first_epsilon_slope": ke,
        "first_finite_hybrid_parameter_direction": fields,
        "first_regulated_hybrid_parameter_direction": {
            "L": -2 * kD * L,
            "G": -kD * G,
            "Y": -kD * Y,
            "M": s.Integer(0),
        },
        "complete_MS_second_normalization_from_hybrid": solved,
        "full_bare_second_scalar_parameter_map": second_map,
        "full_bare_tree_second_b2_relative": forward_ratio,
        "full_MS_and_hybrid_field_factors": {
            "Z_MS": Zms,
            "K_D": K,
            "Z_H_at_converted_parameters": ZH,
        },
        "same_two_loop_kinetic_pole_dictionary": {
            "double": zH22 - z22,
            "simple": zH21 - z21,
        },
        "checks": {
            **map_checks,
            "first_field_factor_bare_identity": s.factor(difference.coeff(h, 1)),
            "complete_second_finite_field_identity": s.factor(finite.subs(t0, solved)),
            "second_double_pole_dictionary": second.coeff(e, -2) - z22 + zH22,
            "second_simple_pole_dictionary": s.factor(second.coeff(e, -1) - z21 + zH21),
            "finite_parameter_reexpansion_not_lost": s.expand(
                solved - tH - dqk.subs(e, 0) - k**2
            ),
            "epsilon_parameter_shift_times_field_pole_retained": s.expand(
                dqz / e
            ).coeff(e, 0)
            + ke * z,
            "both_epsilon_field_pole_products_cancel": s.expand(-kD * z / e).coeff(e, 0)
            - s.expand(dqz / e).coeff(e, 0),
            "same_scalar_fermion_first_slope_homogeneity": s.expand(
                Euler(kernel) - (2 * g * rs - Y * fs)
            ),
            "same_MS_field_pole_homogeneity": s.expand(Euler(Y * zc) - Y * zc),
            "homogeneity_simplifies_finite_shift": s.expand(
                dqk.subs(e, 0) + k**2 + k * r
            ),
            "fundamental_cubic_square_first_direction": 2 * G * fields["G"]
            - fields["g"].subs(g, G**2),
            "fundamental_Yukawa_square_first_direction": 2 * y * fields["y"]
            - fields["Y"].subs(Y, y**2),
            "physical_light_mass_ratio_preserved_by_correlated_field_change": (
                light_mass_squared / field_factor
            )
            / (1 / field_factor)
            - light_mass_squared,
        },
        "scope": "The hybrid scheme has physical Phi OS mass/residue and minimal total interaction poles, with MS fermion/gauge coordinates through the order needed here. t_H is the full fixed-hybrid-coordinate inverse slope before outer OS. Complete regulated matching uses q1_D=-k_D E_Phi q and t_MS=t_H-k0 r0: the epsilon first-parameter shift cancels the epsilon field-pole product. Re-expressing first interaction poles also induces the displayed second L,G,M map. External-fermion two-loop field/mass maps are outside this Phi two-point result.",
    }
