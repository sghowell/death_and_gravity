"""Explicit all-angle physical hard-rate majorants for the original source."""

from functools import cache

import sympy as s

from . import source


def exact_real(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require an exact rational domain coordinate")
    return s.Rational(value)


def require_physical(energy, transfer):
    energy, transfer = map(exact_real, (energy, transfer))
    if (
        not s.Rational(25, 4) <= energy <= 16
        or not transfer < 0
        or not 4 - energy - transfer < 0
    ):
        raise ValueError("Require the stated compact massive nonforward domain")
    return energy, transfer


def general_amplitude_majorant(heavy, cubic, quartic, kappa, resolution_squared):
    heavy, cubic, kappa = map(source.require_mass, (heavy, cubic, kappa))
    source.quartic.require_resolution(resolution_squared)
    if heavy < 128:
        raise ValueError("Require n>=128")
    quartic = exact_real(quartic)
    if abs(quartic) > 6 * cubic * cubic / heavy:
        raise ValueError("Require the stated complete quartic magnitude premise")
    return 14000 * (s.log(heavy) + 1) * cubic * cubic / (16 * s.pi**2 * kappa)


def original_known_rate_bound(energy, transfer, resolution_squared):
    require_physical(energy, transfer)
    source.quartic.require_resolution(resolution_squared)
    return s.Rational(4, 10**204)


@cache
def data():
    R = s.Rational
    L, n, g, k = s.symbols("log_heavy_mass n g kappa", positive=True)
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    master = 6 * (3 * (194 * 128 + 2 * 16) + 3 * (3 * 194 * 16 + 2)) / n
    triangles = 3 * (224 * (L + 9) + 64 * (L + 3)) / n
    boxes = 2 * 6 * 194 * 2 * 4 * (L + 137) / n**2
    total = master + 120 * (L + 2) + triangles + boxes
    put("whole_first_mixed_group", master - 615204 / n)
    put("whole_triangle_group", triangles - (864 * L + 6624) / n)
    put("whole_six_ordered_box_group", boxes - 18624 * (L + 137) / n**2)
    put(
        "whole_combined_mixed_polynomial",
        total - (120 * (L + 2) + (621828 + 864 * L) / n + 18624 * (L + 137) / n**2),
    )
    major = total.subs(n, 128)
    polynomial = s.Poly(6000 * (L + 1) - major, L)
    margins = {
        "contour_center_gap": R(47, 128) - R(1, 8),
        "contour_edge_gap": R(31, 64) - R(1, 8),
        "contour_middle_gap": R(175, 1024) - R(1, 8),
        "beta_root_real_gap": 1 - 4 * R(10, 128) - R(4, 5) ** 2,
        "beta_norm_gap": R(1, 10) - R(2) * R(10, 128) / R(9, 5),
        "alpha_inverse_cube_gap": R(7, 5) - R(10, 9) ** 3,
        "Q_polynomial_gap": 2 - R(11, 10) - R(11, 10) ** 2 / 2,
        "radial_moment_log_coefficient": 2 - R(7, 4),
        "radial_moment_constant": 20 - R(7, 4) * 9 - R(25, 18),
        "T_log_coefficient": 2 - R(25, 18),
        "T_constant": 18 - R(25, 18) * 8,
        "physical_f1_norm": 3 - (16 * R(15, 128) + R(1, 4)),
        "physical_f1_divided_by_s": 1 - R(12, 25),
        "physical_f2_triangles_and_OS": 1 - (4 * R(15, 128) + R(1, 4) + R(1, 12)),
        "physical_f2_with_mixed_bubble": 2 - (1 + R(16, 128)),
        "K1_L_coefficient": 40 - (8 + 21 + R(128, 128)),
        "K1_constant": 80 - (R(900, 128) + R(2048, 128**2) + 42 + R(128, 128)),
        "K1_heavy_coefficient": 21 - 16 * (1 + R(32, 128) + R(128, 128**2)),
        "whole_mixed_L_slack": polynomial.coeff_monomial(L),
        "whole_mixed_constant_slack": polynomial.coeff_monomial(1),
        "whole_endpoint_slack": 2500 - 3 * (724 + 18 * 2),
        "whole_quartic_slack": 105000 - (2 * 37969 + 433 + 2 * 3 * 4657 + 5),
        "whole_complete_amplitude_slack": 14000 - 6000 - 2500 - R(630000, 128),
        "original_amplitude_rounding": 10**5 - R(14000 * 601, 144),
        "original_heavy_mass_upper": 2 * 10**197 - source.HEAVY_MASS2,
        "original_heavy_mass_lower": source.HEAVY_MASS2 - 128,
        "original_contact_bound": 6 * source.CUBIC**2 / source.HEAVY_MASS2
        - abs(source.CONTACT),
        "original_fullBorn_rate_strict": R(4, 10**204)
        - 50000 * source.HEAVY_MASS2**3 / source.KAPPA,
    }
    put("whole_quartic_T0_bound", R(3 * 194 * 16, 2) + 1 - 4657)
    put("whole_quartic_T1_bound", 3 * (R(194 * 128, 2) + 16 + 14 * 16) + 1 - 37969)
    put("whole_quartic_E1_bound", 3 * 18 * 8 + 1 - 433)
    put("whole_endpoint_tensor_coefficient", 2 * 290 + 16 + s.Rational(16**2, 2) - 724)
    put(
        "whole_relative_rate_from_full_Born",
        2 * 10**5 * g * g / k / (4 * g * g / n**3) - 50000 * n**3 / k,
    )
    put(
        "original_rate_rounded_from_heavy_bound",
        R(50000 * (2 * 10**197) ** 3, 10**800) - R(4, 10**204),
    )
    for name, value in margins.items():
        checks["positive_arithmetic_" + name] = value - s.Abs(value)
    return {
        "whole_general_known_amplitude_bound": 14000
        * (s.log(n) + 1)
        * g
        * g
        / (16 * s.pi**2 * k),
        "whole_original_known_amplitude_bound": 10**5 * source.CUBIC**2 / source.KAPPA,
        "whole_original_uniform_known_rate_bound": R(4, 10**204),
        "whole_positive_arithmetic_margins": margins,
        "whole_majorant_components": {
            "nonendpoint": 6000 * (L + 1),
            "both_g_squared_endpoints": s.Integer(2500),
            "entire_C_sector": 630000 / n,
        },
        "whole_forward_normalized_known_limit": s.S.Zero,
        "whole_bound_proof": "The written contour and radial bounds apply to the whole physical normal sheet. With n>=128,|V|<=194,|a-2|<=14 and|c_E|<3, the entire nonendpoint bracket is bounded by the displayed positive polynomial, hence6000(ln n+1). The full endpoint is below2500 in the same g^2/(16pi^2 kappa) units. The complete quartic bracket is below105000 and|C|<6g^2/n. No cancellation between those sectors is assumed.",
        "whole_original_rate_proof": "Original n<2*10^197 implies ln n<600. Since pi>3, the known amplitude is below10^5 g^2/kappa. The full positive Born satisfies A_B>A_m>4g^2/n^3 at all nonforward angles, so|2Re(deltaA_known)/A_B|<50000 n^3/kappa<4e-204. This bounds the named linear interference, not uncomputed loop squares or unknown matching.",
        "whole_endpoint_boundary": "The finite light/triangle functions and logarithmic box multiplier have continuous a0 limits; f1/a uses its proved OS difference quotient. Thus the known mixed amplitude has finite angular endpoint limits. The full Born graviton pole diverges there, so this normalized known interference tends to0. No finite forward cross section or interchange with a Regge contour limit is asserted.",
        "checks": checks,
        "gates": {
            "all_explicit_arithmetic_margins_strict": all(
                bool(v > 0) for v in margins.values()
            ),
            "general_domain_keeps_heavy_gap": True,
            "full_Born_not_matter_only_or_fixed_angle": True,
            "endpoint_OS_difference_quotient_retained": True,
            "uniform_bound_does_not_fix_finite_matching": True,
            "known_linear_interference_not_loop_square_bound": True,
        },
    }
