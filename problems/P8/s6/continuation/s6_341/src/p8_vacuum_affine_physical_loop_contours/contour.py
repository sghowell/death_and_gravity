"""Full Feynman-sign homotopy and a joint fixed-contour complex gap."""

from functools import cache

import sympy as s

from . import source

RHO = s.Rational(1, 4096)
LIGHT_GAP = s.Rational(1, 32)
GAP_DIVISOR = s.Integer(64)


def light(x, h=1, timelike=True):
    return x + s.I * h * x * (1 - x) * (1 - 2 * x) if timelike else x


def heavy(t, h=1):
    return t - s.I * h * t * (1 - t)


@cache
def data():
    checks = {}
    x, t, h, n, v, b = s.symbols("x t h n v b", real=True)
    R = s.Rational
    q = x * (1 - x)
    d = 1 - 2 * x
    Q = t * (1 - t)
    xh = x + s.I * h * q * d
    zh = t - s.I * h * Q
    lh = 1 - v * xh * (1 - xh)
    er, el = s.symbols("eR eL", real=True)
    ch = 1 - xh * el - (1 - xh) * er
    checks["deformed_light_product"] = s.expand(
        xh * (1 - xh) - (q + h * h * q * q * d * d + s.I * h * q * d * d)
    )
    checks["light_real_part"] = s.expand(
        s.re(lh) - (1 - v * (q + h * h * q * q * d * d))
    )
    checks["light_negative_imaginary_part"] = s.expand(s.im(lh) + v * h * q * d * d)
    checks["light_derivative_square"] = s.expand(
        abs(s.diff(xh, x)) ** 2 - (1 + h * h * (1 - 6 * q) ** 2)
    )
    checks["heavy_derivative_square"] = s.expand(
        abs(s.diff(zh, t)) ** 2 - (1 + h * h * (1 - 2 * t) ** 2)
    )
    checks["heavy_coordinate_square"] = s.expand(
        abs(zh) ** 2 - t * t * (1 + h * h * (1 - t) ** 2)
    )
    checks["complement_coordinate_square"] = s.expand(
        abs(1 - zh) ** 2 - (1 - t) ** 2 * (1 + h * h * t * t)
    )
    checks["light_coordinate_square"] = s.expand(
        abs(xh) ** 2 - x * x * (1 + h * h * (1 - x) ** 2 * d * d)
    )
    checks["light_complement_square"] = s.expand(
        abs(1 - xh) ** 2 - (1 - x) ** 2 * (1 + h * h * x * x * d * d)
    )
    checks["heavy_complement_squared_real"] = s.expand(
        s.re((1 - zh) ** 2) - (1 - t) ** 2 * (1 - h * h * t * t)
    )
    checks["heavy_complement_squared_imag"] = s.expand(
        s.im((1 - zh) ** 2) - 2 * h * Q * (1 - t)
    )
    checks["mixed_weight_real"] = s.expand(s.re(zh * (1 - zh)) - (Q + h * h * Q * Q))
    checks["mixed_weight_imag"] = s.expand(s.im(zh * (1 - zh)) - h * Q * (2 * t - 1))
    checks["convex_virtuality_real"] = s.expand(s.re(ch) - (1 - x * el - (1 - x) * er))
    checks["virtuality_imaginary_part"] = s.expand(s.im(ch) + h * q * d * (el - er))
    checks["heavy_transfer_imaginary_part"] = s.expand(
        s.im(-b * zh * zh) - 2 * b * t * h * Q
    )
    checks["n_term_imaginary_part"] = s.expand(s.im(n * zh) + n * h * Q)
    lre, lim, cre, cim = s.symbols("Lreal Limag creal cimag", real=True)
    delta_complex = (
        (1 - zh) ** 2 * (lre + s.I * lim)
        + n * zh
        + (cre + s.I * cim) * zh * (1 - zh)
        - b * zh * zh
    )
    imformula = (
        (1 - t) ** 2 * (1 - h * h * t * t) * lim
        + 2 * h * Q * (1 - t) * lre
        - n * h * Q
        + (Q + h * h * Q * Q) * cim
        + h * Q * (2 * t - 1) * cre
        + 2 * b * t * h * Q
    )
    checks["complete_homotopy_imaginary_part"] = s.expand(
        s.im(delta_complex) - imformula
    )
    checks["three_q_ranges_cover_discriminant"] = s.expand(d * d - (1 - 4 * q))

    gates = {}
    gates.update(
        {
            "light_small_q_margin": 1 - 16 * (R(1, 32) + R(1, 32) ** 2) == R(31, 64),
            "light_large_q_margin": 1 - R(45, 8) * R(3, 16) == -R(7, 128),
            "light_middle_imaginary_gap": R(45, 8) * R(1, 32) * R(1, 4) > R(1, 32),
            "small_real_gap": R(31, 64) > R(1, 32),
            "large_real_gap": R(7, 128) > R(1, 32),
            "light_product_modulus": R(1, 4) + R(1, 16) + R(1, 4) == R(9, 16) < 1,
            "light_uniform_modulus": 1 + 16 * R(9, 16) < 20,
            "timelike_homotopy_budget": 2 + R(13, 8) + 8 < 12,
            "crossed_homotopy_budget": 8 + 1 + 8 == 17 and s.Integer(17) < 20,
            "complex_c_bound": 1 + 4 * (2 + R(1, 4096)) < 10,
            "complex_b_bound": (16 + R(1, 4096)) / 4 < 5,
            "small_t_remaining_budget": 120 + 40 + 10 < 180,
            "light_neighborhood_relative_budget": 4 * R(1, 4096) == R(1, 32) / 32,
            "small_t_ell_margin": R(1, 10) - R(1, 32) > R(1, 16),
            "small_t_W_margin": R(1, 10) - R(1, 64) - R(180, 10**6) > R(1, 16),
            "large_t_remaining_budget": 4 * (20 + R(1, 4096)) + 40 + 20 < 260,
            "large_t_mass_margin": R(127, 128) - R(260, 500000) > R(1, 2),
            "large_t_joint_gap": 31 * 500000 > 20,
            "original_mass_lower": source.HEAVY_MASS2 > 10**6,
            "original_mass_upper": source.HEAVY_MASS2 < 10**198,
            "ln10_elementary_bound": sum(
                R(7, 3) ** j / s.factorial(j) for j in range(7)
            )
            > 10,
            "logarithm_500_budget": 200 * R(7, 3) < 500,
            "Cbar_full_budget": 8 * 64 * 500 == 256000,
            "Dbar_full_budget": 16 * 64**2 * 500 == 32768000,
        }
    )

    gates["subthreshold_light_gap"] = 1 - R(4, 3) / 4 == R(2, 3) > R(1, 32)
    return {
        "checks": {key: s.factor(value) for key, value in checks.items()},
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_homotopy": {
            "xi": xh,
            "z": zh,
            "light_factor": lh,
            "imaginary_part": imformula,
        },
        "whole_common_Feynman_sign": "Im Delta<=-(n-20)h t(1-t)<0 for0<h<=1,0<t<1. Timelike budget2+13/8+8<12; subthreshold budget8+1+8=17. Both endpoint faces remain nonzero for h>0.",
        "whole_uniform_light_gap": LIGHT_GAP,
        "whole_uniform_joint_gap": "|Delta|>=(1/32+n0*t)/64 on the fixed final contour. Every invariant/virtuality varies by<=1/4096 and |n-n0|<=n0/128,10^6<=n0<10^198.",
        "whole_continuation_boundary": "Cauchy/Stokes homotopy of the holomorphic parameter top form preserves the negative Feynman prescription; fixed-coordinate side faces vanish. The final integral defines a local continued physical boundary germ. This is NOT a single first-sheet disk across the cut, or an absolute-value real-parameter bound across a pole.",
        "whole_written_gap_proof": "Three light-q ranges give|L0|>1/32. For t<=1/2, |L0+n0*z|>=(|L0|+n0*t)/10; all remaining terms<=180t+4rho+n0*t/64 leave>(|L0|+n0*t)/16. For t>=1/2, |nz|>=(127/128)n0*t and the rest<260. The common1/64 lower bound follows.",
    }
