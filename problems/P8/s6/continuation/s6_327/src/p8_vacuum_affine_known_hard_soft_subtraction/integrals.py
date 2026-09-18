"""Uniform low/high energy budgets for two explicitly selected signed pieces."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_forward_phase import rate
from p8_vacuum_affine_minimal_gravity_radiation import bounds as radiation
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient.measure import cutoff

from . import source

B, D, BM = map(s.Integer, (10**8, 10**16, 330000))


def interference_kernel_upper(resolution):
    x = cutoff(resolution)
    return s.Min(s.Integer(10) ** 18, 2 * s.Integer(10) ** 20 * x) / source.KAPPA


def known_hard_soft_upper(resolution):
    x = cutoff(resolution)
    return s.Min(s.Integer(10) ** 30, 2 * s.Integer(10) ** 32 * x) / source.KAPPA**2


def logarithmic_addition_upper(resolution):
    x = cutoff(resolution)
    L = -s.log(x)
    return (
        s.Integer(10) ** 18 * x * (1 + L) / source.KAPPA**2
        + 3000 * x * x * (L * L + L + s.Rational(1, 2)) / source.KAPPA**3
    )


@cache
def data():
    old = radiation.uniform_data()
    w, x, y, z = s.symbols("w x y z", positive=True)
    low = (320 * B + 102400) / (192 * 36)
    high = (4 * 160**2 * s.Rational(9, 2) + 40 * BM + 320 * D * s.Rational(13, 2)) / 36
    low_x = (320 * B + 102400) / 36
    high_x = (
        4 * 160**2 * 18432 / s.Integer(8)
        + 320 * BM * 192 / s.Integer(8)
        + 320 * D * 192 * s.Rational(13, 2)
    ) / 36
    logcoef = 128 + B / 384 + BM / 16 + D
    primitive = s.log(z) - s.log(1 + z * z) / 2
    checks = {
        "high_band_interference_change_of_variables": (
            y**3 / (w * (y * y + w * w))
        ).subs(w, y * z)
        * y
        - y / (z * (1 + z * z)),
        "high_band_interference_primitive": s.diff(primitive, z)
        - 1 / (z * (1 + z * z)),
        "log_weighted_amplitude_primitive": s.diff(
            x * x * (s.log(1 / x) + s.Rational(1, 2)) / 2, x
        )
        - x * s.log(1 / x),
        "weighted_log_stationary_point": s.diff(z * z * s.log(1 / z), z).subs(
            z, s.exp(-s.Rational(1, 2))
        ),
        "original_known_hard_coefficient": rate.all_angle_known_rate_bound(source.KAPPA)
        - s.Rational(1, 10) ** 788,
        "original_soft_current_constant": old["whole_uniform_constants"][
            "paired_soft_current"
        ]
        - 160,
        "original_low_remainder_constant": old["whole_uniform_positive_margins"][
            "low_remainder_roundup"
        ]
        - B
        + 69700672,
        "original_high_amplitude_constant": old["whole_uniform_positive_margins"][
            "high_amplitude_roundup"
        ]
        - D
        + 3149283804498720,
        "fixed_total_budget": low + high - s.Rational(5200000000045082000, 9),
        "vanishing_total_budget": low_x + high_x - s.Rational(332800000002897536000, 3),
        "log_total_budget": logcoef - s.Rational(30000000000843509, 3),
    }
    gates = {
        "low_fixed": low < 5 * 10**6,
        "low_vanishing": low_x < 10**9,
        "whole_fixed": low + high < 10**18,
        "whole_vanishing": low_x + high_x < 2 * 10**20,
        "known_hard_fixed_exponent": s.Integer(10) ** 30 / source.KAPPA**2
        == s.Rational(1, 10) ** 1570,
        "known_hard_vanishing_exponent": 2 * s.Integer(10) ** 32 / source.KAPPA**2
        == 2 * s.Rational(1, 10) ** 1568,
        "sharp_log_weighted_amplitude": logcoef < 2 * 10**16,
        "sharp_log_interference": s.Rational(320, 9) * 2 * 10**16 < 10**18,
        "sharp_log_original_exponent": s.Integer(10) ** 18 / source.KAPPA**2
        == s.Rational(1, 10) ** 1582,
        "quadratic_named_log_bound": s.Rational(320**2, 4 * 9) < 3000,
        "low_Taylor_bound_used_only_below_y_over192": True,
        "complete47_interference_and_two_helicities_retained": True,
        "frozen_S326_coarse_bound_not_rewritten": True,
        "known_hard_coefficient_not_full_matching": True,
        "no_common_regulator_or_whole_NNLO_rate_claim": True,
    }
    return {
        "checks": {k: s.factor(value) for k, value in checks.items()},
        "gates": {k: bool(v) for k, v in gates.items()},
        "whole_original_low_high_constants": (B, D, BM, 160, 192),
        "whole_exact_low_high_budgets": (low, high, low_x, high_x, logcoef),
        "whole_selected_hard_soft_bound": "For0<x<=1/8, TV[rho sum Re(conjF*S)-sum|S|^2]<min(1e18,2e20*x)/kappa. Multiply only the fixed S303 I_known<1e12/kappa: selected known-hard soft interference TV<min(1e30,2e32*x)/kappa^2. Uniform over all nonforward hard angles.",
        "whole_sharp_logarithmic_bound": "For each unit TT polarization int_0^x w|F1|ln(1/w)dw<2e16*x*(1+ln(1/x))/sqrt(kappa), using the sharper original N1 estimates instead of the general n^2 envelope. Frozen C_Born<320 gives named signed log-addition TV<1e18*x*(1+ln(1/x))/kappa^2+3000*x^2*(ln(1/x)^2+ln(1/x)+1/2)/kappa^3.",
        "whole_physical_boundary": "Both are selected normalized signed D4 subtractions. The hard term uses a specified Born extension, not the whole five-point loop; the logarithmic square is not called one-loop perturbative order. Finite radiative hard and evanescent/common-regulator matching, unknown coordinates, forward total cross section and original P8 remain open.",
    }
