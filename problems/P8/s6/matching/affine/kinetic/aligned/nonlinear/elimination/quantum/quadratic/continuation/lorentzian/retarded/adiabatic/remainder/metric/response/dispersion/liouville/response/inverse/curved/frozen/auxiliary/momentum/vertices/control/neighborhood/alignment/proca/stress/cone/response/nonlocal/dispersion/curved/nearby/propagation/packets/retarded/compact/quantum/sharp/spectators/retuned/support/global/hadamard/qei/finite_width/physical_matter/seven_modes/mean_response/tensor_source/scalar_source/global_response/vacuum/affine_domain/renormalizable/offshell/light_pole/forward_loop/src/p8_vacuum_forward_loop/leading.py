"""Renormalized angular-zero radial loop and its explicit second-derivative bound."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model


@cache
def data():
    s = sp.Symbol("forward_s", real=True)
    a, y, M, c, lam = sp.symbols(
        "parameter_a radial_y M small_c quartic_L", positive=True
    )
    g = c * M * M
    ell = 1 - a * s
    A = y + M - ell
    L = y + ell
    C = -lam + g / (M - s)
    V = C + 2 * g / A
    delta = lam - 3 * g / M
    alternative = -delta + g * s / (M * (M - s)) - 2 * g * (y - ell) / (M * A)
    vp = g / (M - s) ** 2 - 2 * g * a / A**2
    vpp = 2 * g / (M - s) ** 3 + 4 * g * a * a / A**3
    derivative = (
        2 * (vp * vp + V * vpp) / L**2
        + 8 * a * V * vp / L**3
        + 6 * a * a * V * V / L**4
    )
    f = V * V / L**2 - C * C / (y + 1) ** 2
    h1 = 4 * g * C / (A * L**2)
    h2 = 4 * g * g / (A * A * L * L)
    hd = C * C * (1 / L**2 - 1 / (y + 1) ** 2)
    T = 1 / (A * L * L)
    Tp = -a / (A * A * L * L) + 2 * a / (A * L**3)
    Tpp = (
        2 * a * a / (A**3 * L * L) - 4 * a * a / (A * A * L**3) + 6 * a * a / (A * L**4)
    )
    U = 1 / (A * A * L * L)
    Upp = (
        6 * a * a / (A**4 * L * L)
        - 8 * a * a / (A**3 * L**3)
        + 6 * a * a / (A * A * L**4)
    )
    gap = 1 / L**2 - 1 / (y + 1) ** 2
    low_log = sp.Integer(72) + 960 + 2400 + 192
    low_const = sp.Integer(2880)
    high_first = sp.Integer(64) + sp.Rational(24, 2) + sp.Rational(15, 3)
    high_second = sp.Rational(5, 4)
    high_third = sp.Integer(288) + 40 + sp.Rational(75, 16)
    high_total = high_first + high_second + high_third
    exp_three_lower = sum(sp.Rational(3) ** k / sp.factorial(k) for k in range(4))
    q = y + sp.Rational(1, 4)
    p = model.data()["actual_parameters"]
    Mp = p["heavy_mass_squared"]
    gp = p["cubic_coupling_squared"]
    cp = gp / Mp**2
    return {
        "s": s,
        "a": a,
        "radial_y": y,
        "M": M,
        "small_c": c,
        "radial_C": C,
        "radial_V0": V,
        "renormalized_angular_zero_integrand_without_radial_measure": f,
        "exact_UV_cancellation_terms": [h1, h2, hd],
        "low_region_log_coefficient": low_log,
        "low_region_constant": low_const,
        "high_region_constant": high_total,
        "all_channel_angular_zero_b2_upper": "13 c^2 [1+log(1+4M)] < 14000 c^2",
        "actual_angular_zero_b2_upper": 14000 * cp * cp,
        "checks": {
            "exact_low_energy_vertex_cancellation": sp.factor(V - alternative),
            "actual_first_radial_vertex_derivative": sp.factor(sp.diff(V, s) - vp),
            "actual_second_radial_vertex_derivative": sp.factor(sp.diff(V, s, 2) - vpp),
            "literal_second_derivative_of_unsubtracted_radial_integrand": sp.factor(
                sp.diff(V * V / L**2, s, 2) - derivative
            ),
            "exact_high_momentum_UV_cancellation_rewrite": sp.factor(f - h1 - h2 - hd),
            "first_high_momentum_kernel_derivative": sp.factor(sp.diff(T, s) - Tp),
            "second_high_momentum_kernel_derivative": sp.factor(sp.diff(T, s, 2) - Tpp),
            "second_squared_heavy_kernel_derivative": sp.factor(sp.diff(U, s, 2) - Upp),
            "first_denominator_difference_derivative": sp.factor(
                sp.diff(gap, s) - 2 * a / L**3
            ),
            "second_denominator_difference_derivative": sp.factor(
                sp.diff(gap, s, 2) - 6 * a * a / L**4
            ),
            "low_region_ratio_positive_gap": sp.factor(
                4 - (y + 1) / q - sp.Rational(3) * y / q
            ),
            "low_region_radial_ratio_positive_gap": sp.factor(
                1 - y / q - sp.Rational(1, 4) / q
            ),
            "low_region_log_primitive_derivative": sp.factor(
                sp.diff(sp.log(1 + 4 * y), y) - 1 / q
            ),
            "low_region_log_primitive_anchor": sp.log(1 + 4 * y).subs(y, 0),
            "low_region_exact_log_constant": low_log - 3624,
            "high_region_exact_rational_constant": high_total - sp.Rational(6639, 16),
            "finite_positive_Taylor_lower_for_exp_three": sum(
                sp.Rational(3) ** k / sp.factorial(k) for k in range(4)
            )
            - 13,
            **{
                "high_radial_power_integral_" + str(k): sp.integrate(
                    y ** (-k), (y, M, sp.oo)
                )
                - 1 / ((k - 1) * M ** (k - 1))
                for k in range(2, 6)
            },
        },
        "bounds": {
            "high_region_constant_below_415": high_total < 415,
            "constant_and_log_fit_one_common_3624_majorant": low_const + high_total
            < 3624,
            "two_crossed_channels_half_derivative_and_radial_factor_below_thirteen": sp.Rational(
                3624, 288
            )
            < 13,
            "actual_one_plus_four_M_below_ten_to_198": 1 + 4 * Mp < 10**198,
            "exp_three_positive_Taylor_lower_above_ten": exp_three_lower > 10,
            "log_of_ten_to_198_below_one_thousand": 3 * 198 < 1000,
            "leading_coefficient_rounded_up_to_14000": 13 * (1 + 1000) < 14000,
        },
    }
