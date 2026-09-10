"""Literal stress triangles, on-shell Ward normalization and finite spin-two slope."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_light_pole import kernel


@cache
def data():
    x, z, t = sp.symbols(
        "heavy_parameter unit_pair_parameter momentum_transfer", real=True
    )
    M, g = sp.symbols("heavy_mass_squared cubic_squared", positive=True)
    d1 = x * M + (1 - x) ** 2
    b = x * (1 - x) / d1
    wL = x * x * (1 - x)
    wH = x * (1 - x) ** 2
    hL = (1 - x) ** 2 * z * (1 - z)
    hH = x * x * z * (1 - z)
    termL = wL / (d1 - hL * t)
    termH = wH / (d1 - hH * t)
    pref = g / (16 * sp.pi**2)
    f = pref * (termL + termH - (wL + wH) / d1)
    slope = sp.integrate(sp.diff(f, t).subs(t, 0), (z, 0, 1))
    chi, zz, mA, mB = sp.symbols(
        "single_mass_parameter pair_shift_parameter equal_line_mass_squared single_line_mass_squared",
        real=True,
    )
    pp, qq, pq = sp.symbols(
        "Euclidean_p_squared Euclidean_q_squared Euclidean_p_dot_q", real=True
    )
    qv = sp.Matrix(sp.symbols("q0:4", real=True))
    pv = sp.Matrix(sp.symbols("p0:4", real=True))
    kv = sp.Matrix(sp.symbols("k0:4", real=True))
    dot = lambda a, b: (a.T * b)[0]
    general = (
        (1 - chi) * mA
        + chi * mB
        + zz * (1 - zz) * qq
        + chi * (1 - chi) * pp
        + 2 * zz * chi * pq
    )
    shift = kv + zz * qv - chi * pv
    weighted = (
        (1 - chi - zz) * (dot(kv, kv) + mA)
        + zz * (dot(kv + qv, kv + qv) + mA)
        + chi * (dot(kv - pv, kv - pv) + mB)
    )
    on_shell = (1 - chi) * mA + chi * mB - chi * (1 - chi) - zz * (1 - chi - zz) * t
    a0 = kernel.data()
    parentM = model.data()["heavy_mass_squared"]
    parentG = model.data()["cubic_coupling_squared"]
    actual_map = {M: parentM, g: parentG}
    old_weight = a0["normalized_parameter_weight"].subs(a0["Feynman_parameter"], x)
    old_second = (
        sp.diff(a0["self_energy_momentum_kernel"], a0["Minkowski_invariant"], 2)
        .subs(a0["Minkowski_invariant"], 1)
        .subs(a0["Feynman_parameter"], x)
    )
    pa = model.data()["actual_parameters"]
    Mp = pa["heavy_mass_squared"]
    gp = pa["cubic_coupling_squared"]
    slope_upper = gp / (2592 * Mp * Mp)
    y, D, K = sp.symbols(
        "positive_radial_y positive_Delta positive_transfer_weight", positive=True
    )
    radial = sp.integrate(y / (y + D) ** 3, (y, 0, sp.oo))
    equal_mass_integral = sp.integrate((x * (1 - x) / (1 - x + x * x)) ** 2, (x, 0, 1))
    return {
        "x": x,
        "z": z,
        "t": t,
        "M": M,
        "g": g,
        "on_shell_two_point_denominator": d1,
        "light_stress_triangle": termL,
        "heavy_stress_triangle": termH,
        "once_Ward_subtracted_unit_square_integrand": f,
        "spin_two_slope_parameter_integrand": pref * b * b / 6,
        "actual_spin_two_slope_rational_upper": slope_upper,
        "actual_heavy_mass_squared": Mp,
        "actual_cubic_squared": gp,
        "equal_mass_two_distinct_field_integral": equal_mass_integral,
        "analytic_transfer_disc": "|t|<4, conservatively; the light-pair threshold is distinct from the mixed two-point threshold",
        "scope": "Actual one-loop spin-two matter stress form factor in the minimally curved polynomial model. Not the full gravitational amplitude, its all-channel coefficient, pure graviton loops, Regge remainder or bounce matching.",
        "checks": {
            "literal_four_dimensional_three_propagator_shift": sp.expand(
                weighted
                - dot(shift, shift)
                - general.subs({pp: dot(pv, pv), qq: dot(qv, qv), pq: dot(pv, qv)})
            ),
            "on_shell_triangle_parameter_polynomial": sp.expand(
                general.subs({pp: -1, qq: -t, pq: t / 2}) - on_shell
            ),
            "stress_on_light_parameter_denominator": sp.expand(
                on_shell.subs({chi: x, mA: 1, mB: M, zz: (1 - x) * z}) - (d1 - hL * t)
            ),
            "stress_on_heavy_parameter_denominator": sp.expand(
                on_shell.subs({chi: 1 - x, mA: M, mB: 1, zz: x * z}) - (d1 - hH * t)
            ),
            "actual_parent_two_point_denominator_match": sp.factor(
                d1.subs(M, parentM)
                - a0["on_shell_Delta"].subs(a0["Feynman_parameter"], x)
            ),
            "sum_of_two_stress_triangle_Ward_weights": sp.expand(wL + wH - x * (1 - x)),
            "once_subtracted_stress_normalization_at_zero": sp.factor(f.subs(t, 0)),
            "actual_once_fixed_light_residue_Ward_normalization": sp.factor(
                (pref * (wL + wH) / d1).subs(actual_map)
                - a0["self_energy_prefactor"] * old_weight
            ),
            "actual_two_point_second_derivative_spin_two_slope_relation": sp.factor(
                slope.subs(actual_map) - a0["self_energy_prefactor"] * old_second / 6
            ),
            "integrated_slope_is_one_sixth_two_point_second_derivative": sp.factor(
                slope - pref * b * b / 6
            ),
            "light_transfer_ratio_positive_gap": sp.factor(
                z * (1 - z) - hL / d1 - z * (1 - z) * x * M / d1
            ),
            "heavy_transfer_ratio_positive_gap": sp.factor(
                x * z * (1 - z) / M
                - hH / d1
                - x * z * (1 - z) * (1 - x) ** 2 / (M * d1)
            ),
            "quarter_bound_parameter_gap": sp.expand(
                sp.Rational(1, 4) - z * (1 - z) - (z - sp.Rational(1, 2)) ** 2
            ),
            "exact_subtracted_resolvent_quadratic_remainder": sp.factor(
                1 / (D - K * t)
                - 1 / D
                - K * t / D**2
                - K * K * t * t / (D * D * (D - K * t))
            ),
            "finite_four_dimensional_radial_triangle_integral": radial - 1 / (2 * D),
            "parameter_weight_square_upper_integral": sp.integrate(
                (1 - x) ** 2, (x, 0, 1)
            )
            - sp.Rational(1, 3),
            "equal_mass_independent_closed_integral": sp.simplify(
                equal_mass_integral - (45 - 8 * sp.pi * sp.sqrt(3)) / 27
            ),
        },
        "bounds": {
            "actual_M_above_one_for_common_transfer_disc": Mp > 1,
            "actual_positive_spin_two_slope_majorant": slope_upper > 0,
            "actual_spin_two_slope_below_one_e_minus_405": slope_upper
            < sp.Rational(1, 10**405),
        },
    }
