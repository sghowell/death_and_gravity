"""First-sheet routings, anchored logarithm bounds and independent radial integrals."""

from functools import cache

import sympy as sp
from p8_vacuum_two_loop_insertions import kernel as inherited


@cache
def data():
    s = sp.Symbol("complex_Minkowski_invariant")
    re, im = sp.symbols("real_s imaginary_s", real=True)
    x, a = sp.symbols("unit_heavy_parameter light_bubble_weight", real=True)
    y, t = sp.symbols("nonnegative_radial nonnegative_sqrt_radial", nonnegative=True)
    M = sp.Symbol("positive_heavy_mass_squared", positive=True)
    p = sp.Matrix([sp.I * (s + 3) / (2 * sp.sqrt(3)), (s - 3) / (2 * sp.sqrt(3)), 0, 0])
    pc = p.subs(s, re + sp.I * im)
    norm = sp.expand((sp.conjugate(pc).T * pc)[0])
    delta = 1 + x * (M - 1) - x * (1 - x) * s
    delta0 = 1 + x * (M - 1)
    D = (16 * sp.sqrt(y) + 12) / (y + 4)
    P_integer = -1 / (y + 1) + 1 / (2 * (y + 1) ** 2)
    P_half = 3 * sp.atan(t) / 4 - 5 * t / (4 * (1 + t**2)) + t / (2 * (1 + t**2) ** 2)
    inherited_data = inherited.data()
    z = inherited_data["z"]
    b = inherited_data["b"]
    pref = inherited_data["prefactor"]
    R_real = -pref * sp.Integral(sp.log(1 + b * z) / z, (inherited_data["x"], 0, 1))
    v, radius = sp.symbols(
        "nonnegative_disc_distance positive_outer_radius", positive=True
    )
    majorant = v**2 / (radius * (radius - v))
    checks = {
        "literal_entire_lightlike_coordinate_mass_shell": sp.expand((p.T * p)[0] + s),
        "exact_Hermitian_routing_norm": sp.expand(norm - (re**2 + im**2 + 9) / 6),
        "routing_norm_bound_from_abs_s_at_most_three": sp.Rational(9 + 9, 6) - 3,
        "combined_light_heavy_denominator": sp.expand(delta - delta0 + x * (1 - x) * s),
        "combined_denominator_lower_gap": sp.expand(
            delta.subs(s, 3) - (1 + x * (M - 4)) - 3 * x**2
        ),
        "heavy_parameter_product_upper": sp.expand(
            sp.Rational(1, 4) - x * (1 - x) - (x - sp.Rational(1, 2)) ** 2
        ),
        "shifted_log_positive_quarter_gap": sp.expand(
            1
            + a * (y / 2 - 3)
            - (1 + a * y) / 4
            - (3 * (sp.Rational(1, 4) - a) + a * y / 4)
        ),
        "log_shift_norm_rounding_positive_margin": 4**2 - (2 * sp.sqrt(3)) ** 2 - 4,
        "parameter_fraction_monotone": sp.factor(
            sp.diff(a / (1 + a * y), a) - 1 / (1 + a * y) ** 2
        ),
        "log_difference_uniform_seven_SOS": sp.factor(
            (7 - D.subs(y, t**2)) * (t**2 + 4)
            - 7 * (t - sp.Rational(8, 7)) ** 2
            - sp.Rational(48, 7)
        ),
        "log_below_square_root_derivative": sp.factor(
            sp.diff(t - sp.log(1 + t**2), t) - (t - 1) ** 2 / (1 + t**2)
        ),
        "log_below_square_root_anchor": (t - sp.log(1 + t**2)).subs(t, 0),
        "shifted_log_argument_product_majorant": 2 * (1 + y) - (y + 2) - y,
        "heavy_radial_denominator_margin": sp.expand(
            y / 2 + M - 3 - (y + M) / 2 - (M / 2 - 3)
        ),
        "integer_radial_primitive": sp.factor(sp.diff(P_integer, y) - y / (y + 1) ** 3),
        "integer_radial_infinity_anchor": sp.limit(P_integer, y, sp.oo),
        "integer_radial_zero_anchor": P_integer.subs(y, 0) + sp.Rational(1, 2),
        "half_power_radial_primitive": sp.factor(
            sp.diff(P_half, t) - 2 * t**4 / (1 + t**2) ** 3
        ),
        "half_power_radial_zero_anchor": P_half.subs(t, 0),
        "half_power_radial_infinity_anchor": sp.limit(P_half, t, sp.oo) - 3 * sp.pi / 8,
        "mixed_shift_log_integral_majorant": 16 * (3 * sp.pi / 8)
        + 12 * sp.Rational(1, 2)
        - (6 * sp.pi + 6),
        "mixed_reference_denominator_integral_majorant": sp.Rational(3, 2)
        * (3 * sp.pi / 8)
        - 9 * sp.pi / 16,
        "mixed_total_rational_slack_using_pi_below_four": 33
        - (6 * 4 + 6 + sp.Rational(9, 16) * 4)
        - sp.Rational(3, 4),
        "same_heavy_integral_majorant": sp.limit(
            P_half + 7 * P_integer.subs(y, t**2), t, sp.oo
        )
        - (P_half + 7 * P_integer.subs(y, t**2)).subs(t, 0)
        - (3 * sp.pi / 8 + sp.Rational(7, 2)),
        "same_heavy_rational_endpoint": sp.Rational(3, 8) * 4 + sp.Rational(7, 2) - 5,
        "nested_radial_rational_endpoint": sp.Rational(1, 2)
        + sp.Rational(3, 8) * 4
        - 2,
        "heavy_squared_Feynman_weight": sp.integrate(2 * x, (x, 0, 1)) - 1,
        "same_actual_inherited_real_axis_remainder": sp.simplify(
            R_real - inherited_data["decaying_remainder"]
        ),
        "OS_Cauchy_geometric_tail": sp.factor(
            majorant - ((1 - v / radius) ** -1 - 1 - v / radius)
        ),
        "OS_Cauchy_unit_inner_radius": (majorant / v**2).subs({radius: 2, v: 1})
        - sp.Rational(1, 2),
        "OS_line_Taylor_integral": sp.integrate(1 - x, (x, 0, 1)) - sp.Rational(1, 2),
    }
    return {
        "s": s,
        "M": M,
        "x": x,
        "y": y,
        "entire_external_routing": p,
        "Hermitian_norm_squared": norm,
        "combined_light_heavy_Delta": delta,
        "constant_reference_Delta": delta0,
        "anchored_log_difference_majorant": D,
        "integer_radial_primitive": P_integer,
        "half_power_radial_primitive_in_sqrt_y": P_half,
        "integer_radial_integral": sp.Rational(1, 2),
        "half_power_radial_integral": 3 * sp.pi / 8,
        "mixed_reference_subtracted_integral_upper_without_loop_measure": 33,
        "same_heavy_integral_upper_without_loop_measure": 5,
        "nested_R_integral_upper_without_couplings_and_loop_measure": 4,
        "inherited_real_axis_R": R_real,
        "real_axis_R_absolute_upper": pref
        * sp.log(1 + z / inherited_data["heavy_mass_squared"])
        / z,
        "OS_Cauchy_tail_factor": majorant,
        "scope": "For |s-1|<=2, all displayed integrated representations are continued from the massive Euclidean domain with a strict first-sheet denominator gap. The mixed integral is first subtracted at zero external momentum. The estimates use pi<4 and log(1+y)<=sqrt(y), not a momentum cutoff or a large-heavy-mass expansion.",
        "checks": checks,
    }
