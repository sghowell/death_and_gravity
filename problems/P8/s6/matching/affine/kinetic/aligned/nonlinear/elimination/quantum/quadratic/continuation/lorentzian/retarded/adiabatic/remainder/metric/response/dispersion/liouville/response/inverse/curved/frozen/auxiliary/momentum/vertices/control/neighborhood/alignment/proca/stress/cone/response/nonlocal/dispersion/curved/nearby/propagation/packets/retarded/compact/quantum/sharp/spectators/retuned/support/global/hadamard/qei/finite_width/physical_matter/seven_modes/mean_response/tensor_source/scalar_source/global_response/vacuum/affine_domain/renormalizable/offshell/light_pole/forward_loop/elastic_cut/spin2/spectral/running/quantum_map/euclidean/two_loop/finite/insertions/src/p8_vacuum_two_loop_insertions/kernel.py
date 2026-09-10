"""Actual one-loop multiplier split with a decaying complex-strip remainder."""

from functools import cache

import sympy as sp
from p8_vacuum_light_pole import kernel as parent


@cache
def data():
    d = parent.data()
    x = d["Feynman_parameter"]
    b = d["normalized_parameter_weight"]
    pref = d["self_energy_prefactor"]
    z = sp.Symbol("complex_one_minus_line_invariant")
    t = sp.Symbol("unit_auxiliary_parameter", real=True)
    re, im = sp.symbols("real_increment imaginary_increment", real=True)
    y = sp.Symbol("nonnegative_centered_radial", nonnegative=True)
    M = d["heavy_mass_squared"]
    general_mass = sp.Symbol("positive_general_mass", positive=True)
    bs = sp.Symbol("positive_parameter_weight", positive=True)
    q = b - sp.log(1 + b * z) / z
    remainder = -sp.log(1 + b * z) / z
    primitive = -sp.log(1 + b * t * z) / z
    return {
        "x": x,
        "z": z,
        "t": t,
        "b": b,
        "prefactor": pref,
        "heavy_mass_squared": M,
        "relative_kernel": pref * sp.Integral(q, (x, 0, 1)),
        "fixed_asymptotic_multiplier": pref * sp.Integral(b, (x, 0, 1)),
        "decaying_remainder": pref * sp.Integral(remainder, (x, 0, 1)),
        "decaying_remainder_auxiliary_integrand": -b / (1 + b * t * z),
        "radial_remainder_majorant": 2 * pref * sp.log(1 + y / M) / y,
        "radial_majorant_continuous_at_zero": 2 * pref / M,
        "scope": "The actual fixed on-shell one-loop kernel, analytically continued in the stated complex strip. The constant/remainder split must be recombined with the inherited outer subtraction; it is not an unregulated finite bubble.",
        "checks": {
            "same_parent_OS_kernel_divided_by_inverse": sp.simplify(
                d["local_analytic_on_shell_kernel"].subs(
                    d["Minkowski_invariant"], 1 - z
                )
                / z
                - q
            ),
            "constant_plus_decaying_split": sp.simplify(q - b - remainder),
            "auxiliary_remainder_primitive": sp.simplify(
                sp.diff(primitive, t) + b / (1 + b * t * z)
            ),
            "auxiliary_remainder_anchor": primitive.subs(t, 0),
            "auxiliary_remainder_endpoint": sp.simplify(
                primitive.subs(t, 1) - remainder
            ),
            "continuous_kernel_at_light_shell": sp.limit(q, z, 0),
            "continuous_remainder_at_light_shell": sp.limit(
                -sp.log(1 + bs * z) / z, z, 0
            )
            + bs,
            "complex_fraction_modulus_gap": sp.expand(
                (1 + re) ** 2 + im**2 - (re**2 + im**2) - (1 + 2 * re)
            ),
            "radial_majorant_at_zero": sp.limit(sp.log(1 + y / general_mass) / y, y, 0)
            - 1 / general_mass,
            "strip_denominator_radial_gap": sp.expand(
                1
                + bs * t * (y / 2 - 2)
                - (1 + bs * t * y) / 2
                - (sp.Rational(1, 2) - 2 * bs * t)
            ),
        },
    }
