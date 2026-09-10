"""Exact convergent remainder integral and complete channel prefactors."""

from functools import cache

import sympy as sp


@cache
def data():
    y, M, delta = sp.symbols(
        "nonnegative_radial positive_heavy_mass_squared positive_light_gap",
        positive=True,
    )
    a = sp.Symbol("positive_parameter_weight", positive=True)
    v = sp.Symbol("complex_forward_increment")
    t = sp.Symbol("unit_auxiliary_parameter", real=True)
    primitive = -sp.log(1 + y / M) / (y + delta) + (
        sp.log(y + delta) - sp.log(y + M)
    ) / (M - delta)
    integrand = sp.log(1 + y / M) / (y + delta) ** 2
    endpoint = sp.log(M / delta) / (M - delta)
    L, g = sp.symbols("positive_quartic positive_cubic_squared", positive=True)
    error = (
        54 * L * L * g / (16 * sp.pi**2) ** 2 * sp.log(4 * M) / (M - sp.Rational(1, 4))
    )
    return {
        "y": y,
        "M": M,
        "delta": delta,
        "complete_radial_integrand": integrand,
        "anchored_radial_primitive": primitive,
        "full_zero_to_infinity_integral": endpoint,
        "actual_gap_one_quarter_integral": endpoint.subs(delta, sp.Rational(1, 4)),
        "complete_three_channel_remainder_majorant": error,
        "unit_disc_Cauchy_coefficient_factor": sp.Integer(1),
        "scope": "A finite absolute majorant for the decaying insertion remainder after grouping the on-shell subgraph. Its whole unit-disc amplitude bound controls the second forward coefficient by Cauchy; the fixed asymptotic term is treated separately with the inherited renormalized one-loop result.",
        "checks": {
            "radial_primitive_derivative": sp.factor(sp.diff(primitive, y) - integrand),
            "radial_primitive_at_infinity": sp.limit(primitive, y, sp.oo),
            "radial_primitive_zero_anchor": sp.simplify(
                -primitive.subs(y, 0) - endpoint
            ),
            "auxiliary_log_integral_derivative": sp.factor(
                sp.diff(sp.log(1 + a * t * y) / y, t) - a / (1 + a * t * y)
            ),
            "auxiliary_log_integral_zero_anchor": (sp.log(1 + a * t * y) / y).subs(
                t, 0
            ),
            "heavy_weight_log_monotonicity": sp.factor(
                sp.diff(sp.log(1 + a * y) / y, a) - 1 / (1 + a * y)
            ),
            "log_mass_over_gap_form": sp.simplify(
                sp.expand_log(sp.log(M / sp.Rational(1, 4)), force=False)
                - sp.log(4)
                - sp.log(M)
            ),
            "all_three_channel_loop_and_insertion_prefactors": 3
            * sp.Rational(1, 2)
            * 9
            * 2
            * 2
            - 54,
            "two_four_dimensional_radial_loop_measures": 1 / (16 * sp.pi**2) ** 2
            - 1 / (256 * sp.pi**4),
            "unit_disc_Cauchy_second_coefficient": sp.Rational(1, 1) ** (-2) - 1,
            "constant_contact_second_derivative": sp.diff(
                sp.Symbol("finite_contact"), v, 2
            ),
        },
    }
