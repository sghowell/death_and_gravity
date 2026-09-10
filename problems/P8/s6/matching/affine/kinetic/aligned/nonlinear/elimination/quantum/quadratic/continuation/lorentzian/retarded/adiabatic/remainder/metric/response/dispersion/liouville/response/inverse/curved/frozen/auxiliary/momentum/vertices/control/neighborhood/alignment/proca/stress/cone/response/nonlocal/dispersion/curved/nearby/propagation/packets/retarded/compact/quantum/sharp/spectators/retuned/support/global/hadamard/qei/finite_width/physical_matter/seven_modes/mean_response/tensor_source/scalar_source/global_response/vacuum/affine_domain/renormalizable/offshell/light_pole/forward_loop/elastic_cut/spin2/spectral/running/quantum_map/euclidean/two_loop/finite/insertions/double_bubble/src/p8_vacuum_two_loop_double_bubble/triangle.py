"""Unexpanded finite heavy triangles and the convergent all-radius majorant."""

from functools import cache

import sympy as sp
from p8_vacuum_two_loop_insertions import routing


@cache
def data():
    y, M, delta = sp.symbols(
        "nonnegative_radial positive_heavy_mass_squared positive_light_gap",
        positive=True,
    )
    integrand = y / ((y + delta) ** 2 * (y + M))
    primitive = M * (sp.log(y + delta) - sp.log(y + M)) / (M - delta) ** 2 + delta / (
        (M - delta) * (y + delta)
    )
    exact = M * sp.log(M / delta) / (M - delta) ** 2 - 1 / (M - delta)
    g = sp.Symbol("positive_cubic_squared", positive=True)
    gap = sp.Rational(1, 4)
    checks = {
        "radial_primitive_derivative": sp.factor(sp.diff(primitive, y) - integrand),
        "radial_primitive_at_infinity": sp.limit(primitive, y, sp.oo),
        "radial_primitive_zero_anchor": sp.simplify(-primitive.subs(y, 0) - exact),
        "same_light_gap": routing.data()["centered_light_real_mass_gap"] - gap,
        "heavy_radial_denominator_margin": sp.expand(
            y / 2 + M - 9 - (y + M) / 2 - (M / 2 - 9)
        ),
        "two_unexpanded_heavy_channels_times_denominator_half": 2 * 2 - 4,
        "conservative_mass_squared_ratio_margin_at_thirty_two": 2 * (32 - gap) ** 2
        - 32**2
        - sp.Rational(7937, 8),
        "mass_squared_ratio_margin_increases_above_thirty_two": sp.expand(
            2 * (M - gap) ** 2 - M**2 - (M**2 - M + sp.Rational(1, 8))
        ),
        "exact_triangle_bound_remainder_identity": sp.expand(
            exact - M * sp.log(M / delta) / (M - delta) ** 2 + 1 / (M - delta)
        ),
    }
    return {
        "y": y,
        "M": M,
        "delta": delta,
        "all_radius_majorant_integrand": integrand,
        "anchored_radial_primitive": primitive,
        "exact_zero_to_infinity_radial_integral": exact,
        "gap_one_quarter_radial_integral": exact.subs(delta, gap),
        "each_summed_two_channel_heavy_triangle_upper": 4
        * g
        * exact.subs(delta, gap)
        / (16 * sp.pi**2),
        "simplified_each_triangle_modulus_upper": 8
        * g
        * sp.log(4 * M)
        / (16 * sp.pi**2 * M),
        "same_channel_routings": routing.data()["channels"],
        "scope": "Each outer full vertex contributes two genuinely finite heavy triangles. The same complex routing has |w|<3, so Re((q+w)^2+M)>=(q^2+M)/2 at M>32. The loop integral runs to infinity; no heavy derivative expansion or hard cutoff is used.",
        "checks": checks,
    }
