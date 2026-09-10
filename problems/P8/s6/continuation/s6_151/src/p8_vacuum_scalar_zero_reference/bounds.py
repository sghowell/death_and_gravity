"""All-radius bounds for the added full-propagator scale term."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_finite_contact.calibration import rational
from p8_vacuum_two_loop_wineglass import integrals as previous


def enclosure(L, tree_b2, ell=1200, ellH=600, Q=144):
    L, tree_b2, ell, ellH, Q = map(rational, (L, tree_b2, ell, ellH, Q))
    if (
        L <= 0
        or tree_b2 <= 0
        or not 0 <= ell <= 1200
        or not 0 <= ellH <= 600
        or not 0 < Q <= 144
    ):
        raise ValueError("Need L,tree>0; 0<=ell<=1200; 0<=ellH<=600; 0<Q_lower<=144")
    old = (240 + 40 * ellH * (ellH + 4)) * L**3 / Q**2
    extra = ell * (12 + 40 * ellH) * L**3 / Q**2
    F = (ell**2 / 2 + ell + 22 + s.Rational(4, 3)) / Q**2
    ref = tree_b2 * L**2 * F
    return {
        "old_fully_subtracted_wineglass_upper": old,
        "added_complete_scale_term_upper": extra,
        "finite_proper_MS_reference_absolute_upper": F,
        "local_reference_b2_absolute_upper": ref,
        "local_reference_b2_relative_upper": L**2 * F,
        "complete_wineglass_MS_upper": old + extra + ref,
        "complete_wineglass_MS_relative_upper": (old + extra + ref) / tree_b2,
    }


@cache
def data():
    L, g, M, Q, ell, ellH = s.symbols("L g M Q ell ellH", positive=True)
    z, x = s.symbols("z x", real=True)
    previous_triangle = previous.data()["same_full_triangle_radial_integral"]
    a = x * (1 - x)
    # The inherited routing has |z|<=3, Delta>=1/4 and g/M<L/3.
    core = 3 * (2 * L) * L**2 * (2 / Q) * ell / Q
    heavy = 3 * (20 * L**2 * g / Q) * (2 * ellH / M) * ell / Q
    return {
        "outer_IR": -s.Integral(s.log(1 - a * z), (x, 0, 1)) / Q,
        "outer_IR_modulus_upper": 2 / Q,
        "same_full_triangle_integral": previous_triangle,
        "heavy_numerator_upper": 20 * L**2 * g / (s.Symbol("y", nonnegative=True) + M),
        "all_channel_scale_majorant": ell * (12 + 40 * ellH) * L**3 / Q**2,
        "checks": {
            "outer_parameter_weight_max": a.subs(x, s.Rational(1, 2))
            - s.Rational(1, 4),
            "outer_log_series_radius": 3 * s.Rational(1, 4) - s.Rational(3, 4),
            "core_all_channel_prefactor": s.expand(core - 12 * ell * L**3 / Q**2),
            "heavy_all_channel_prefactor": s.expand(
                heavy.subs(g, L * M / 3) - 40 * ell * ellH * L**3 / Q**2
            ),
            "scale_sum_prefactor": s.expand(
                core + heavy.subs(g, L * M / 3) - ell * (12 + 40 * ellH) * L**3 / Q**2
            ),
            "pi_squared_finite_anchor_majorant": s.Rational(4**2, 12)
            - s.Rational(4, 3),
        },
        "scope": "The same all-radius full-heavy routing bounds and complex radius-one forward disc as S6.124. No propagator expansion or loop-momentum cutoff is used. The local reference b2 is evaluated by its exact tree parameter dependence, not by a whole-amplitude bound on the large constant quartic term.",
    }
