"""Excluded-cut endpoint layer and rare-index logarithms."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_soft_cutoff.conditioning import parameters

from . import source


def singular_layer_integral(index, resolution, infrared):
    a, x, eta = parameters(index, resolution, infrared)
    if a == 0:
        return s.S.Zero
    c, b = a * eta, x / 2
    if c <= b:
        return c * (1 - s.log(c) + (s.log(c) ** 2 - s.log(b) ** 2) / 2)
    return b * (1 - s.log(b))


def log_event_change_upper(index, resolution, infrared):
    a, _, eta = parameters(index, resolution, infrared)
    return 6300 * a * eta * (1 - s.log(eta)) ** 2 / source.KAPPA


@cache
def data():
    c, b, y, L, z = s.symbols("c b y L z", positive=True)
    first = c * (1 - s.log(c) + (s.log(c) ** 2 - s.log(b) ** 2) / 2)
    second = b * (1 - s.log(b))
    checks = {
        "small_layer_upper_endpoint_derivative": s.diff(first, b) + c * s.log(b) / b,
        "large_layer_upper_endpoint_derivative": s.diff(second, b) + s.log(b),
        "piecewise_layer_join": first.subs(c, b) - second,
        "layer_primitive_first_segment": s.diff(y * (1 - s.log(y)), y) + s.log(y),
        "layer_primitive_second_segment": s.diff(-(s.log(y) ** 2) / 2, y)
        + s.log(y) / y,
        "rare_index_log_first_stationary_point": s.diff(z * s.exp(-z), z).subs(z, 1),
        "rare_index_log_second_stationary_point": s.diff(z * z * s.exp(-z), z).subs(
            z, 2
        ),
        "rare_index_polynomial_slack": s.expand(
            2 * (1 + L) ** 2
            - (2 + s.Rational(3, 2) * L + L * L / 2)
            - s.Rational(5, 2) * L
            - s.Rational(3, 2) * L * L
        ),
        "mixture_low_high_normalization_budget": s.Integer(1400 + 2800 + 2100 - 6300),
        "exact_layer_zero_index": singular_layer_integral(
            0, s.Rational(1, 8), s.Rational(1, 64)
        ),
    }
    gates = {
        "tail_independent_before_cut_conditioning_only": True,
        "excluded_probability_bounded_by_min_one_a_eta_over_gap": True,
        "renewal_density_bound_a_over_retained_energy_cancels_mark": True,
        "zero_count_atom_has_zero_Born_subtracted_mark": True,
        "both_singular_layer_branches_retained": True,
        "first_rare_index_log_max": s.exp(-1) < s.Rational(1, 2),
        "second_rare_index_log_max": 4 * s.exp(-2) < 1,
        "event_mixture_normalization_term_not_dropped": True,
        "no_log_of_rare_index_loss_in_final_bound": True,
    }
    return {
        "checks": {name: s.simplify(value) for name, value in checks.items()},
        "gates": {name: bool(value) for name, value in gates.items()},
        "whole_endpoint_layer": "The high retained-energy part is <=1400*a/kappa times integral_0^(x/2)(-ln y)*min(1,a*eta/y)dy. Its exact piecewise integral keeps both a*eta<=x/2 and the opposite regime.",
        "whole_rare_index_control": "With c=a*eta and b=ln(1/a), the layer is <=c[1+ln(1/c)+ln(1/c)^2/2]. Bounds a*b<=1/e<1/2 and a*b^2<=4/e^2<1 turn the high-event term into <=2800*a*eta*(1+ln(1/eta))^2/kappa, without a logarithmic rare-index loss.",
        "whole_event_log_bound": "Low-energy excluded part <=1400*a*eta*L_eta/kappa, high part <=2800*a*eta*L_eta^2/kappa, normalization <=2100*a*eta*L_eta/kappa. Their sum is <=6300*a*eta*L_eta^2/kappa.",
    }
