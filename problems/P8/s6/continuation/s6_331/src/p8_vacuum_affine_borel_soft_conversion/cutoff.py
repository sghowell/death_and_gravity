"""Quantitative conditioned infrared cutoff for the original a and Delta marks."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_soft_cutoff.conditioning import parameters

from . import source


def index_mean_cutoff_upper(index, resolution, infrared):
    a, _, eta = parameters(index, resolution, infrared)
    return 4500 * a * eta / source.KAPPA


def conversion_mean_cutoff_upper(index, resolution, infrared):
    a, _, eta = parameters(index, resolution, infrared)
    return 42000 * a * eta * (1 - s.log(eta)) / source.KAPPA


def original_index_mean_cutoff_upper(resolution, infrared):
    _, _, eta = parameters(0, resolution, infrared)
    return 4000 * eta / source.KAPPA**2


def original_conversion_mean_cutoff_upper(resolution, infrared):
    _, _, eta = parameters(0, resolution, infrared)
    return 34000 * eta * (1 - s.log(eta)) / source.KAPPA**2


@cache
def data():
    t, u, eta, z = s.symbols("t u eta z", positive=True)
    f = lambda v: v * (1 - s.log(v))
    checks = {
        "positive_tail_entropy_subadditivity_identity": s.expand(
            f(t)
            + f(u)
            - f(t + u)
            - ((t + u) * s.log(t + u) - t * s.log(t) - u * s.log(u))
        ),
        "missing_energy_log_integral": s.integrate(1 - s.log(z), (z, 0, eta))
        - eta * (2 - s.log(eta)),
        "conditional_index_same_event_and_mixture_budget": s.Integer(
            1700 + 2 * 1400 - 4500
        ),
        "conditional_conversion_same_event_and_mixture_constant": s.Integer(
            2 * 11000 + 20000 - 42000
        ),
        "conditional_conversion_log_coefficient_margin": s.Integer(
            42000 - 11000 - 20000 - 11000
        ),
    }
    gates = {
        "original_index_cutoff_roundup": s.Rational(4, 5) * 4500 < 4000,
        "original_conversion_cutoff_roundup": s.Rational(4, 5) * 42000 < 34000,
        "one_distinguished_emission_not_conditioned_Poisson_independence": True,
        "subadditivity_controls_tail_log_without_log_of_rare_index": True,
        "same_Born_mark_subtracted_before_event_mixture": True,
        "relative_conditioning_loss_not_loose_absolute_error": True,
        "uniform_ordered_eta_and_x": True,
        "no_cutoff_rate_for_remaining_energy_log_is_claimed": True,
    }
    return {
        "checks": {name: s.simplify(value) for name, value in checks.items()},
        "gates": {name: bool(value) for name, value in gates.items()},
        "whole_conditional_tail_log_bound": "OnR<=x, f(T_eta)<=sum_(w<eta)f(w), f(t)=t(1+ln1/t). The positive distinguished-emission identity gives E[f(T_eta)|cut]<=a*eta*(2+ln1/eta), without an extra ln1/a loss.",
        "whole_conditioned_cutoff_moduli": "At0<eta<=x<=1/8, actual finite-cutoff versus limiting conditioned a means differ by<=4500a*eta/kappa and Delta means by<=42000a*eta*(1+ln1/eta)/kappa. Original physical bounds are4000eta/kappa^2 and34000eta*(1+ln1/eta)/kappa^2.",
        "whole_cutoff_boundary": "This compares a and Delta marks. It does not claim a quantitative cutoff rate after also changing the remaining-energy logarithm inZ, nor any new radiation dynamics or hard matching.",
    }
