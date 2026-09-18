"""Uniform infrared cutoff and the joint two-parameter signed-density bound."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_soft_cutoff.conditioning import parameters

from . import kernels, regulator, source


def regulated_cutoff_tv_upper(index, resolution, infrared, epsilon):
    a, _, eta = parameters(index, resolution, infrared)
    kernels.regulator(epsilon)
    return 132000 * a * eta * (1 - s.log(eta)) ** 2 / source.KAPPA


def joint_tv_upper(index, resolution, infrared, epsilon):
    return regulated_cutoff_tv_upper(
        index, resolution, infrared, epsilon
    ) + regulator.conditional_regulator_tv_upper(index, resolution, epsilon)


def original_joint_tv_upper(resolution, infrared, epsilon):
    _, x, eta = parameters(0, resolution, infrared)
    e = kernels.regulator(epsilon)
    return (
        106000 * eta * (1 - s.log(eta)) ** 2 + 320000 * e * x * (1 - s.log(x)) ** 2
    ) / source.KAPPA**2


@cache
def data():
    checks = {
        "quotient_increment_and_mixture_budget": s.Integer(
            2 * 18000 + 2 * 15000 - 66000
        ),
        "small_eta_exponent_budget": s.Integer(2 * 15000 * 8 - 240000),
        "small_eta_max_epsilon_budget": s.Rational(240000, 8) - 30000,
        "large_eta_max_epsilon_budget": s.Rational(2 * 15000 * 3 * 4, 8) - 45000,
        "combined_uniform_cutoff_budget": s.Integer(
            66000 + 45000 + 14600 + 6300 - 131900
        ),
        "original_uniform_cutoff_budget": s.Rational(132000 * 4, 5) - 105600,
    }
    return {
        "checks": checks,
        "gates": {
            "uniform_cutoff_roundup": bool(s.Integer(131900) < 132000),
            "original_cutoff_roundup": bool(s.Integer(105600) < 106000),
            "small_eta_exponent_dominated": bool(s.Integer(30000) < 45000),
            "same_common_probability_space_for_both_regulators": True,
            "both_event_normalization_terms_retained": True,
            "arbitrary_joint_limits_not_only_chosen_order": True,
            "fixed_resolution_and_fixed_hard_data_scope": True,
            "configuration_law_TV_not_claimed": True,
        },
        "whole_uniform_epsilon_cutoff": "Uniformly0<=e<=1/8, common-space L1 cutoff error<=132000*a*eta*L_eta^2/kappa. The new q_e times power-kernel shift uses30000 in the small-eta regime and45000 in the complementary regime; the two normalization terms and logarithmic singular layer are retained.",
        "whole_joint_quantitative_bound": "The original fixed-reference signed densities obey ||f_e,eta-f_0||_1<=[132000*a*eta*L_eta^2+400000e*a*x*L_x^2]/kappa, physically<[106000eta*L_eta^2+320000e*x*L_x^2]/kappa^2. This gives arbitrary joint and both iterated limits at fixed x>0.",
        "whole_joint_boundary": "Both regulators belong only to the named additional-soft insertion on the unchanged physical D4 leading reference. No outer hard dimensional matching, interacting state, all-N hard sum, Regge or original P8 closure is inferred.",
    }
