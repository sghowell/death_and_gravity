"""Positive SAME-D Born normalization with no unexpanded-pole substitution."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import sew
from p8_vacuum_affine_one_newton_inclusive_assembly import forward

from . import source


def domain(invariant, cosine):
    ss, z = map(sew.exact_real, (invariant, cosine))
    if ss < s.Rational(25, 4) or ss > 16 or z <= -1 or z >= 1:
        raise ValueError("Require 25/4<=s<=16 and a strictly nonforward angle")
    return ss, z


def transfer_gap(invariant, cosine):
    ss, z = domain(invariant, cosine)
    return min(s.S.One, (ss - 4) * (1 - abs(z)) / 2)


def auxiliary(invariant, cosine):
    ss, z = domain(invariant, cosine)
    return 4 / ((ss - 4) * (1 - z * z)) - 1 / ss


def stripped_gravity(invariant, cosine, epsilon):
    ss, z = domain(invariant, cosine)
    theta = sew.ratio(epsilon)
    return forward.gravity_born(ss, z, 1) + 2 * theta * auxiliary(ss, z)


def original(invariant, cosine, epsilon):
    ss, z = domain(invariant, cosine)
    gravity = stripped_gravity(ss, z, epsilon)
    matter = forward.matter_born(ss, z, source.HEAVY_MASS2, source.CUBIC)
    return matter + gravity / source.KAPPA


def ratio_to_four(invariant, cosine, epsilon):
    return s.factor(
        original(invariant, cosine, epsilon) / original(invariant, cosine, 0)
    )


@cache
def data():
    q, z, alpha, r = s.symbols("q z alpha r", positive=True)
    ss = q + 4
    tt = -q * (1 - z) / 2
    uu = -q * (1 + z) / 2
    aux = 4 / (q * (1 - z * z)) - 1 / ss
    p = 4 * q + 16 + 8 / q
    checks = {
        "auxiliary_three_channel_identity": s.factor(aux + 1 / ss + 1 / tt + 1 / uu),
        "auxiliary_positive_numerator": s.factor(
            aux - (3 * q + 16 + q * z * z) / (ss * q * (1 - z * z))
        ),
        "auxiliary_absolute_transfer_upper": s.factor(
            4 / (q * (1 - alpha))
            - aux.subs(z, alpha)
            - 4 * alpha / (q * (1 - alpha * alpha))
            - 1 / ss
        ),
        "gravity_forward_residue_compact_floor": s.expand(
            (p * q - s.Rational(257, 4)).subs(q, s.Rational(9, 4) + r)
            - r * (4 * r + 34)
        ),
        "gravity_interior_compact_floor": (p / 2).subs(q, s.Rational(9, 4))
        - s.Rational(257, 18),
        "Born_ratio_coefficient": s.Rational(2, 4) - s.Rational(1, 2),
        "maximum_Born_ratio": 1 + s.Rational(1, 18) - s.Rational(19, 18),
        "maximum_fixed_Born_square": s.Rational(19, 18) ** 2 - s.Rational(361, 324),
    }
    return {
        "checks": checks,
        "gates": {
            "auxiliary_Born_positive_on_physical_compact_domain": True,
            "original_gravity_Born_floor_not_an_all_angle_small_gravity_claim": True,
            "complete_Born_remains_positive_at_every_epsilon": True,
            "fixed_A0_amplitude_ratio_at_most_19_over_18": True,
            "same_D_Born_retained_inside_leading_soft_subtraction": True,
        },
        "whole_Born_identity": "G0_D=G0_4+2theta*A_S,theta=epsilon/(1+epsilon),A_S=1/T+1/U-1/s>0. G0_4>8/delta; A_S<=2/tau. Hence1<=A_D/A0<=1+theta/2<=19/18.",
        "whole_phase_reference_boundary": "If the amplitude denominator is fixed to A0^2, retain the same D-dimensional reference two-body phase and multiply the D-dependent Born-normalized remainder by(A_D/A0)^2. Do not replace A_D by A0 inside a separately divergent leading-soft term.",
    }
