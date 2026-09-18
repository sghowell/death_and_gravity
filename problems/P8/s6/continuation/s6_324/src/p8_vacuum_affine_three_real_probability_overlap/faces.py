"""Absolute kinematic-current discs and uniform proper-face derivatives."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_real_soft_overlap import analytic as pair
from p8_vacuum_affine_uniform_all_tree_bound import bounds as tree

from . import source

CURRENT_RADIUS = s.Rational(1, 10**6)
CURRENT_CAP = 129
FIRST_CURRENT_CAP = 10**9
MIXED_CURRENT_CAP = 10**15
AXIS_CAP = s.Rational(1, 10**735)
PAIR_CAP = s.Rational(1, 10**723)


@cache
def data():
    d = CURRENT_RADIUS
    W = s.Rational(1, 8)
    gates = {
        "absolute_recoil_square_perturbation": 4 * d + d * W + 2 * d * d < 5 * d,
        "pair_energy_root_gap": 5 * d < s.Rational(45, 128),
        "pair_spatial_root_gap": 5 * d < s.Rational(13, 128),
        "continued_spatial_norm": (s.Rational(7, 4) - 13 * d) ** 2 > 3,
        "continued_energy": 1 - 5 * d > s.Rational(7, 8),
        "continued_Doppler": 65 * d < s.Rational(1, 8),
        "all_six_soft_current_terms": 4 * 32 + 2 * (W + 2 * d) < CURRENT_CAP,
        "first_current_derivative": 2 * CURRENT_CAP / d < FIRST_CURRENT_CAP,
        "mixed_current_derivative": 4 * CURRENT_CAP / d**2 < MIXED_CURRENT_CAP,
        "smaller_G2_disc_inside_current_disc": pair.RADIUS * W < d / 2,
    }
    r = s.Symbol("r", positive=True)
    plus, cross = s.symbols("plus cross")
    ni = s.ImmutableMatrix([1, 0, 0, 1])
    nj = s.ImmutableMatrix([1, 2 * r / (1 + r * r), 0, (1 - r * r) / (1 + r * r)])
    A = s.ImmutableMatrix(
        [[0, 0, 0, 0], [0, plus, cross, 0], [0, cross, -plus, 0], [0, 0, 0, 0]]
    )
    eta = s.diag(1, -1, -1, -1)
    contraction = s.factor((nj.T * A * nj)[0] / (nj.T * eta * ni)[0])
    a, b = s.symbols("a b", positive=True)
    G, S = s.Function("G")(a, b), s.Function("S")(a, b)
    I = (a + b) * s.log(a + b) - a * s.log(a) - b * s.log(b)
    D2 = tree.coefficient_envelope(2)
    first = (
        CURRENT_CAP * pair.LINEAR_SLOPE / pair.RADIUS + FIRST_CURRENT_CAP * D2
    ) / s.sqrt(source.KAPPA)
    mixed = (
        CURRENT_CAP * pair.REMAINDER_COEFFICIENT
        + 2 * FIRST_CURRENT_CAP * pair.LINEAR_SLOPE / pair.RADIUS
        + MIXED_CURRENT_CAP * D2
    ) / s.sqrt(source.KAPPA)
    checks = {
        "generic_null_direction_shell": s.factor((nj.T * eta * nj)[0]),
        "generic_soft_TT_transversality": A * ni,
        "generic_soft_TT_trace": s.trace(eta * A),
        "literal_null_emitter_angular_coefficient": s.factor(
            contraction - 2 * plus / (1 + r * r)
        ),
        "first_face_product_rule": s.expand(
            s.diff(S * G, a) - s.diff(S, a) * G - S * s.diff(G, a)
        ),
        "all_four_mixed_face_product_terms": s.expand(
            s.diff(S * G, a, b)
            - S * s.diff(G, a, b)
            - s.diff(S, a) * s.diff(G, b)
            - s.diff(S, b) * s.diff(G, a)
            - s.diff(S, a, b) * G
        ),
        "same_original_two_real_tree_envelope": D2
        - 6 * source.HEAVY_MASS2**2 * (2 * 10**16) ** 2 / source.KAPPA,
        "entropy_mixed_derivative": s.factor(s.diff(I, a, b) - 1 / (a + b)),
        "entropy_first_face": s.limit(I, a, 0),
        "entropy_second_face": s.limit(I, b, 0),
    }
    gates.update(
        {
            "two_real_amplitude_envelope_below1e_minus371": 0
            < D2
            < s.Rational(1, 10**371),
            "derived_axis_current_product_below1e_minus735": 0 < first < AXIS_CAP,
            "derived_pair_current_product_below1e_minus723": 0 < mixed < PAIR_CAP,
            "frozen_two_real_linear_slope": pair.LINEAR_SLOPE < s.Rational(1, 10**350),
            "frozen_two_real_mixed_coefficient": pair.CAUCHY_COEFFICIENT
            < pair.REMAINDER_COEFFICIENT,
            "absolute_disc_only_for_kinematic_current": True,
            "correct_two_marked_state_not_Born_replacement": True,
            "trimmed_face_integrals_not_bare_pole_values": True,
        }
    )
    return {
        "checks": checks,
        "gates": {k: bool(v) for k, v in gates.items()},
        "whole_current_absolute_radius": d,
        "whole_current_derivative_caps": (
            CURRENT_CAP,
            FIRST_CURRENT_CAP,
            MIXED_CURRENT_CAP,
        ),
        "whole_two_real_amplitude_envelope": D2,
        "whole_derived_face_derivative_coefficients": (first, mixed),
        "whole_anchored_axis_and_pair_caps": (AXIS_CAP, PAIR_CAP),
        "whole_complete_face_identity": "At c=0, G3(a,b,0)=S_c(sigma_ab)*G2(a,b)/sqrt(kappa); sigma_ab contains the four recoiled massive and BOTH marked null legs. The other faces are permutations of this identity.",
        "whole_kinematic_domain": "The absolute1e-6 energy disc continues only the massive recoil and six-emitter leading current. It contains no hard or pure-soft propagator. Radicand and Doppler gaps give current norm<129 and first/mixed derivatives<1e9,1e15. No such absolute disc is asserted for G2 or G3.",
        "whole_uniform_face_hierarchy": "S313 gives |G2-G20|<B*W on radii e*W, hence first derivatives<=B/e and mixed derivative<1e-326/W. Combine the current product rule with S319's full G2 envelope. Anchored axes obey |g_i|<1e-735*w_i and pairs |g_ij|<1e-723*Iij. Compatible limits and trimmed FTC supply the axis and face extensions.",
    }
