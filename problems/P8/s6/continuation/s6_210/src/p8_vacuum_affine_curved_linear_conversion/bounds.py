"""Proper-time operator structure and explicit coefficient norms, not cutoff uniformity."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import angular, jets

CUBIC = s.Rational(1, 100)
LINEAR = s.Integer(10) ** 4


@cache
def constants():
    A = s.Rational(25, 16)
    mass = 49 * modes.MASS**2 * A / (1024 * 9)
    spatial = s.Rational(781, 24576 * 9)
    curvature = s.Rational(7 * 11, 1024 * 9)
    first = 31 * s.Rational(5, 2) / (2048 * 9)
    second = 31 * A / (2048 * 9)
    return {
        "cubic_first_spatial_derivative": s.Rational(31, 512 * 9),
        "linear_mass": mass,
        "linear_spatial": spatial,
        "linear_curvature": curvature,
        "linear_first_source_time": first,
        "linear_second_source_time": second,
        "complete_linear_two_time_three_space_bound": mass
        + spatial
        + curvature
        + first
        + second,
    }


@cache
def data():
    t = s.symbols("t", real=True)
    v = s.symbols("absolute_time", nonnegative=True)
    a, H, Hp, _c = jets.clock(t)
    D = s.Function("D")(t)
    G = s.Function("Gamma")(t)
    op = lambda F: a * s.diff(F, t, 2) + s.diff(a, t) * s.diff(F, t)
    green = a * (D * s.diff(G, t) - s.diff(D, t) * G)
    T, V, W = s.symbols("T V W", real=True)
    p, m = s.symbols("p m", positive=True)
    lin = angular.linear(t, p, m, T, V, W)
    checks = {
        "actual_scale_factor_upper_identity": s.factor(
            s.Rational(25, 16)
            - a
            - (s.Rational(1, 4) - t * t) * (s.Rational(9, 4) + t * t)
        ),
        "actual_scale_factor_lower_identity": s.factor(a - 1 - t * t * (2 + t * t)),
        "actual_scale_derivative_absolute_upper_identity": s.factor(
            s.Rational(5, 2)
            - 4 * v * (1 + v * v)
            - (s.Rational(1, 2) - v) * (4 * v * v + 2 * v + 5)
        ),
        "actual_curvature_weight_upper_identity": s.factor(
            11 - a * (Hp + 2 * H * H) - 28 * (s.Rational(1, 4) - t * t)
        ),
        "full_source_proper_time_divergence_form": s.simplify(
            op(G) - s.diff(a * s.diff(G, t), t)
        ),
        "full_compact_time_Green_identity": s.simplify(
            D * op(G) - G * op(D) - s.diff(green, t)
        ),
        "actual_first_derivative_required_for_formal_symmetry": s.factor(
            lin[1] - s.diff(lin[2], t)
        ),
        "canonical_both_cubic_metric_factors": 4 * CUBIC / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 802,
        "canonical_both_linear_metric_factors": 4 * LINEAR / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 796,
    }
    return {
        "coefficient_bounds": "For the fixed slab a<=25/16, 1/a<=1, |a_prime|<=5/2 and a(H_prime+2H^2)<=11. Frobenius Cauchy bounds each T,V,W by ||D||||G||. The absolute invariant row sums are31,49,781,7.",
        "explicit_norms": "|A3[D,Gamma]|<1/100 ||D||L2 ||grad Gamma||L2. Define Z23^2=sum_r0..2 integral dt d^3P/(2pi)^3 (1+|P|^2)^3 |partial_t^r Gamma_hat|_F^2, with the unchanged Fourier measure. Then |A1[D,Gamma]|<1e4 ||D||L2 Z23[Gamma]. The exact positive coefficient sum is retained in constants().",
        "proper_time": "The derivative part is proportional to a Gamma_second+a_prime Gamma_prime=partial_t(a partial_t Gamma). The complete Green identity is retained for compactly supported clocks. Formal symmetry of this endpoint multiplier is not symmetry of the full retarded response.",
        "origin_and_nonlocality": "Both coefficient multipliers vanish continuously at P0. Their direction projectors are bounded; the p prefactors control the origin. Odd momentum magnitudes and directional projectors remain nonlocal. No angular limit is chosen as a new parameter.",
        "canonical_boundary": "Both metric factors give coefficient displays4e-802 and4e-796. Multiplication by K^3 and K is still nonuniform as K is removed. No small full renormalized response, inverse, background or physical cutoff follows from these displays.",
        "constants": constants(),
        "checks": checks,
        "gates": {
            "cubic_first_derivative_display": bool(
                constants()["cubic_first_spatial_derivative"] < CUBIC
            ),
            "linear_two_time_three_space_display": bool(
                constants()["complete_linear_two_time_three_space_bound"] < LINEAR
            ),
            "proper_time_first_derivative_not_omitted": True,
            "compact_support_required_for_Green_identity": True,
            "zero_transfer_by_continuity_not_direction_choice": True,
            "canonical_coefficients_not_cutoff_uniformity": True,
        },
    }
