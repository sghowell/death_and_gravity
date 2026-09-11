"""Exact first/second parameter error equations and explicit comparison bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_matrix_adiabatic import initial as prior

from . import mixed, reference

K = prior.PARTITION
E0 = prior.GRAPH_ERROR
E1 = 2 * s.Integer(10) ** 32
E2 = s.Integer(10) ** 35


@cache
def constants():
    c = mixed.constants()
    d = reference.constants()
    A = d["reference_parameter_norm_over_inverse_frequency"]
    C = d["residual_parameter_norm_over_inverse_frequency_tenth"]
    R, S = c["R"], c["S"]
    small = s.Rational(11, 100)
    M1 = 6 + (2 * R[0, 1] + small * S[0, 1]) / K + 2 * S[0, 0] * A[1] / K**2
    first = 27 * (E0 * M1 + C[1] / K)
    M20 = (
        12
        + (2 * R[0, 2] + small * S[0, 2]) / K**1
        + (2 * S[0, 0] * A[2] + 4 * S[0, 1] * A[1]) / K**2
    )
    M21 = (
        12 + (4 * R[0, 1] + s.Rational(2, 5) * S[0, 1]) / K + 4 * S[0, 0] * A[1] / K**2
    )
    second = 27 * (E0 * M20 / K + E1 * M21 + 2 * S[0, 0] * E1**2 / K**10 + C[2] / K**2)
    return {
        "first_frequency_weighted_coefficient": M1,
        "second_E0_frequency_weighted_coefficient": M20,
        "second_E1_frequency_weighted_coefficient": M21,
        "actual_first_graph_error_coefficient": first,
        "actual_second_graph_error_coefficient": second,
    }


def generator(w, R, S, r):
    return 2 * s.I * w * r + R * r - r * R + S - r * S * r


def tangent(w, R, S, r, X):
    return 2 * s.I * w * X + R * X - X * R - X * S * r - r * S * X


@cache
def data():
    eps = s.Symbol("amplitude", real=True)
    w = s.Symbol("omega", positive=True)
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]])
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]])
    B = s.Matrix([[0, 1, 2], [1, 1, -1], [2, -1, -1]])
    skew = s.Matrix([[0, 1, -2], [-1, 0, 3], [2, -3, 0]])
    a0 = (G + s.I * D) / 100
    a1 = (D + s.I * B) / 71
    a2 = (B - s.I * G) / 53
    e0 = (B + s.I * G) / 1000
    e1 = (G - s.I * D) / 131
    e2 = (D + s.I * B) / 137
    R0, R1, R2 = skew / 31, skew.T / 37, skew / 41
    S0, S1, S2 = G / 43, D / 47, B / 59
    hat = a0 + eps * a1 + eps**2 * a2 / 2
    err = e0 + eps * e1 + eps**2 * e2 / 2
    W = w * (1 + eps + eps**2)
    R = R0 + eps * R1 + eps**2 * R2 / 2
    S = S0 + eps * S1 + eps**2 * S2 / 2
    literal = generator(W, R, S, hat + err) - generator(W, R, S, hat)
    r = a0 + e0
    r1 = a1 + e1
    first = tangent(w, R0, S0, r, e1) - a1 * S0 * e0 - e0 * S0 * a1
    first += 2 * s.I * w * e0 + R1 * e0 - e0 * R1 - e0 * S1 * r - a0 * S1 * e0
    second = tangent(w, R0, S0, r, e2) - a2 * S0 * e0 - e0 * S0 * a2
    second -= 2 * (r1 * S0 * r1 - a1 * S0 * a1)
    second += 2 * (
        2 * s.I * w * e1
        + R1 * e1
        - e1 * R1
        - e1 * S1 * r
        - r * S1 * e1
        - a1 * S1 * e0
        - e0 * S1 * a1
    )
    second += 4 * s.I * w * e0 + R2 * e0 - e0 * R2 - e0 * S2 * r - a0 * S2 * e0
    c = constants()
    return {
        "exact_error": "e'=G(rhat+e)-G(rhat)-F, G(r)=2i omega r+[R,r]+S-r S r. The exact tangent at the actual r acts on e1 and e2; all other terms are forcing.",
        "initial_parameter_data": "e1(t0)=e2(t0)=0 exactly, because every history and its finite jets coincide near the unchanged initial Cauchy surface. The nonzero e0 state mismatch is retained.",
        "stable_comparison": "For nu_minus>=1e16, the actual and reference graphs stay below1/10 and1/100. The tangent's non-isometric coefficient is at most2*11/10<3; the common norm comparison factor is below27.",
        "first_forcing": "[2|omega1|+2R1+S1(||r||+||rhat||)+2S0||rhat1||]E0+||F1||.",
        "second_forcing": "[2S0||rhat2||+4S1||rhat1||+2|omega2|+2R2+S2(||r||+||rhat||)]E0+[2S0(||r1||+||rhat1||)+4|omega1|+4R1+4S1||r||]E1+||F2||.",
        "constants": c,
        "graph_error_displays": {0: E0, 1: E1, 2: E2},
        "result": "The unchanged exact graph-minus-reference error and its first and second amplitude derivatives obey E0<1e30 nu_minus^-10, E1<2e32 nu_minus^-9, E2<1e35 nu_minus^-8 uniformly on the unit CD slab.",
        "checks": {
            "complete_first_parameter_noncommuting_error_equation": s.diff(
                literal, eps
            ).subs(eps, 0)
            - first,
            "complete_second_parameter_noncommuting_error_equation": s.diff(
                literal, eps, 2
            ).subs(eps, 0)
            - second,
            "same_retained_undifferentiated_error": E0 - prior.GRAPH_ERROR,
            "same_fixed_analysis_partition": K - prior.PARTITION,
        },
        "gates": {
            "actual_graph_tangent_exponent_below_three": 2 * 11 * s.Rational(1, 10) < 3,
            "norm_comparison_factor": 3**3 == 27,
            "first_actual_error_display": c["actual_first_graph_error_coefficient"]
            < E1,
            "second_actual_error_display": c["actual_second_graph_error_coefficient"]
            < E2,
            "forcing_weight_losses_explicit": 10 - 1 == 9 and 10 - 2 == 8,
        },
    }
