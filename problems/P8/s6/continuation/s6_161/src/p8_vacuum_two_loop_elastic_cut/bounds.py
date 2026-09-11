"""Exact rational full-amplitude and finite-window cut majorants."""

from fractions import Fraction
from functools import cache

import sympy as s


def exact(v):
    if isinstance(v, bool) or not isinstance(v, (int, Fraction, s.Rational)):
        raise TypeError("Require a finite exact rational")
    return s.Rational(v)


def scalar_enclosure(L, g, M, ell, log_cap, Q=144):
    L, g, M, ell, Q = map(exact, (L, g, M, ell, Q))
    if min(L, g, ell) < 0 or M < 32 or Q <= 0:
        raise ValueError(
            "Require nonnegative interactions/log allowance, M>=32 and Q>0"
        )
    if type(log_cap) is not int or log_cap < 1:
        raise TypeError("Require a positive native integer logarithm power cap")
    if M > 2**log_cap:
        raise ValueError("The declared heavy logarithm power cap fails")
    C = L + g / (M - 6)
    tri = (3 * log_cap + 22) / M
    box = (6 * log_cap + 52) / M**2
    bubble = ell + 6
    return {
        "channel_C_absolute_upper": C,
        "bubble_MS_absolute_upper_without_Q": bubble,
        "triangle_absolute_upper_without_Q": tri,
        "box_absolute_upper_without_Q": box,
        "full_scalar_MS_amplitude_absolute_upper": 3
        * (C * C * bubble + 4 * C * g * tri + 4 * g * g * box)
        / (2 * Q),
    }


def fermion_enclosure(m, Y, N=6, Q=144):
    m, Y, N, Q = map(exact, (m, Y, N, Q))
    if m < 24 or Y <= 0 or N <= 0 or Q <= 0:
        raise ValueError("Require m>=24, Y>0 and positive multiplicity/loop factor")
    return N * Y * Y * (64 + 3072 / m) / Q


def cut_enclosure(lam, B):
    lam, B = map(exact, (lam, B))
    if lam <= 0 or B < 0:
        raise ValueError("Require lambda>0 and a nonnegative full-amplitude allowance")
    return {
        "density_second_absolute_upper": 73 * lam * B / 80,
        "integrated_second_absolute_upper": 73 * lam * B / 1280,
        "integrated_first_upper": 3 * lam * lam,
    }


def improved_band(lam, computed_relative_error, B):
    lam, E, B = map(exact, (lam, computed_relative_error, B))
    if lam <= 0 or E < 0 or B < 0:
        raise ValueError("Require positive tree scale and nonnegative error allowances")
    cut = cut_enclosure(lam, B)
    total = (
        4 * lam * E
        + cut["integrated_first_upper"]
        + cut["integrated_second_absolute_upper"]
    )
    return {
        "tree": 4 * lam,
        "total_absolute_error": total,
        "total_relative_error": total / (4 * lam),
        "formal_uniform_lower": 4 * lam - total,
        "positive_formal_uniform_lower": bool(4 * lam - total > 0),
    }


@cache
def data():
    lam, B, E = s.symbols("lambda B1 Erelative", positive=True)
    C, g, ell, T, F, Q = s.symbols("C g ell triangle box Q", positive=True)
    u, v = s.symbols("u v", real=True)
    raw = (C + g * u + g * v) ** 2
    scalar = 3 * (C * C * (ell + 6) + 4 * C * g * T + 4 * g * g * F) / (2 * Q)
    return {
        "scalar_uniform_bound": scalar,
        "exact_log_cap": "M<=2^K implies log(M)<K because log(2)<1. The separately supplied ell is the already proved upper bound on log(mu^2); its actual value is imported from S6.133.",
        "full_one_loop_bound": "B1=Bscalar+Bferm+146 k0_upper lambda, because |A0|<73 lambda and the sole full first field correction is -2 k0 A0.",
        "checks": {
            "four_mixed_vertex_terms": s.diff(raw, C).subs({C: 0, u: 1, v: 1}) - 4 * g,
            "four_heavy_vertex_terms": raw.subs({C: 0, u: 1, v: 1}) - 4 * g * g,
            "three_channel_half_bubble_weight": scalar
            - 3 * (C * C * (ell + 6) + 4 * C * g * T + 4 * g * g * F) / (2 * Q),
            "two_loop_cut_tree_relative": s.Rational(73, 1280) * lam * B / (4 * lam)
            - s.Rational(73, 5120) * B,
            "cut_subtracted_formal_relative_budget": (
                4 * lam * E + 3 * lam * lam + s.Rational(73, 1280) * lam * B
            )
            / (4 * lam)
            - E
            - 3 * lam / 4
            - s.Rational(73, 5120) * B,
            "canonical_field_correction_absolute_multiplier": 2 * 73 - 146,
        },
        "scope": "All quantities are bounds on computed coefficients or their finite-window integrals. A positive formal band does not bound the omitted physical higher orders or establish a dispersive equality.",
    }
