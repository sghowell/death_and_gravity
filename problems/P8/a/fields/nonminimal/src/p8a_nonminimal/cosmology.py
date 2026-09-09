"""Fresh scalar, finite-scheme and field-strength costs in a verified history tube."""

from fractions import Fraction
from functools import cache

import sympy as sp

from . import stress

RATIO = sp.Rational(1, 100)
PAST = (sp.Rational(21, 10), sp.Integer(9), sp.Integer(70), sp.Integer(800))
FUTURE = (sp.Integer(4), sp.Integer(128), sp.Integer(16384), sp.Integer(1048576))
DELTA_MAX = sp.Rational(1, 10**8)
ZETA_MAX = sp.Rational(1, 5000)
SIGMA_MAX = sp.Integer(5)
COST_UPPER = sp.Integer(5000000)


def exact_nonnegative(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require an exact nonnegative rational, not binary floats")
    value = sp.Rational(value)
    if value < 0:
        raise ValueError("Require a nonnegative exact value")
    return value


@cache
def costs():
    r = RATIO
    d, c = PAST, FUTURE
    qd = sp.Rational(3, 4) * d[0] ** 2 + sp.Rational(3, 2) * d[1]
    qc = sp.Rational(3, 4) * c[0] ** 2 + sp.Rational(3, 2) * c[1]
    p = (
        sp.Rational(7, 2)
        + sp.Rational(11, 5) * d[0] * r
        + sp.Rational(2, 3) * qd * r * r
    ) ** 2 / r**3
    f = (sp.Rational(7, 2) * (1 + 2 * c[0]) + sp.Rational(25, 12) * qc) ** 2
    C0 = (
        sp.Rational(7, 18) * (p + f)
        + sp.Rational(13, 12600) * r * stress.loss_polynomial(d, 0)
        + sp.Rational(13, 1080) * stress.loss_polynomial(c, 0)
    )
    Cb = (
        sp.Rational(13, 12600)
        * r
        * (stress.loss_polynomial(d, 1) - stress.loss_polynomial(d, 0))
    )
    Cb += sp.Rational(13, 1080) * (
        stress.loss_polynomial(c, 1) - stress.loss_polynomial(c, 0)
    )
    Wp = (sp.Rational(11, 10) + d[0] * r) ** 2 / r
    Wf = (sp.Rational(11, 10) + sp.Rational(27, 20) * c[0]) ** 2
    return {
        "scalar_C0": C0,
        "scalar_Cbeta": Cb,
        "past_Wick_square_cost": Wp,
        "future_Wick_square_cost": Wf,
        "Wick_square_cost": Wp + Wf,
        "past_caps": PAST,
        "future_caps": FUTURE,
        "actual_past_ratio": r,
        "positive_C0_margin": COST_UPPER - C0,
        "positive_Cbeta_margin": COST_UPPER - Cb,
    }


def gate(weighted_delta, zeta, sigma):
    weighted_delta, zeta, sigma = map(exact_nonnegative, (weighted_delta, zeta, sigma))
    if weighted_delta > DELTA_MAX or zeta > ZETA_MAX or sigma > SIGMA_MAX:
        raise ValueError(
            "The input exceeds the named quantum, field-strength or source gate"
        )
    h = sp.Rational(19, 10)
    gain = 3 * h + sp.Rational(39, 35) * h * h * RATIO
    margin = (
        gain
        - sp.Rational(18, 5)
        - COST_UPPER * weighted_delta
        - costs()["Wick_square_cost"] * zeta
        - sp.Rational(13, 35) * (1 + RATIO) * sigma
    )
    return {
        "weighted_delta": weighted_delta,
        "zeta_upper_kappa_Phi_square_over_three": zeta,
        "sigma": sigma,
        "strict_focusing_margin": margin,
        "margin_above_nine_fiftieths": margin - sp.Rational(9, 50),
        "future_geometry_and_state_amplitude_are_conditional_hypotheses_not_verified_from_past": True,
    }


@cache
def checks():
    s = sp.Symbol("s", real=True)
    p = 3 * s * s - 2 * s**3
    t = sp.Symbol("t", real=True)
    H = sp.Function("H")(t)
    f = sp.Function("f")(t)
    Ric = 3 * (sp.diff(H, t) + H * H)
    r, d0, c0 = sp.symbols("ratio past_H_cap future_H_cap", positive=True)
    wp = (sp.Rational(11, 10) + d0 * r) ** 2 / r
    wf = (sp.Rational(11, 10) + sp.Rational(27, 20) * c0) ** 2
    exact_past = (
        sp.Rational(6, 5) / r
        + sp.Rational(3, 2) * d0
        + sp.Rational(117, 140) * d0 * d0 * r
    )
    exact_future = (
        sp.Rational(6, 5) + sp.Rational(27, 10) * c0 + sp.Rational(9, 5) * c0 * c0
    )
    return {
        "geometric_square_identity_without_pointwise_SEC": sp.expand(
            Ric * f * f
            - 3 * (sp.diff(f, t) - H * f) ** 2
            + 3 * sp.diff(f, t) ** 2
            - sp.diff(3 * H * f * f, t)
        ),
        "cubic_state_weight_relative_endpoint_moment": sp.integrate(
            (p / s) ** 2, (s, 0, 1)
        )
        - sp.Rational(4, 5),
        "cubic_cross_relative_endpoint_moment": sp.integrate(
            sp.diff(p, s) * p / s, (s, 0, 1)
        )
        - sp.Rational(9, 10),
        "past_state_cost_dominates_exact_polynomial_upper": sp.factor(
            wp
            - exact_past
            - (
                sp.Rational(1, 100) / r
                + sp.Rational(7, 10) * d0
                + sp.Rational(23, 140) * d0 * d0 * r
            )
        ),
        "future_state_cost_dominates_exact_polynomial_upper": sp.factor(
            wf
            - exact_future
            - (
                sp.Rational(1, 100)
                + sp.Rational(27, 100) * c0
                + sp.Rational(9, 400) * c0 * c0
            )
        ),
        "fresh_scalar_C0_constant": costs()["scalar_C0"]
        - sp.Rational(373108140471263, 75000000),
        "fresh_scalar_Cbeta_constant": costs()["scalar_Cbeta"]
        - sp.Rational(3345309847373, 10500000),
        "fresh_Wick_square_cost_constant": costs()["Wick_square_cost"]
        - sp.Rational(1679141, 10000),
        "actual_nonminimal_focusing_margin": gate(DELTA_MAX, ZETA_MAX, SIGMA_MAX)[
            "strict_focusing_margin"
        ]
        - sp.Rational(63325013, 350000000),
        "actual_nonminimal_strict_margin_above_point_eighteen": gate(
            DELTA_MAX, ZETA_MAX, SIGMA_MAX
        )["margin_above_nine_fiftieths"]
        - sp.Rational(325013, 350000000),
    }
