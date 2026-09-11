"""Finite matrix Riccati reference, uniform residual and exact-state comparison."""

from functools import cache
from math import comb

import sympy as s

from . import jets


def _three_derivative_terms(order):
    for a in range(order + 1):
        for b in range(order - a + 1):
            c = order - a - b
            yield a, b, c, comb(order, a) * comb(order - a, b)


def reciprocal_jet(omega):
    v = [1 / omega[0]]
    for n in range(1, len(omega)):
        v.append(-v[0] * sum(comb(n, j) * omega[j] * v[n - j] for j in range(1, n + 1)))
    return v


def reference_jets(omega, rotation, squeeze, order=jets.ORDER):
    if type(order) is not int or not 1 <= order <= jets.ORDER:
        raise ValueError("Require a native finite Riccati order1..10")
    if min(len(omega), len(rotation), len(squeeze)) < order + 1:
        raise ValueError(
            "The complete raw derivative jets through the order are required"
        )
    inv = reciprocal_jet(omega[: order + 1])
    zero = s.zeros(squeeze[0].rows)
    ref = {
        1: [
            s.I
            * sum((comb(d, a) * inv[a] * squeeze[d - a] for a in range(d + 1)), zero)
            / 2
            for d in range(order + 1)
        ]
    }
    for n in range(1, order):

        def right(d, n=n):
            value = ref[n][d + 1]
            value -= sum(
                (
                    comb(d, a)
                    * (rotation[a] * ref[n][d - a] - ref[n][d - a] * rotation[a])
                    for a in range(d + 1)
                ),
                zero,
            )
            for j in range(1, n):
                l = n - j
                value += sum(
                    (
                        factor * ref[j][a] * squeeze[b] * ref[l][c]
                        for a, b, c, factor in _three_derivative_terms(d)
                    ),
                    zero,
                )
            return value

        ref[n + 1] = [
            -s.I
            * sum((comb(d, a) * inv[a] * right(d - a) for a in range(d + 1)), zero)
            / 2
            for d in range(order - n + 1)
        ]
    return ref


@cache
def constants():
    c = jets.constants()
    v, R, S = (
        c[key]
        for key in ("omega_reciprocal_relative", "R_derivatives", "S_derivatives")
    )
    N = jets.ORDER
    bounds = {
        1: [
            sum(comb(d, a) * v[a] * S[d - a] for a in range(d + 1)) / 2
            for d in range(N + 1)
        ]
    }
    for n in range(1, N):

        def right(d, n=n):
            value = bounds[n][d + 1] + 2 * sum(
                comb(d, a) * R[a] * bounds[n][d - a] for a in range(d + 1)
            )
            for j in range(1, n):
                l = n - j
                value += sum(
                    factor * bounds[j][a] * S[b] * bounds[l][c]
                    for a, b, c, factor in _three_derivative_terms(d)
                )
            return value

        bounds[n + 1] = [
            sum(comb(d, a) * v[a] * right(d - a) for a in range(d + 1)) / 2
            for d in range(N - n + 1)
        ]
    size = sum(bounds[n][0] / jets.MASS**n for n in bounds)
    residual = bounds[N][1] + 2 * R[0] * bounds[N][0]
    residual += S[0] * sum(
        bounds[j][0] * bounds[l][0] / jets.MASS ** (j + l - N)
        for j in bounds
        for l in bounds
        if j + l >= N
    )
    return {
        "coefficient_derivative_bounds": bounds,
        "reference_norm_at_mass_floor": size,
        "residual_over_inverse_frequency_tenth": residual,
    }


@cache
def data():
    w = [s.Integer(1000), s.Integer(2), s.Integer(-3), s.Integer(5)]
    R0 = s.Matrix([[0, 1, -2], [-1, 0, 3], [2, -3, 0]]) / 7
    S0 = s.Matrix([[1, 2, -1], [2, -3, 1], [-1, 1, 2]]) / 11
    D = s.diag(1, -1, 0) / 13
    R = [R0, 2 * R0, -R0, 3 * R0]
    S = [S0, D, 2 * S0 - D, -S0 + D]
    ref = reference_jets(w, R, S, 3)
    total = sum((ref[n][0] for n in ref), s.zeros(3))
    derivative = sum((ref[n][1] for n in ref), s.zeros(3))
    direct = (
        derivative
        - 2 * s.I * w[0] * total
        - R[0] * total
        + total * R[0]
        - S[0]
        + total * S[0] * total
    )
    tail = ref[3][1] - R[0] * ref[3][0] + ref[3][0] * R[0]
    tail += sum(
        (ref[j][0] * S[0] * ref[l][0] for j in ref for l in ref if j + l >= 3),
        s.zeros(3),
    )
    c = constants()
    return {
        "formal_recurrence": "r1=-S/(2i omega); r_(n+1)=(r_n'-[R,r_n]+sum_(j+l=n)r_j S r_l)/(2i omega), with all noncommuting products and raw time derivative jets retained.",
        "finite_order": jets.ORDER,
        "uniform_reference_constants": c,
        "residual": "F_N=r_Ntotal'-2i omega r_Ntotal-[R,r_Ntotal]-S+r_Ntotal S r_Ntotal. Its norm is bounded by the displayed C times omega^-10; the statement concerns a finite reference, not convergence of the formal series.",
        "exact_comparison": "The scalar fast phase and skew rotation generate norm-preserving two-sided propagation. For the actual pure Gaussian graph r, ||r||<1; if ||r_Ntotal||<=1/100 then ||r-r_Ntotal||(t)<=1e6(||initial mismatch||+C nu_gamma^-10) on the unit slab, where nu_gamma=min_I omega.",
        "checks": {
            "complete_nonlinear_noncommuting_reference_residual": direct - tail,
            "first_reference_is_symmetric": ref[1][0] - ref[1][0].T,
            "second_reference_is_symmetric": ref[2][0] - ref[2][0].T,
            "third_reference_is_symmetric": ref[3][0] - ref[3][0].T,
        },
        "gates": {
            "reference_small_at_every_actual_momentum": c[
                "reference_norm_at_mass_floor"
            ]
            < s.Rational(1, 100),
            "finite_residual_coefficient_positive": c[
                "residual_over_inverse_frequency_tenth"
            ]
            > 0,
            "explicit_residual_coefficient_display": c[
                "residual_over_inverse_frequency_tenth"
            ]
            < 37 * s.Integer(10) ** 22,
            "comparison_exponent_below_twelve": (1 + s.Rational(1, 100)) * 11 < 12,
            "exact_exponential_display": 3**12 < 10**6,
        },
    }
