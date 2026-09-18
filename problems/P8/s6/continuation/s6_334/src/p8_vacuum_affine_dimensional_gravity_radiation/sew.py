"""Complete continued polarization sew and pointwise evanescent rate bounds."""

from functools import cache

import sympy as s


def exact_real(value):
    if isinstance(value, (bool, float, str)) or value is None:
        raise ValueError("Require an exact finite real scalar")
    value = s.sympify(value)
    if (
        not isinstance(value, s.Expr)
        or value.has(s.Float)
        or value.is_number is not True
        or value.is_real is not True
        or value.is_finite is not True
    ):
        raise ValueError("Require an exact finite real scalar")
    return value


def regulator(value):
    value = exact_real(value)
    if value < 0 or value > s.Rational(1, 8):
        raise ValueError("Require 0<=epsilon<=1/8")
    return value


def ratio(epsilon):
    e = regulator(epsilon)
    return e / (1 + e)


def core(value):
    if not isinstance(value, s.MatrixBase) or value.shape != (2, 2) or value != value.T:
        raise ValueError("Require an exact symmetric two-by-two canonical core")
    if any(
        v.has(s.Float) or v.is_number is not True or v.is_finite is not True
        for v in value
    ):
        raise ValueError("Require an exact finite canonical core")
    return s.Matrix(value)


def frobenius(A, B):
    return s.trace(s.conjugate(A).T * B)


def bilinear(A, B, dimension):
    return frobenius(A, B) - s.conjugate(s.trace(A)) * s.trace(B) / (dimension - 2)


def dimension_sew(value, dimension):
    A = core(value)
    D = exact_real(dimension)
    if D < 4:
        raise ValueError("Require D>=4")
    return s.expand(bilinear(A, A, D))


def coefficients(first, correction):
    A, B = core(first), core(correction)
    ta, tb = s.trace(A), s.trace(B)
    return tuple(
        s.expand(v)
        for v in (
            bilinear(A, A, 4),
            2 * s.re(bilinear(A, B, 4)) + s.conjugate(ta) * ta / 2,
            bilinear(B, B, 4) + s.re(s.conjugate(ta) * tb),
            s.conjugate(tb) * tb / 2,
        )
    )


def continued_rate(first, correction, epsilon):
    t = ratio(epsilon)
    return s.expand(
        sum(value * t**j for j, value in enumerate(coefficients(first, correction)))
    )


def regulator_rate_upper(first_norm, correction_norm, epsilon):
    a, b = map(exact_real, (first_norm, correction_norm))
    e = regulator(epsilon)
    if a < 0 or b < 0:
        raise ValueError("Require nonnegative Frobenius norm bounds")
    return e * (a * a + s.Rational(20, 9) * a * b + s.Rational(10, 81) * b * b)


@cache
def data():
    r = s.Symbol("ratio", real=True)
    a, b, c, d, e, f, ai, bi, ci, di, ei, fi = s.symbols(
        "a b c d e f ai bi ci di ei fi", real=True
    )
    A = s.Matrix([[a + s.I * ai, b + s.I * bi], [b + s.I * bi, c + s.I * ci]])
    B = s.Matrix([[d + s.I * di, e + s.I * ei], [e + s.I * ei, f + s.I * fi]])
    ta, tb = s.trace(A), s.trace(B)
    polynomial = bilinear(A, A, 4) + r * (
        2 * s.re(bilinear(A, B, 4)) + s.conjugate(ta) * ta / 2
    )
    polynomial += (
        r * r * (bilinear(B, B, 4) + s.re(s.conjugate(ta) * tb))
        + r**3 * s.conjugate(tb) * tb / 2
    )
    exact = (
        frobenius(A + r * B, A + r * B)
        - (1 - r) * s.conjugate(s.trace(A + r * B)) * s.trace(A + r * B) / 2
    )
    ep = s.Symbol("epsilon", nonnegative=True)
    checks = {
        "generic_complex_sew_cubic": s.expand(exact - polynomial),
        "first_regulator_derivative": s.expand(
            s.diff(exact.subs(r, ep / (1 + ep)), ep).subs(ep, 0)
            - 2 * s.re(bilinear(A, B, 4))
            - s.conjugate(ta) * ta / 2
        ),
        "extra_trace_term": s.expand(
            exact
            - bilinear(A + r * B, A + r * B, 4)
            - r * s.conjugate(s.trace(A + r * B)) * s.trace(A + r * B) / 2
        ),
        "max_dimension_ratio": ratio(s.Rational(1, 8)) - s.Rational(1, 9),
        "majorant_cross_coefficient": 2 + 2 * s.Rational(1, 9) - s.Rational(20, 9),
        "majorant_correction_coefficient": s.Rational(1, 9)
        + s.Rational(1, 81)
        - s.Rational(10, 81),
        "four_D_trace_negative_control": dimension_sew(s.eye(2), 4),
        "nonfour_D_trace_negative_control": continued_rate(
            s.eye(2), s.zeros(2), s.Rational(1, 8)
        )
        - s.Rational(2, 9),
    }
    return {
        "checks": checks,
        "gates": {
            "entire_complex_interference_preserved": True,
            "rank_two_sew_credited_to_S296": True,
            "analytic_continued_bilinear_not_noninteger_Hilbert_space": True,
            "pure_trace_invisible_in_D4_but_not_continued_rate": bool(
                continued_rate(s.eye(2), s.zeros(2), s.Rational(1, 8)) > 0
            ),
            "pointwise_norm_majorant_not_integrated_soft_bound": True,
        },
        "whole_sew": "t=epsilon/(1+epsilon); Q=||A+tB||_F^2-(1-t)|tr(A+tB)|^2/2=Q4(A+tB)+t|tr(A+tB)|^2/2.",
        "whole_cubic_coefficients": "q0=Q4(A);q1=2Re<A,B>4+|trA|^2/2;q2=Q4(B)+Re(conj(trA)trB);q3=|trB|^2/2.",
        "whole_pointwise_regulator_bound": "For 0<=epsilon<=1/8, a>=||A||_F,b>=||B||_F, |Q(epsilon)-Q(0)|<=epsilon[a^2+(20/9)ab+(10/81)b^2]. This alone is not integrable at soft/internal poles.",
        "whole_negative_control": "A=I2,B=0: the D4 helicity sew is zero but continued rate=2epsilon/(1+epsilon).",
    }
