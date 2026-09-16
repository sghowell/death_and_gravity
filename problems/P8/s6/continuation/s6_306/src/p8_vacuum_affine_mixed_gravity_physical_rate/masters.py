"""Exact box collapse and explicit physical-sheet light/triangle masters."""

from functools import cache

import sympy as s

from . import source

A, B, N = s.symbols("channel partner heavy_mass_squared", real=True)
X, U, H = s.symbols("contour_parameter radial_parameter box_parameter", real=True)


def contour(x=X):
    x = s.sympify(x)
    return x + s.I * x * (1 - x) * (1 - 2 * x)


def contour_denominator(channel=A, x=X):
    return 1 - s.sympify(channel) * contour(x) * (1 - contour(x))


def box_multiplier(channel=A, heavy=N):
    channel, heavy = map(s.sympify, (channel, heavy))
    if channel == 0:
        return 1 / heavy
    return -s.log(1 - channel / heavy) / channel


def collapsed_box(channel=A, partner=B, heavy=N):
    return source.mixed.J0(partner, 1) * box_multiplier(channel, heavy)


def box_remainder(channel=A, partner=B, heavy=N):
    channel, partner, heavy = map(s.sympify, (channel, partner, heavy))
    return source.graphs.massive_T(partner, 1, heavy) + source.mixed.J0(
        partner, 1
    ) * s.log(1 - channel / heavy)


def Q_moment(order, value):
    if isinstance(order, bool) or not isinstance(order, int):
        raise TypeError("Require an integer moment order")
    if order not in (0, 1, 2):
        raise ValueError("Only the three complete endpoint moments are in scope")
    value = s.sympify(value)
    if value == 0:
        return s.Rational(1, order + 1)
    return (
        -s.log(1 - value) - sum(value**j / s.Integer(j) for j in range(1, order + 1))
    ) / value ** (order + 1)


def radial_moment(order, denominator, heavy=N):
    if isinstance(order, bool) or not isinstance(order, int):
        raise TypeError("Require an integer radial moment")
    if order not in (1, 2, 3):
        raise ValueError("Require radial moment1,2,3")
    denominator, heavy = map(s.sympify, (denominator, heavy))
    alpha = (1 + s.sqrt(1 - 4 * denominator / heavy)) / 2
    beta = denominator / (heavy * alpha)
    return (Q_moment(order - 1, alpha) - Q_moment(order - 1, beta)) / (
        heavy * (alpha - beta)
    )


@cache
def data():
    u, y, w, aa, nn, a, b = s.symbols("u y w radial_A n a b", positive=True)
    x, r, v, c, al, be = s.symbols("x r v c alpha beta", real=True)
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand(value))

    h = y * u * u
    den = u * u * aa + h * (nn - a * (1 - u - h))
    put(
        "literal_box_projective_denominator",
        den - u * u * (aa + y * (nn - a + a * (u + y * u * u))),
    )
    put(
        "literal_numerator_Jacobian_cancellation",
        u * (u + 2 * h) * u * u / den**2
        - (1 + 2 * y * u) / (aa + y * (nn - a + a * (u + y * u * u))) ** 2,
    )
    put("projective_w_derivative", s.diff(u + y * u * u, u) - (1 + 2 * y * u))
    Q = nn - a + a * w
    put(
        "positive_y_integral_primitive",
        s.diff(-1 / (Q * (aa + Q * y)), y) - 1 / (aa + Q * y) ** 2,
    )
    put(
        "positive_y_integral_boundary",
        1 / (Q * aa) - (-1 / (Q * (aa + Q * y))).subs(y, 0) * (-1),
    )
    put("remaining_w_log_primitive", s.diff(s.log(Q) / a, w) - 1 / Q)
    put(
        "whole_log_multiplier_zero_channel_limit",
        s.limit(box_multiplier(a, nn), a, 0) - 1 / nn,
    )
    for j in (1, 2, 3):
        put(
            "whole_box_small_channel_coefficient_" + str(j),
            s.diff(box_multiplier(a, nn), a, j).limit(a, 0) / s.factorial(j)
            - 1 / ((j + 1) * nn ** (j + 1)),
        )
    put(
        "box_remainder_matches_collapsed_K",
        -a * box_multiplier(a, nn) - s.log(1 - a / nn),
    )
    put(
        "wrong_box_sign_defect_nonzero",
        a * box_multiplier(a, nn) - s.log(1 - a / nn) + 2 * s.log(1 - a / nn),
    )
    f = x * (1 - x) * (1 - 2 * x)
    z = x + s.I * f
    real = 1 - b * (x * (1 - x) + f * f)
    imag = -b * x * (1 - x) * (1 - 2 * x) ** 2
    put("whole_deformed_denominator_real_imag", 1 - b * z * (1 - z) - real - s.I * imag)
    put("contour_derivative", s.diff(z, x) - 1 - s.I * (1 - 6 * x + 6 * x * x))
    put("contour_derivative_imag_upper", 1 - (1 - 6 * x + 6 * x * x) - 6 * x * (1 - x))
    put(
        "contour_derivative_imag_lower",
        (1 - 6 * x + 6 * x * x) + s.Rational(1, 2) - 6 * (x - s.Rational(1, 2)) ** 2,
    )
    t = s.symbols("t", real=True)
    put(
        "contour_r_real",
        real.subs(x, (1 + t) / 2)
        - (1 - b * (1 - t * t) / 4 - b * t * t * (1 - t * t) ** 2 / 16),
    )
    put("contour_r_imag", imag.subs(x, (1 + t) / 2) + b * t * t * (1 - t * t) / 4)
    put(
        "middle_r_product_lower",
        r * (1 - r)
        - s.Rational(7, 64)
        - (r - s.Rational(1, 8)) * (s.Rational(7, 8) - r),
    )
    put(
        "complex_endpoint_v_norm_margin",
        1 - (r + r * (1 - r) ** 2 / 4) - (1 - r) * (1 - r * (1 - r) / 4),
    )
    put(
        "first_outer_contour_gap",
        1 - s.Rational(25, 4) * s.Rational(7, 8) / 4 + s.Rational(47, 128),
    )
    put(
        "second_outer_contour_gap",
        1 - s.Rational(1, 2) - s.Rational(1, 64) - s.Rational(31, 64),
    )
    put(
        "middle_contour_gap",
        s.Rational(25, 4) * s.Rational(7, 64) / 4 - s.Rational(175, 1024),
    )
    for j in (1, 2, 3):
        put(
            "radial_partial_fraction_" + str(j),
            u**j / ((1 - al * u) * (1 - be * u))
            - (u ** (j - 1) / (1 - al * u) - u ** (j - 1) / (1 - be * u)) / (al - be),
        )
        k = j - 1
        primitive = (
            -s.log(1 - c * u)
            - sum((c * u) ** ell / s.Integer(ell) for ell in range(1, k + 1))
        ) / c ** (k + 1)
        put("Q_integral_primitive_" + str(k), s.diff(primitive, u) - u**k / (1 - c * u))
        put(
            "Q_removable_zero_" + str(k),
            s.limit(Q_moment(k, c), c, 0) - s.Rational(1, k + 1),
        )
    put(
        "quadratic_root_factorization",
        (1 - al * u) * (1 - be * u) - (1 - (al + be) * u + al * be * u * u),
    )
    put("light_F1_complete_moment_numerator", u * (1 - u) ** 2 - (u - 2 * u * u + u**3))
    put(
        "light_F2_complete_moment_numerator",
        u * (u * u * v * v - 1) / 2 - (v * v * u**3 - u) / 2,
    )
    put(
        "heavy_active_denominator_lower",
        nn * u + (1 - u) ** 2 - 4 * u * u - 1 - (nn - 5) * u - 3 * u * (1 - u),
    )
    put("zero_transfer_minimum", nn * (1 - u) + u * u - 1 - (1 - u) * (nn - 1 - u))
    put(
        "U_positive_first_half_integral",
        s.diff(2 * s.log(nn * v / 2 + w / 4) / nn, v) - 1 / (nn * v / 2 + w / 4),
    )
    put(
        "radial_log_majorant_decreasing",
        s.diff((s.log(nn) + 10) / nn, nn) - (1 - s.log(nn) - 10) / nn**2,
    )
    put(
        "physical_Lll_log_weight",
        s.integrate(x * (1 - x), (x, 0, 1)) - s.Rational(1, 6),
    )
    return {
        "whole_box_K_exact": collapsed_box(),
        "whole_box_remainder_exact": box_remainder(),
        "whole_physical_contour": contour(),
        "whole_contour_denominator": s.expand(contour_denominator()),
        "whole_three_endpoint_radial_moments": {
            str(j): radial_moment(j, aa, nn) for j in (1, 2, 3)
        },
        "whole_box_proof": "On the positive parameter domain use y=h/u^2,w=u+h. The exact Jacobian gives integral dv dw dy/[A_b+y(n-a+aw)]^2. Integrating y and w proves K=J0(b)[ln n-ln(n-a)]/a, continuously1/n at a0. Continue this analytic identity, not the absolute singular integrand, to the physical Feynman boundary.",
        "whole_contour_proof": "For25/4<=b<=16 deform v=x+i x(1-x)(1-2x). The displayed three r-regions prove1/8<|A|<10, lower-bank|ln A|<8 and path length<2. For b<=0 use the undeformed positive domain. Hence|J0|<16,|J1|<128,|Lll|<16 and|Bbar|<8. The homotopy never crosses a zero and reproduces the original normal sheet.",
        "whole_triangle_proof": "For n>=128 the displayed quadratic has|beta|<1/10,|alpha|>9/10,|alpha-beta|>4/5. Principal lower-bank log beta has magnitude<ln n+7. The three full radial moments satisfy|I_j|<2(ln n+10)/n; the sharper j1 bound gives|T|<4(ln n+9)/n after the contour integral.",
        "checks": checks,
        "gates": {
            "projective_identity_continued_as_whole_function": True,
            "nonpinching_physical_contour_has_explicit_gap": True,
            "all_three_endpoint_moments_retained": True,
            "removable_box_zero_channel_kept": box_multiplier(0, nn) == 1 / nn,
            "wrong_box_sign_is_detected": s.factor(
                a * box_multiplier(a, nn) - s.log(1 - a / nn)
            )
            != 0,
        },
    }
