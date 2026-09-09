"""Native null mass-shell measure, oscillator normalization and actual stress."""

from functools import cache

import sympy as sp

EPSILON = sp.Rational(1, 1000)
WIDTH = sp.Rational(1, 10000)


@cache
def data():
    q, p2 = sp.symbols("q transverse_squared", positive=True)
    energy = (p2 + q * q) / (2 * q)
    kz = (p2 - q * q) / (2 * q)
    jac = -sp.diff(kz, q)
    e = sp.Symbol("epsilon", positive=True)
    s = sp.Symbol("affine_clock", real=True)
    C1, C2 = 4 - sp.I / e, -2 + sp.I / e
    v = C1 * sp.exp(-sp.I * e * s) + C2 * sp.exp(-2 * sp.I * e * s)
    x, y, xp, yp, xpp, ypp, xi = sp.symbols(
        "Re_v Im_v Re_v_prime Im_v_prime Re_v_second Im_v_second xi", real=True
    )
    n, m = sp.Rational(1, 3), -sp.Rational(2, 3)
    uv, duv = x + sp.I * y, xp + sp.I * yp
    wick = sp.expand(
        2 * n * uv * sp.conjugate(uv) + m * uv**2 + m * sp.conjugate(uv) ** 2
    )
    derivative = sp.expand(
        2 * n * duv * sp.conjugate(duv) + m * duv**2 + m * sp.conjugate(duv) ** 2
    )
    total_derivative = lambda f: (
        sp.diff(f, x) * xp
        + sp.diff(f, y) * yp
        + sp.diff(f, xp) * xpp
        + sp.diff(f, yp) * ypp
    )
    wick_second = sp.expand(total_derivative(total_derivative(wick)))
    null = sp.expand(derivative - xi * wick_second)
    return {
        "q": q,
        "transverse_squared": p2,
        "energy": energy,
        "longitudinal_momentum": kz,
        "positive_measure_jacobian": jac,
        "epsilon": e,
        "clock": s,
        "first_frequency_coefficient": C1,
        "second_frequency_coefficient": C2,
        "unmollified_profile": v,
        "occupation": n,
        "anomalous_occupation": m,
        "relative_Wick_square_per_mode_amplitude_squared": wick,
        "null_derivative_square_per_mode_amplitude_squared": derivative,
        "relative_Wick_second": wick_second,
        "actual_nonminimal_null_stress": null,
        "real_jets": (x, y, xp, yp, xpp, ypp),
        "xi": xi,
        "checks": {
            "actual_positive_mass_shell": sp.factor(energy**2 - kz**2 - p2),
            "actual_affine_null_frequency": sp.factor(energy - kz - q),
            "actual_null_coordinate_measure": sp.factor(
                jac / (2 * energy) - 1 / (2 * q)
            ),
            "positive_squeezed_covariance_uncertainty_saturated": (
                n + sp.Rational(1, 2)
            )
            ** 2
            - m * m
            - sp.Rational(1, 4),
            "actual_positive_position_variance": n
            + sp.Rational(1, 2)
            + m
            - sp.Rational(1, 6),
            "actual_positive_momentum_variance": n
            + sp.Rational(1, 2)
            - m
            - sp.Rational(3, 2),
            "profile_value_is_two": sp.simplify(v.subs(s, 0) - 2),
            "profile_first_derivative_is_one": sp.simplify(
                sp.diff(v, s).subs(s, 0) - 1
            ),
            "disjoint_bump_mode_norm_keeps_both_frequencies": sp.simplify(
                e * sp.expand_complex(C1 * sp.conjugate(C1))
                + 2 * e * sp.expand_complex(C2 * sp.conjugate(C2))
                - (24 * e + 3 / e)
            ),
            "actual_Wick_square_keeps_negative_quadrature": sp.expand(
                wick - sp.Rational(2, 3) * (-x * x + 3 * y * y)
            ),
            "actual_null_derivative_square_keeps_same_covariance": sp.expand(
                derivative - sp.Rational(2, 3) * (-xp * xp + 3 * yp * yp)
            ),
            "full_nonminimal_null_stress_keeps_improvement": sp.expand(
                null
                - sp.Rational(2, 3) * (1 - 2 * xi) * (-xp * xp + 3 * yp * yp)
                - sp.Rational(4, 3) * xi * (x * xpp - 3 * y * ypp)
            ),
            "real_affine_profile_negative_minimal_stress": sp.expand(
                null.subs({x: 2 + s, y: 0, xp: 1, yp: 0, xpp: 0, ypp: 0, xi: 0})
                + sp.Rational(2, 3)
            ),
            "real_affine_profile_negative_conformal_stress": sp.expand(
                null.subs(
                    {
                        x: 2 + s,
                        y: 0,
                        xp: 1,
                        yp: 0,
                        xpp: 0,
                        ypp: 0,
                        xi: sp.Rational(1, 6),
                    }
                )
                + sp.Rational(4, 9)
            ),
        },
    }


def profile_bounds():
    e, width = EPSILON, WIDTH
    # Triangle coefficients, then a real Taylor integral on every |s|<=1.
    second = e * e * (4 + 1 / e) + 4 * e * e * (2 + 1 / e)
    discrete_value, discrete_first = 3 + second / 2, 1 + second
    error0 = second / 2 + width**2 * discrete_value / 2
    error1 = second + width**2 * discrete_first / 2 + width**2 * discrete_value
    error2 = second + 2 * width**2 * discrete_first + width**2 * discrete_value
    err = sp.Rational(1, 100)
    negative_derivative = (1 - err) ** 2 - 3 * err * err
    improvement = (3 + err) * err + 3 * err * err
    wick_upper = -sp.Rational(2, 3) * negative_derivative
    uniform_null_upper = -sp.Rational(1, 3) * negative_derivative + improvement / 3
    conformal_null_upper = (
        -sp.Rational(4, 9) * negative_derivative + sp.Rational(2, 9) * improvement
    )
    return {
        "epsilon": e,
        "bump_width": width,
        "discrete_second_derivative_bound": second,
        "mollified_value_error": error0,
        "mollified_first_error": error1,
        "mollified_second_upper": error2,
        "common_jet_error_bound": err,
        "strict_error_margins": [err - v for v in (error0, error1, error2)],
        "Wick_square_upper_per_mode_amplitude_squared": wick_upper,
        "null_stress_upper_all_xi_zero_to_quarter": uniform_null_upper,
        "null_stress_upper_at_conformal_xi": conformal_null_upper,
        "negative_Wick_margin_below_minus_three_fifths": -sp.Rational(3, 5)
        - wick_upper,
        "negative_uniform_null_margin_below_minus_three_tenths": -sp.Rational(3, 10)
        - uniform_null_upper,
        "negative_conformal_null_margin_below_minus_two_fifths": -sp.Rational(2, 5)
        - conformal_null_upper,
        "positive_low_frequency_support_margin": e - width,
        "disjoint_frequency_support_margin": e - 2 * width,
    }


@cache
def checks():
    out = dict(data()["checks"])
    b = profile_bounds()
    out.update(
        {
            "exact_second_derivative_bound": b["discrete_second_derivative_bound"]
            - sp.Rational(1253, 250000),
            "negative_Wick_square_known_answer": b[
                "Wick_square_upper_per_mode_amplitude_squared"
            ]
            + sp.Rational(1633, 2500),
            "uniform_nonminimal_negative_null_known_answer": b[
                "null_stress_upper_all_xi_zero_to_quarter"
            ]
            + sp.Rational(4747, 15000),
            "conformal_negative_null_known_answer": b[
                "null_stress_upper_at_conformal_xi"
            ]
            + sp.Rational(9646, 22500),
        }
    )
    return out
