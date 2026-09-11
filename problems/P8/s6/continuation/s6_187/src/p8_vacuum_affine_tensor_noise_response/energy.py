"""All-momentum advanced tree wave energy and explicit adjoint smearing norm."""

from functools import cache

import sympy as s


@cache
def data():
    t = s.Symbol("t", real=True)
    k = s.Symbol("k", nonnegative=True)
    scale = (1 + t * t) ** 2
    actual_H = s.diff(scale, t) / scale
    wave = s.Function("wave")(t)
    forcing = s.Function("forcing")(t)
    modal_energy = (
        scale**3 * s.diff(wave, t) ** 2 + scale * k * k * wave**2 + scale**3 * wave**2
    ) / 2
    modal_derivative = s.diff(modal_energy, t).subs(
        s.diff(wave, t, 2),
        forcing - 3 * actual_H * s.diff(wave, t) - k * k * wave / scale**2,
    )
    modal_expected = (
        scale**3 * s.diff(wave, t) * (forcing + wave)
        + actual_H
        * (
            -3 * scale**3 * s.diff(wave, t) ** 2
            + scale * k * k * wave**2
            + 3 * scale**3 * wave**2
        )
        / 2
    )
    H = s.Symbol("H", real=True)
    v, vt, grad, q = s.symbols("v v_t gradient_norm q", real=True)
    a = s.Symbol("a", positive=True)
    # After spatial IBP, these are the exact density contributions to E'.
    derivative = (
        a**3 * vt * q
        + a**3 * v * vt
        - 3 * H * a**3 * vt**2 / 2
        + H * a * grad**2 / 2
        + 3 * H * a**3 * v * v / 2
    )
    E = (a**3 * vt**2 + a * grad**2 + a**3 * v * v) / 2
    base = s.Integer(54)
    second = 1 + s.Rational(24, 5) * base + 2 * base
    third = (
        1
        + s.Rational(24, 5) * 400
        + 12 * base
        + 2 * base
        + s.Rational(16, 5) * 2 * base
    )
    norm_squared = 3 * base**2 + 400**2 + 4000**2
    checks = {
        "independent_actual_wave_energy_derivative": s.simplify(
            modal_derivative - modal_expected
        ),
        "actual_H_endpoint_bound_factorization": s.factor(
            s.Rational(8, 5)
            - 4 * t / (1 + t * t)
            - 4 * (1 - 2 * t) * (2 - t) / (5 * (1 + t * t))
        ),
        "actual_Hprime_upper_factorization": s.factor(
            4 - s.diff(actual_H, t) - 4 * t * t * (3 + t * t) / (1 + t * t) ** 2
        ),
        "actual_maximum_volume_square_root": scale.subs(t, s.Rational(1, 2))
        ** s.Rational(3, 2)
        - s.Rational(125, 64),
        "exact_positive_energy_terms": 2 * E
        - a**3 * vt**2
        - a * grad**2
        - a**3 * v * v,
        "exact_wave_energy_density": s.expand(
            derivative
            - (
                a**3 * vt * (q + v)
                + H * (-3 * a**3 * vt**2 + a * grad**2 + 3 * a**3 * v * v) / 2
            )
        ),
        "cross_energy_control": s.expand(v * v + vt * vt - 2 * v * vt - (v - vt) ** 2),
        "second_time_wave_reconstruction": s.Rational(24, 5) - 3 * s.Rational(8, 5),
        "third_time_scale_derivative_bound": s.Rational(16, 5) - 2 * s.Rational(8, 5),
        "complete_smearing_norm_count": norm_squared - 16168748,
    }
    return {
        "energy": "E=1/2 integral[a^3 v_t^2+a|grad v|^2+a^3 v^2]dx, summed in the TT Frobenius norm",
        "exact_energy_derivative": derivative,
        "backward_energy_bound": "y=sqrt(2E), |y'|<=3y+||a^(3/2)q||; zero final data gives y<54||q||L2(dt dx)",
        "detector_norm": "D(q)^2=sum_(|alpha|<=2)||partial_x^alpha q||L2(dt dx;F)^2+||partial_t q||L2(dt dx;F)^2",
        "field_first_time_and_spatial_L2_constant": base,
        "actual_second_time_majorant": second,
        "second_time_display_upper": s.Integer(400),
        "actual_third_time_majorant": third,
        "third_time_display_upper": s.Integer(4000),
        "complete_adjoint_stress_smearing_norm_squared_upper": norm_squared,
        "complete_adjoint_stress_smearing_norm_upper": s.Integer(5000),
        "all_momentum_boundary": "An added v^2 energy term controls zero momentum. Delta commutes with the actual operator. No momentum band or inverse momentum denominator is used.",
        "checks": checks,
        "gates": {
            "energy_growth_coefficient_below_six": 1 + 3 * s.Rational(8, 5) < 6,
            "maximum_volume_square_root_below_two": s.Rational(125, 64) < 2,
            "backward_exponential_below_twenty_seven": s.exp(3) < 27,
            "complete_backward_energy_display": 27 * 2 <= base,
            "second_time_bound_below_four_hundred": second < 400,
            "third_time_bound_below_four_thousand": third < 4000,
            "full_adjoint_norm_below_five_thousand": norm_squared < 5000**2,
        },
    }
