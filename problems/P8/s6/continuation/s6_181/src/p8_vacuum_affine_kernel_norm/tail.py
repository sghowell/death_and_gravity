"""Continuous high-frequency derivative majorants, not asymptotic samples."""

from functools import cache

import sympy as s

r, v, x = s.symbols("r v x", positive=True)
C_TAIL = s.Integer(20000000)


@cache
def data():
    A = s.Rational(16, 15) + v * v - 3 * v**4
    B = v * (3 - 2 * v * v + 3 * v**4)
    z = s.symbols("z", positive=True)
    q = s.Function("q")(r)
    f = s.Function("a")(r)
    q1, q2 = s.Integer(2040), s.Integer(29016)
    ar = 12 + 2 * q1
    arr = 60 + 2 * 12 * q1 + 2 * q2 + 4 * q1**2
    scaled2 = ar + arr
    partial_e2 = sum(s.Rational(2) ** j / s.factorial(j) for j in range(6))
    checks = {
        "tail_A_lower_completion": s.expand(
            A + s.Rational(14, 15) - (1 - v * v) * (3 * v * v + 2)
        ),
        "tail_B_monotone_completion": s.expand(
            s.diff(B, v) - (15 * (v * v - s.Rational(1, 5)) ** 2 + s.Rational(12, 5))
        ),
        "tail_B_lower_completion": s.expand(
            B - s.Rational(8, 3) * v - 3 * v * (v * v - s.Rational(1, 3)) ** 2
        ),
        "tail_B_endpoint": B.subs(v, 1) - 4,
        "tail_D_r_triangle": s.Integer(14) + 4 - 18,
        "tail_D_rr_triangle": s.Integer(66) + 2 * 24 - 114,
        "tail_A_rr_triangle": s.Integer(38) + 2 * 14 - 66,
        "tail_B_rr_triangle": s.Integer(72) + 2 * 24 - 120,
        "tail_Q_r_ratio": s.Integer(2) * 60 + 2 * 10 * 2 * 12 * 4 - q1,
        "tail_Q_rr_ratio": s.Integer(2) * 60**2
        + 2 * 348
        + 2 * 10 * (12**2 + 2 * 60) * 4
        - q2,
        "tail_a_r_envelope": ar - 4092,
        "tail_a_rr_envelope": arr - 16753452,
        "tail_x2_a_second_envelope": scaled2 - 16757544,
        "tail_exp_partial_sum": partial_e2 - s.Rational(109, 15),
        "tail_tanh_lower_fraction": (s.Integer(7) - 1) / (7 + 1) - s.Rational(3, 4),
        "tail_frequency_cosine": s.expand((1 + x / 2) ** 2 - x * (x + 4) / 4 - 1),
        "tail_q_derivative": s.factor(
            s.diff(s.tanh(r / 2), r) - (1 - s.tanh(r / 2) ** 2) / 2
        ),
        "tail_x_first_chain": s.cancel(
            ((2 * s.cosh(r) - 2) / (2 * s.sinh(r)) - s.tanh(r / 2)).rewrite(s.exp)
        ),
        "tail_second_chain": s.expand(
            q * s.diff(q * s.diff(f, r), r)
            - q * s.diff(f, r)
            - (q * q * s.diff(f, r, 2) + q * (s.diff(q, r) - 1) * s.diff(f, r))
        ),
        "tail_log2_positive_integrand": s.factor(
            1 / z - s.Rational(1, 2) - (2 - z) / (2 * z)
        ),
    }
    return {
        "interval": "x>=1; r=arcosh(1+x/2)>=log(1+x)>=log(2)>1/2",
        "same_density": "a(x)=U/Q, U=B(tanh r)/2, D=A(tanh r)+B(tanh r)*r, Q=D^2+pi^2*U^2",
        "A": A,
        "B": B,
        "D_lower": "D>=r: use the common 16/15 gap for r<=1, and tanh(1)>3/4 for r>=1",
        "Dr_over_r_upper": s.Integer(60),
        "Drr_over_r_upper": s.Integer(348),
        "U_derivative_envelopes": (s.Integer(2), s.Integer(12), s.Integer(60)),
        "Q_ratio_derivative_envelopes": (q1, q2),
        "a_r_derivative_envelopes_times_r_squared": (s.Integer(2), ar, arr),
        "weighted_x_derivative_envelopes_times_r_squared": (s.Integer(2), ar, scaled2),
        "common_tail_constant": C_TAIL,
        "checks": checks,
        "gates": {
            "exp_two_partial_sum_exceeds_seven": partial_e2 > 7,
            "large_r_D_lower_exceeds_r": 2 - s.Rational(14, 15) > 1,
            "Dr_uniform_ratio": 18 * 2 + 24 <= 60,
            "Drr_uniform_ratio": 114 * 2 + 120 <= 348,
            "all_scaled_tail_constants_below_common": all(
                y < C_TAIL for y in (2, ar, scaled2)
            ),
        },
    }
