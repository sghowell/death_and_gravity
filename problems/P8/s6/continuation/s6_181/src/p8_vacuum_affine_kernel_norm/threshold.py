"""Continuous endpoint derivative bounds, with exact rational arithmetic."""

from functools import cache

import sympy as s
from p8_proca_rank_one_inverse import cut

x, z, t = s.symbols("x z t", positive=True)
C_SMALL = s.Integer(200000000)


@cache
def data():
    C = 3 - 2 * z + 3 * z * z
    A = s.Rational(16, 15) + z - 3 * z * z
    G = z * C
    zz = 1 - 4 / (x + 2) ** 2
    c = s.sqrt(x + 4) / (x + 2)
    lower = s.Rational(16, 15) + z * (
        3 * (z - s.Rational(5, 6)) ** 2 + s.Rational(23, 12)
    )
    f = 1 / (1 - z * t * t)
    Q1 = s.Integer(3168)
    Q2 = s.Integer(50066)
    inv2 = 2 * Q1**2 + Q2
    b0, b1, b2 = s.Integer(2), 6 + 2 * Q1, 21 + 12 * Q1 + 2 * inv2
    scaled = (b0, b0 / 2 + 2 * b1, b0 / 4 + 2 * b1 + 4 * b2)
    pi_quotient, pi_remainder = s.div(t**4 * (1 - t) ** 4, 1 + t * t, t)
    checks = {
        "same_parent_density_D": s.factor(
            cut.data()["real_positive_H_on_upper_cut"].subs(cut.z, z)
            - (A + s.sqrt(z) * C * cut.ell)
        ),
        "same_parent_density_U": cut.data()["positive_imaginary_H_over_pi"].subs(
            cut.z, z
        )
        - s.sqrt(z) * C / 2,
        "threshold_fraction": s.factor(zz - x * (x + 4) / (x + 2) ** 2),
        "threshold_z_first": s.factor(s.diff(zz, x) - 8 / (x + 2) ** 3),
        "threshold_z_second": s.factor(s.diff(zz, x, 2) + 24 / (x + 2) ** 4),
        "threshold_sqrt_factor": s.factor(c * c - zz / x),
        "threshold_c_upper_completion": s.expand((x + 2) ** 2 - (x + 4) - x * (x + 3)),
        "positive_D_gap_completion": s.expand(
            lower - (s.Rational(16, 15) + 4 * z - 5 * z * z + 3 * z**3)
        ),
        "f_first_integrand": s.diff(f, z) - t * t / (1 - z * t * t) ** 2,
        "f_second_integrand": s.diff(f, z, 2) - 2 * t**4 / (1 - z * t * t) ** 3,
        "C_pair_positive_completion": s.expand(
            C - (s.Rational(8, 3) + 3 * (z - s.Rational(1, 3)) ** 2)
        ),
        "G_first_product": s.diff(G, z) - C - z * s.diff(C, z),
        "G_second_product": s.diff(G, z, 2) - 2 * s.diff(C, z) - z * s.diff(C, z, 2),
        "V_first_product": s.diff(z * C * C / 4, z)
        - (C * C + 2 * z * C * s.diff(C, z)) / 4,
        "V_second_product": s.diff(z * C * C / 4, z, 2)
        - (4 * C * s.diff(C, z) + 2 * z * (s.diff(C, z) ** 2 + C * s.diff(C, z, 2)))
        / 4,
        "D0_envelope": s.Integer(4) + 3 * 4 - 16,
        "Dz_envelope": s.Integer(6) + 10 * 4 + 3 * 16 - 94,
        "Dzz_envelope": s.Integer(6) + 21 * 4 + 2 * 10 * 16 + 3 * 128 - 794,
        "Dxx_envelope": s.Integer(794) + 2 * 94 - 982,
        "Vx_envelope": (s.Integer(16) + 2 * s.Rational(3, 4) * 4 * 8) / 4 - 16,
        "Vzz_envelope": (s.Integer(4) * 4 * 8 + 2 * s.Rational(3, 4) * (8**2 + 4 * 6))
        / 4
        - 65,
        "Vxx_envelope": s.Integer(65) + 2 * 16 - 97,
        "Qx_envelope": s.Integer(2) * 16 * 94 + 10 * 16 - Q1,
        "Qxx_envelope": s.Integer(2) * 94**2 + 2 * 16 * 982 + 10 * 97 - Q2,
        "inverse_Q_second_envelope": inv2 - 20122514,
        "b_first_envelope": b1 - 6342,
        "b_second_envelope": b2 - 40283065,
        "weighted_a_first_envelope": scaled[1] - 12685,
        "weighted_a_second_envelope": scaled[2] - s.Rational(322289889, 2),
        "positive_pi_integrand_division": s.expand(
            t**4 * (1 - t) ** 4 - (1 + t * t) * pi_quotient - pi_remainder
        ),
        "pi_integral_identity": s.integrate(pi_quotient, (t, 0, 1))
        + pi_remainder * s.pi / 4
        - (s.Rational(22, 7) - s.pi),
    }
    return {
        "interval": "0<x<=2, z=1-4/(x+2)^2 in [0,3/4]",
        "density": "a(x)=sqrt(x)*b(x), b=c(x)*C(z)/(2Q), Q=D^2+pi^2*z*C(z)^2/4",
        "threshold_c": c,
        "f_integral": "atanh(sqrt(z))/sqrt(z)=integral_0^1 (1-z*t^2)^-1 dt, continued at z=0",
        "f_derivative_envelopes": (s.Integer(4), s.Integer(16), s.Integer(128)),
        "D_envelopes": {"D": 16, "Dz": 94, "Dzz": 794, "Dx": 94, "Dxx": 982},
        "Q_derivative_envelopes": (Q1, Q2),
        "inverse_Q_derivative_envelopes": (s.Integer(1), Q1, inv2),
        "b_derivative_envelopes": (b0, b1, b2),
        "weighted_a_envelopes_over_sqrt_x": scaled,
        "common_threshold_constant": C_SMALL,
        "checks": checks,
        "gates": {
            "all_scaled_endpoint_constants_below_common": all(
                v < C_SMALL for v in scaled
            ),
            "c_first_triangle_below_one": s.Rational(1, 8) + s.Rational(3, 4) < 1,
            "c_second_triangle_below_one": s.Rational(1, 64)
            + s.Rational(1, 8)
            + s.Rational(3, 4)
            < 1,
            "pi_rational_upper_squared_below_ten": s.Rational(22, 7) ** 2 < 10,
            "D_gap_strictly_above_one": s.Rational(16, 15) > 1,
        },
    }
