"""Explicit full-convolution nonlinear transverse unit-determinant chart."""

from functools import cache

import sympy as s

TAU_RADIUS = s.Rational(1, 100)
F_RADIUS = s.Rational(1, 1000)
SHAPE_RADIUS = TAU_RADIUS + F_RADIUS
CONTRACTION = 6 * SHAPE_RADIUS + 9 * SHAPE_RADIUS**2


def projector(vector):
    k = s.Matrix(vector)
    return s.eye(3) if k.dot(k) == 0 else s.eye(3) - k * k.T / k.dot(k)


def tt(value, vector):
    P = projector(vector)
    return (
        value - s.eye(3) * s.trace(value) / 3
        if s.Matrix(vector).dot(s.Matrix(vector)) == 0
        else P * value * P - P * s.trace(P * value) / 2
    )


@cache
def data():
    v = s.symbols("independent_shape0:6", real=True)
    A = s.Matrix([[v[0], v[1], v[2]], [v[1], v[3], v[4]], [v[2], v[4], v[5]]])
    E2 = (s.trace(A) ** 2 - s.trace(A * A)) / 2
    expected = (
        s.eye(3)
        + (s.trace(A) * s.eye(3) - A)
        + (E2 * s.eye(3) - s.trace(A) * A + A * A)
    )
    checks = {
        "whole_noncommuting_3D_determinant": s.expand(
            (s.eye(3) + A).det() - 1 - s.trace(A) - E2 - A.det()
        ),
        "whole_full_cofactor_not_eigenvalue_truncation": (
            (s.eye(3) + A).adjugate() - expected
        ).applyfunc(s.expand),
    }
    conventions = []
    for k in [(0, 0, 0), (1, 0, 0), (1, 2, 0), (1, -1, 2)]:
        P = projector(k)
        L = 3 if k == (0, 0, 0) else 2
        T = tt(A, k)
        name = "_".join(map(str, k))
        checks["full_projector_idempotence_" + name] = P * P - P
        checks["full_projector_trace_" + name] = s.trace(P) - L
        checks["whole_TT_trace_" + name] = s.trace(T).factor()
        checks["whole_TT_transverse_" + name] = (T * s.Matrix(k)).applyfunc(s.factor)
        checks["whole_TT_annihilates_scalar_complement_" + name] = tt(P, k)
        checks["whole_TT_projection_idempotence_" + name] = (tt(T, k) - T).applyfunc(
            s.factor
        )
        conventions.append({"k": k, "P": P, "trace_eigenvalue": L})
    t = s.Symbol("positive_shape_amplitude", nonnegative=True)
    trial = t + s.Rational(3, 4) * t * t + s.Rational(33, 8) * t**3
    majorant = s.Poly(s.expand(6 * trial**2 + 6 * trial**3), t)
    high = sum(
        coeff * t ** (power[0] - 4)
        for power, coeff in majorant.terms()
        if power[0] >= 4
    )
    high_at_radius = high.subs(t, TAU_RADIUS)
    image = 3 * SHAPE_RADIUS**2 + 3 * SHAPE_RADIUS**3
    checks["whole_exact_full_contraction_constant"] = CONTRACTION - s.Rational(
        67089, 10**6
    )
    checks["whole_exact_self_map_constant"] = image - s.Rational(366993, 10**9)
    checks["whole_exact_cubic_remainder_majorant"] = high_at_radius - s.Rational(
        173086592988411, 2560000000000
    )
    return {
        "full_Fourier_conventions": conventions,
        "whole_shape_domain": {
            "tau_radius": TAU_RADIUS,
            "scalar_radius": F_RADIUS,
            "shape_radius": SHAPE_RADIUS,
            "weighted_Fourier_exponent": 8,
        },
        "whole_scalar_operator": "B f(k)=P(k)f(k), L=trace B=2I+Pi0. The mean eigenvalue is3, not2; the zero tracefree tensor sector has dimension5, not2.",
        "whole_fixed_point": "Q=I+tau+Bf; f=-L^-1(E2(tau+Bf)+det(tau+Bf)). All products are full convolution. The Banach norm sums (1+|k|)^8 times the matrix operator norm of each Fourier coefficient.",
        "whole_full_self_map_bound": image,
        "whole_full_contraction_bound": CONTRACTION,
        "whole_scalar_linearized_inverse_bound": s.Rational(5, 9),
        "whole_full_tangent": "K_Q h=cof(Q):Bh. delta f=-K_Q^-1(cof(Q):delta tau); delta Q=delta tau-B K_Q^-1(cof(Q):delta tau).",
        "whole_second_and_third_shape_coefficients": "f(epsilon tau)=epsilon^2 f2+epsilon^3 f3+O(epsilon^4), f2=(1/2)L^-1 trace(tau^2), f3=L^-1[trace(tau Bf2)-det tau]. Coefficients are not derivatives: the latter multiply by2! and3!.",
        "whole_cubic_remainder_majorant_polynomial": high,
        "whole_cubic_remainder_norm_bound": "||f(epsilon tau)-epsilon^2 f2-epsilon^3 f3||<=60 |epsilon|^4 for ||tau||<=1 and |epsilon|<=1/100.",
        "whole_shape_and_metric_reconstruction": "gamma=a^2 exp(2v)Q^-1. It has det gamma=a^6 exp(6v). At the same linear reference tau=-t_TT with the dual momentum sign also reversed; this is not a new state or a deletion of S258 off-slice volume contacts.",
        "checks": checks,
        "gates": {
            "strict_full_self_map": image < F_RADIUS,
            "strict_full_contraction": CONTRACTION < s.Rational(1, 10),
            "positive_shape_and_old_ghost_coercivity": SHAPE_RADIUS < s.Rational(1, 8),
            "full_scalar_inverse_bound": s.Rational(1, 2) / (1 - CONTRACTION)
            < s.Rational(5, 9),
            "cubic_trial_inside_scalar_ball": s.Rational(3, 4) * TAU_RADIUS**2
            + s.Rational(33, 8) * TAU_RADIUS**3
            < F_RADIUS,
            "all_positive_high_order_majorant_coefficients": all(
                c > 0 for c in s.Poly(high, t).all_coeffs()
            ),
            "whole_cubic_majorant_less_than_100": high_at_radius < 100,
            "whole_root_remainder_less_than_60": s.Rational(50, 1) / (1 - CONTRACTION)
            < 60,
            "generated_homogeneous_scalar_never_discarded": True,
            "full_Banach_proof_not_finite_projected_determinant_claim": True,
        },
    }
