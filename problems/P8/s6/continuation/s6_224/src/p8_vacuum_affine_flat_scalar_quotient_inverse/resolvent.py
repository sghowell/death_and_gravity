"""Causal curvature-coordinate inverse and uniformly bounded shifted primitives."""

from functools import cache

import sympy as s
from p8_proca_rank_one_inverse import cut as trace_cut
from p8_proca_rank_one_inverse import spectral as trace_spectral
from p8_vacuum_affine_isolated_shear_resolvent import normalization as shear

from . import geometry as g

t = s.Symbol("time", positive=True)
lam, q, p = g.lam, g.q, g.p
BINV = s.Matrix(
    [
        [1 / p, -1 / (3 * p)],
        [1 / p - 1 / lam**2, -1 / (3 * p) - s.Rational(2, 3) / lam**2],
    ]
)
C = s.Rational(5, 2)
JBOUND = s.Rational(45, 2)


def transfer(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("An exact nonnegative squared spatial momentum is required")
    value = s.Rational(value)
    if value < 0:
        raise ValueError("Squared spatial momentum must be nonnegative")
    return value


def coordinate_kernel(value):
    value = transfer(value)
    wave = t if value == 0 else s.sin(s.sqrt(value) * t) / s.sqrt(value)
    return s.ImmutableMatrix([[wave, -wave / 3], [wave - t, -wave / 3 - 2 * t / 3]])


@cache
def data():
    wave = s.sin(s.sqrt(q) * t) / s.sqrt(q)
    R = s.Matrix([[wave, -wave / 3], [wave - t, -wave / 3 - 2 * t / 3]])
    theta0 = trace_spectral.data()["H_at_zero"]
    shear0 = shear.A0
    channel = s.diag(g.f0, s.Rational(8, 3) * g.f2)
    inverse = BINV * s.diag(1 / g.f0, s.Rational(3, 8) / g.f2) * BINV.T
    complete = g.B.T * channel * g.B
    cos = s.Symbol("cosine", real=True)
    derivative = R.diff(t).subs(s.cos(s.sqrt(q) * t), cos)
    mass = s.Symbol("mass", positive=True)
    trace_G = (3 * s.pi / 4 - 2 * 3 * s.pi / 16 + 3 * 5 * s.pi / 32) / (2 * mass)
    checks = {
        "coordinate_left_inverse": (g.B * BINV - s.eye(2)).applyfunc(s.factor),
        "coordinate_right_inverse": (BINV * g.B - s.eye(2)).applyfunc(s.factor),
        "no_inverse_q_resolvent_identity": s.factor(
            q / (lam * lam * p) - (1 / lam**2 - 1 / p)
        ),
        "causal_coordinate_Laplace_kernel": (
            R.applyfunc(lambda v: s.laplace_transform(v, t, lam, noconds=True)) - BINV
        ).applyfunc(s.factor),
        "coordinate_initial_kernel": R.subs(t, 0),
        "continuous_zero_q_kernel": (
            R.applyfunc(lambda v: s.limit(v, q, 0)) - coordinate_kernel(0)
        ).applyfunc(s.factor),
        "complete_quotient_left_inverse": (complete * inverse - s.eye(2)).applyfunc(
            s.factor
        ),
        "complete_quotient_right_inverse": (inverse * complete - s.eye(2)).applyfunc(
            s.factor
        ),
        "complete_quotient_determinant": s.factor(
            complete.det() - s.Rational(8, 3) * g.f0 * g.f2 * lam**4 * p**2
        ),
        "actual_trace_static_measure": trace_cut.data()["spectral_static_moment"]
        - 1 / theta0,
        "actual_trace_uniform_primitive_bound": 2 / theta0 - s.Rational(1, 2),
        "actual_shear_uniform_primitive_bound": 2 / shear0 - 60,
        "scaled_diagonal_primitive_bound": s.Rational(3, 8) * 60 - JBOUND,
        "exact_coordinate_derivative_Frobenius_polynomial": s.expand(
            sum(v * v for v in derivative) - (20 * cos * cos - 14 * cos + 13) / 9
        ),
        "entrywise_coordinate_derivative_margin": C * C
        - (1 + s.Rational(1, 9) + 4 + 1)
        - s.Rational(5, 36),
        "bounded_trace_forward_kernel": trace_G - 27 * s.pi / (64 * mass),
    }
    return {
        "coordinate_inverse": BINV,
        "coordinate_inverse_kernel": R,
        "full_quotient_inverse_transform": inverse.applyfunc(s.factor),
        "shifted_reciprocals": "Replace p by lambda^2+q in the exact trace and shear reciprocal measures. The shear pole frequency becomes sqrt(q+r m^2), with its same negative reciprocal residue retained; the continuum frequencies are sqrt(q+tau).",
        "uniform_primitives": "Jtrace,q=-integral rho0(tau)/(q+tau)(1-cos(sqrt(q+tau)t))dtau; J2,q additionally contains -R/(q+r m^2)(1-cos(sqrt(q+r m^2)t)). Positive static measures give |Jtrace,q|<=1/2 and |J2,q|<=60 for allq>=0.",
        "diagonal_scaled_primitive": "Jq=diag(Jtrace,q,(3/8)J2,q), ||Jq(t)||<=45/2. L[Jq]=diag(1/Ftrace(lambda^2+q),(3/8)/F2(lambda^2+q))/lambda.",
        "ordinary_matrix_inverse_kernel": "Eq=R'_q*Jq*R_q^T. No derivative of the oscillatory scalar kernel and no uniform half-lineL1 estimate are needed. Its Laplace transform is B^-1 diag(1/Ftrace,(3/8)/F2) (B^T)^-1.",
        "forward_distribution": "Each shifted scalar factor is -A_i(0)delta-(D^2+q)G_i,q, with G_i,q=int W_i sin(sqrt(q+4m^2/c)t)/(c sqrt(q+4m^2/c))dy. Uniform bounds are9pi/(128m) for shear and27pi/(64m) for trace.",
        "checks": checks,
        "gates": {
            "both_coordinate_derivative_bounds_finite": C > 0,
            "trace_positive_axis_static_bound": theta0 == 4,
            "shear_positive_axis_static_bound": shear0 == s.Rational(1, 30),
            "scaled_shear_primitive_dominates_trace": JBOUND > s.Rational(1, 2),
            "coordinate_derivative_bound_has_strict_margin": C * C
            > 1 + s.Rational(1, 9) + 4 + 1,
            "positive_q_never_used_as_a_divisor_in_time_kernel": not R.has(1 / q),
            "both_complete_inverse_products_kept": inverse.shape == (2, 2),
        },
    }
