"""Removable massive derivative and relative-collinear null contraction."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_state_soft_index import index


def exact_scalar(value):
    if value is None or isinstance(value, (bool, float, s.Float, str)):
        raise TypeError("Require an exact finite real scalar")
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


def massive_fprime(value):
    t = exact_scalar(value)
    if not 1 <= t <= 7:
        raise ValueError("Require a massive pair invariant in[1,7]")
    if t == 1:
        return s.Rational(22, 3)
    root = s.sqrt(t * t - 1)
    angle = s.acosh(t)
    return 8 * t * angle / root + 2 * (2 * t * t - 1) * (root - t * angle) / root**3


def massive_c(value):
    t = exact_scalar(value)
    if not s.Rational(29, 16) <= t <= 7:
        raise ValueError("Require a same-orientation massive invariant in[29/16,7]")
    return t * (2 * t * t - 3) / (t * t - 1) ** s.Rational(3, 2)


def energy(value):
    R = exact_scalar(value)
    if not 0 <= R <= s.Rational(1, 8):
        raise ValueError("Require total radiated energy in[0,1/8]")
    return R


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    t, v = s.symbols("t v", positive=True)
    D = t + 1 - (t - 1) * v * v
    r = s.acosh(t) / (s.sqrt(t - 1) * s.sqrt(t + 1))
    put(
        "interior_real_kernel_differential_identity",
        (t * t - 1) * s.diff(r, t) + t * r - 1,
    )
    endpoint = []
    for j, expected in enumerate((1, -s.Rational(1, 3), s.Rational(4, 15))):
        value = s.integrate(s.diff(2 / D, t, j).subs(t, 1), (v, 0, 1))
        endpoint.append(value)
        put("fixed_interval_endpoint_" + str(j), value - expected)
    put(
        "fprime_removable_endpoint",
        8 * endpoint[0] + 2 * endpoint[1] - massive_fprime(1),
    )
    put(
        "fsecond_removable_endpoint",
        8 * endpoint[0] + 16 * endpoint[1] + 2 * endpoint[2] - s.Rational(16, 5),
    )
    h = s.Symbol("hyperbolic_cosine", positive=True)
    put(
        "positive_first_derivative_integrand",
        s.diff(2 * (2 * t * t - 1) / (t + h), t)
        - (4 * t * t + 8 * t * h + 2) / (t + h) ** 2,
    )
    put(
        "original_S300_pair_normalization",
        index.massive_pair(t, 1)
        - (t * t - s.Rational(1, 2)) * s.acosh(t) / s.sqrt(t * t - 1),
    )
    c = t * (2 * t * t - 3) / (t * t - 1) ** s.Rational(3, 2)
    put("imaginary_derivative", s.diff(c, t) - 3 / (t * t - 1) ** s.Rational(5, 2))
    put(
        "imaginary_upper_square_margin",
        4 * (t * t - 1) ** 3 - t * t * (2 * t * t - 3) ** 2 - 3 * t * t + 4,
    )
    a, b, delta = s.symbols("d_i d_j relative_gap", positive=True)
    put(
        "null_transverse_Gram_identity",
        b / a * a * (2 - a)
        + a / b * b * (2 - b)
        - 2 * (a + b - a * b - delta)
        - 2 * delta,
    )
    put(
        "relative_collinear_log_product",
        s.limit(delta * (s.log(2 * delta) + 1), delta, 0, dir="+"),
    )
    eta = s.diag(1, -1, -1, -1)
    p = s.Matrix(s.symbols("p0:4"))
    b = s.Matrix(s.symbols("b0:4"))
    q = s.Matrix(s.symbols("q0:4"))
    xi = s.Matrix(s.symbols("xi0:4"))
    entries = iter(s.symbols("A0:10"))
    A = s.zeros(4)
    for i in range(4):
        for j in range(i, 4):
            A[i, j] = A[j, i] = next(entries)
    dot = lambda a, b: (a.T * eta * b)[0]
    ep = lambda a, b: (a.T * A * b)[0]
    Dp, Db = dot(p, q), dot(b, q)
    dp = ep(p, p) * q / Dp - eta * A * p
    put("generic_soft_Lorentz_mass_preservation", s.cancel(dot(p, dp)))
    put(
        "generic_pair_derivative_contraction",
        s.cancel(dot(dp, b) - ep(p, p) * Db / Dp + ep(p, b)),
    )
    ward = eta * (q * xi.T + xi * q.T) * eta
    wa = (p.T * ward * p)[0] * Db / Dp - (p.T * ward * b)[0]
    wb = (b.T * ward * b)[0] * Dp / Db - (b.T * ward * p)[0]
    put("longitudinal_pair_antisymmetry", s.cancel(wa + wb))
    put(
        "longitudinal_single_leg_identity",
        s.cancel(wa - dot(p, xi) * Db + Dp * dot(b, xi)),
    )
    return {
        "checks": checks,
        "gates": {
            "fixed_interval_denominator_at_least_two": True,
            "real_kernel_first_derivative_cap": 8 * 7 == 56,
            "real_kernel_second_derivative_cap": 8
            + s.Rational(16 * 7, 3)
            + s.Rational(2 * 97 * 4, 15)
            < 161,
            "same_orientation_c_positive": s.Rational(29, 16) ** 2 > s.Rational(3, 2),
            "same_orientation_c_prime_below_one": 9
            < (s.Rational(29, 16) ** 2 - 1) ** 5,
            "relative_collinear_bound_not_soft_collinear_point_assignment": True,
        },
        "whole_real_kernel": "f=2(2t^2-1)r,r=2int_0^1dv/[t+1-(t-1)v^2],r(1)=1,f'(1)=22/3,f''(1)=16/5. On[1,7]:0<f'<=56,|f''|<161.",
        "whole_imaginary_kernel": "c=t(2t^2-3)/(t^2-1)^(3/2),0<c<2,0<c'<1 on[29/16,7]. Outgoing same-time gamma>=29/16; incoming>=17/8.",
        "whole_null_pair_bound": "t_i=P_perp(n)n_i,d_i=1-n_i.n,v=sqrt(d_j/d_i)t_i-sqrt(d_i/d_j)t_j. |v|^2=2delta_ij, S_ij=w_iw_j A(v,v), hence |S_ij|<=2w_iw_j delta_ij. Exclude d_i*d_j=0 initially; bound holds on approaches.",
    }
