"""Exact conformal/proper-time conversion without changing the Hessian measure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

a = s.Symbol("a", positive=True)
H, H1, H2, H3 = s.symbols("H H1 H2 H3", real=True)
q = s.Symbol("q", nonnegative=True)
X = s.symbols("x0:7", real=True)

HELICITY = {
    "tensor": (s.Integer(1), s.Integer(2), s.Integer(1), s.Integer(0)),
    "vector": (s.Integer(0), s.Integer(1), s.Integer(0), s.Integer(0)),
    "scalar": (s.Rational(-1, 3), s.Rational(2, 3), s.Rational(1, 9), s.Rational(2, 3)),
}


def dt(expr):
    out = (
        s.diff(expr, a) * a * H
        + s.diff(expr, H) * H1
        + s.diff(expr, H1) * H2
        + s.diff(expr, H2) * H3
    )
    return s.expand(out + sum(s.diff(expr, X[j]) * X[j + 1] for j in range(6)))


def eta(expr):
    return a * dt(expr)


def eta_power(expr, n):
    for _ in range(n):
        expr = eta(expr)
    return s.expand(expr)


def operators(kind):
    c, b, d, e = HELICITY[kind]
    S = H1 + 2 * H * H
    S1 = H2 + 4 * H * H1
    return {
        "R_old": -(a**3 * X[2] + 3 * a**3 * H * X[1] + a * c * q * X[0]) / 2,
        "R_old_squared": -6
        * (a**3 * S * X[2] + a**3 * (3 * H * S + S1) * X[1] + a * S * c * q * X[0])
        + 2 * e * q * q * X[0] / a,
        "Weyl_squared": a**3
        * (
            X[4]
            + 6 * H * X[3]
            + (4 * H1 + 11 * H * H) * X[2]
            + (H2 + 7 * H * H1 + 6 * H**3) * X[1]
        )
        + b * q * a * (X[2] + H * X[1])
        + d * q * q * X[0] / a,
    }


@cache
def data():
    checks = {}
    r = a * a * (H1 + 2 * H * H)
    for kind, (c, b, d, e) in HELICITY.items():
        converted = {
            "R_old": -(eta(a * a * eta(X[0])) + a * a * c * q * X[0]) / (2 * a),
            "R_old_squared": (
                -6 * (eta(r * eta(X[0])) + r * c * q * X[0]) + 2 * e * q * q * X[0]
            )
            / a,
            "Weyl_squared": (
                eta_power(X[0], 4) + b * q * eta_power(X[0], 2) + d * q * q * X[0]
            )
            / a,
        }
        for name, value in operators(kind).items():
            checks[f"full_proper_time_{kind}_{name}"] = s.cancel(
                converted[name] - value
            )
    checks["fourth_clock_derivative_complete"] = s.expand(
        eta_power(X[0], 4)
        - a**4
        * (
            X[4]
            + 6 * H * X[3]
            + (4 * H1 + 11 * H * H) * X[2]
            + (H2 + 7 * H * H1 + 6 * H**3) * X[1]
        )
    )
    checks["both_canonical_metric_directions"] = (
        4 / modes.KAPPA - 4 * s.Rational(1, 10) ** 800
    )
    return {
        "measure": "Because deta=dt/a and partial_eta=a partial_t, the proper-time density operator is a^-1 H_eta, not H_eta with a mere derivative substitution. This preserves the compact action bilinear and its self-adjoint pairing.",
        "helicity_constants": HELICITY,
        "full_proper_time_operators": {name: operators(name) for name in HELICITY},
        "clock": "The full fourth conformal derivative produces H,Hdot,Hddot and every friction term. No slow-clock or constant-a approximation is used.",
        "scope": "These helicities label the complete tracefree spatial subspace only. No lapse/shift constraint has been solved and no scalar or vector component is identified with a canonically reduced propagating degree of freedom.",
        "checks": checks,
        "gates": {
            "every_helicity_retained": set(HELICITY) == {"tensor", "vector", "scalar"},
            "scalar_gradient_sign_not_suppressed": HELICITY["scalar"][0] < 0,
            "zero_vector_Einstein_gradient_not_vector_dynamics_claim": HELICITY[
                "vector"
            ][0]
            == 0,
            "same_actual_kappa": modes.KAPPA == 10**800,
        },
    }
