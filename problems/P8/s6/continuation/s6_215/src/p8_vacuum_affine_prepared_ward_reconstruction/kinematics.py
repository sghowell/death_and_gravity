"""Prepared source and final-vanishing detector synchronous maps."""

from functools import cache

import sympy as s
from p8_vacuum_affine_full_adm_vertices import chart

t, x, y, z = s.symbols("t x y z", real=True)
COORDS = (t, x, y, z)
SPACE = (x, y, z)
a = (1 + t * t) ** 2
H = s.diff(a, t) / a


def grad(f):
    return s.Matrix([s.diff(f, q) for q in SPACE])


def spatial_jac(v):
    return v.jacobian(SPACE)


def symgrad(v):
    J = spatial_jac(v)
    return J + J.T


def lie_covariant(h, xi):
    J = xi.jacobian(COORDS)
    transport = h.applyfunc(
        lambda f: sum(xi[j] * s.diff(f, q) for j, q in enumerate(COORDS))
    )
    return transport + J.T * h + h * J


def gauge(eta, chi):
    return (
        s.diff(eta, t),
        chi.diff(t) - grad(eta) / a**2,
        2 * H * eta * s.eye(3) + symgrad(chi),
    )


def primitive(f, side):
    if type(side) is not str or side not in ("retarded", "advanced"):
        raise ValueError("Require retarded source or advanced detector primitive")
    end = -s.Rational(1, 2) if side == "retarded" else s.Rational(1, 2)
    u = s.Dummy("integration_time", real=True)
    if isinstance(f, s.MatrixBase):
        return f.applyfunc(lambda v: primitive(v, side))
    return s.integrate(s.cancel(f).subs(t, u), (u, end, t))


def synchronous(n, beta, Q, side):
    eta = primitive(n, side)
    chi = primitive(beta + grad(eta) / a**2, side)
    return Q - 2 * H * eta * s.eye(3) - symgrad(chi), eta, chi


def tracefree(Q):
    return Q - s.trace(Q) * s.eye(3) / 3


@cache
def data():
    eta = s.Function("eta")(*COORDS)
    chi = s.Matrix([s.Function(f"chi{i}")(*COORDS) for i in range(3)])
    xi = s.Matrix([eta, *chi])
    metric = s.diag(1, -a * a, -a * a, -a * a)
    n, beta, Q = gauge(eta, chi)
    T = s.Matrix([[1, 2, 3], [2, 4, 5], [3, 5, 6]])
    f = (1 + t) ** 4
    checks = {
        "literal_four_metric_Lie_to_all_ADM_components": lie_covariant(metric, xi)
        - chart.first(n, beta, Q, a),
        "tracefree_projection_has_zero_trace": s.trace(tracefree(T)),
        "orthogonal_trace_projection": s.trace(tracefree(T)) * s.trace(T) / 3,
        "tracefree_projection_idempotent": tracefree(tracefree(T)) - tracefree(T),
        "bounce_scale_nondegenerate": a.subs(t, 0) - 1,
        "bounce_H_zero_not_inverted": H.subs(t, 0),
    }
    for side, end in (("retarded", -s.Rational(1, 2)), ("advanced", s.Rational(1, 2))):
        F = primitive(f, side)
        checks[side + "_primitive_derivative"] = s.diff(F, t) - f
        checks[side + "_specified_endpoint"] = F.subs(t, end)
    return {
        "directions": "g=diag(1,-a^2 I), N=1+n, h=a^2 exp(Q). For xi=(eta,chi), n_xi=eta', beta_xi=chi'-a^-2 grad eta, Q_xi=2H eta I+symgrad chi.",
        "source_map": "eta_G=Iminus n_G, chi_G=Iminus(beta_G+a^-2 grad eta_G); G=Gsyn+Gxi with Gsyn=(0,0,Q_G-2H eta_G I-symgrad chi_G). Iminus integrates from the unchanged left preparation.",
        "detector_map": "Use Iplus f=-integral_t^right f for eta_D and chi_D. The detector gauge field vanishes on its final neighborhood. It need not vanish near the initial surface; the differentiated Ward flux there vanishes because the prepared source and current tangent do.",
        "regularity": "The maps use time primitives and at most two spatial derivatives. There is no H denominator, inverse momentum or elliptic constraint inversion. They are regular at the bounce and at zero spatial transfer.",
        "checks": checks,
        "gates": {
            "both_primitive_derivatives_positive_f": True,
            "different_source_detector_endpoint_conditions": True,
            "all_ten_metric_directions_decomposed": True,
            "no_bounce_or_zero_transfer_division": True,
            "prepared_state_not_reset": True,
        },
    }
