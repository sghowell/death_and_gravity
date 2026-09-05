"""Exact nonminimal background control OUTSIDE the D matching/tail contract."""

from functools import cache

import sympy as sp

from . import background as bg


@cache
def derive():
    u = sp.Symbol("u", real=True)
    M, tau = sp.symbols("M tau", positive=True)
    d = 1+u**2
    q = sp.exp(-3*u**2)
    physical_H = 2*u/(tau*d)
    HE = sp.factor((physical_H+sp.diff(q, u)/(2*tau*q))/sp.sqrt(q))
    HEd = sp.factor(sp.diff(HE, u)/(tau*sp.sqrt(q)))
    positive = 1+11*u**2+15*u**4+9*u**6
    # phi=t. ZE is the Einstein-frame field metric in the phi coordinate.
    ZE = 2*M**2*positive/(tau**2*d**2)
    VE = sp.factor(M**2*(3*HE**2+HEd))
    Q = M**2*q
    Qdot = sp.diff(Q, u)/tau
    GJ = sp.factor(q*ZE-3*M**2*sp.diff(q, u)**2/(2*tau**2*q))
    VJ = sp.factor(q**2*VE)
    field_velocity_E = 1/sp.sqrt(q)
    kinetic_E = ZE*field_velocity_E**2
    scalar_E = (sp.diff(ZE*field_velocity_E, u)/(tau*sp.sqrt(q))
                +3*HE*ZE*field_velocity_E-sp.diff(ZE, u)*field_velocity_E**2/(2*tau)
                +sp.diff(VE, u)/tau)
    # Replay independently derived physical-frame lapse/scale equations.
    at = {bg.a: d, sp.diff(bg.a, bg.t): sp.diff(d, u)/tau,
          sp.diff(bg.a, bg.t, 2): sp.diff(d, u, 2)/tau**2,
          bg.Q: Q, sp.diff(bg.Q, bg.t): Qdot,
          sp.diff(bg.Q, bg.t, 2): sp.diff(Q, u, 2)/tau**2,
          bg.kinetic: GJ, bg.potential: VJ}
    residuals = {"positive_H_E_decay": sp.simplify(-HEd-sp.exp(3*u**2)*positive/(tau**2*d**2)),
                 "Einstein_kinetic": sp.simplify(kinetic_E+2*M**2*HEd),
                 "Einstein_Friedmann": sp.simplify(3*M**2*HE**2-kinetic_E/2-VE),
                 "Einstein_scalar_equation": sp.simplify(scalar_E),
                 "conformal_field_metric": sp.simplify(ZE-GJ/q-3*M**2*sp.diff(q, u)**2/(2*tau**2*q**2))}
    for key in ("lapse_equation", "scale_equation"):
        residuals[f"Jordan_{key}"] = sp.simplify(bg.equations()[key].subs(at))
    if any(value != 0 for value in residuals.values()):
        raise ValueError(f"Nonminimal escape-control residuals: {residuals}")
    return {"physical_a": str(d), "q": str(q), "Einstein_H": str(HE),
            "Einstein_H_dot": str(HEd), "Einstein_field_metric": str(ZE),
            "Jordan_Q": str(Q), "Jordan_field_metric": str(GJ), "Jordan_potential": str(VJ),
            "Einstein_kinetic_positive_polynomial_coefficients_in_u_squared": [1, 11, 15, 9],
            "bounce_Planck_jets": [str(sp.diff(q, u, n).subs(u, 0)) for n in range(3)],
            "past_q_limit": str(sp.limit(q, u, -sp.oo)), "future_q_limit": str(sp.limit(q, u, sp.oo)),
            "exact_residuals": {key: "0" for key in residuals},
            "warning": "Background-only nonminimal escape control. Q tends to zero, violates P8 tensor tails; not the D extension, no vacuum/UV/perturbation verdict"}
