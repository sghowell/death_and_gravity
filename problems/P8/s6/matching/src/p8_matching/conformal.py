"""Independent conformal geometry and exact near-constant Planck-jet tests."""

from functools import cache

import sympy as sp

from . import background as bg


def exact_rational(value):
    value = sp.sympify(value)
    if value.is_Rational is not True:
        raise ValueError("Use finite exact rational errors")
    return value


def error_bound(epsilon0, epsilon1, epsilon2, *, geometry_error=0):
    """Necessary NEC inequality only; a nonpositive lower bound is no existence proof."""
    e0, e1, e2, eta = map(exact_rational, (epsilon0, epsilon1, epsilon2, geometry_error))
    if not (0 <= e0 < 1 and e1 >= 0 and e2 >= 0 and 0 <= eta < 2):
        raise ValueError("Require 0<=epsilon0<1, nonnegative jet errors and 0<=eta<2")
    lower = sp.factor((4-2*eta)*(1-e0)-e2-3*e1**2/(2*(1-e0)))
    status = "EXCLUDED_BY_CONFORMAL_NEC" if lower > 0 else "NEC_LOWER_BOUND_INCONCLUSIVE"
    return {"lower_margin": str(lower), "status": status,
            "scope": "H(0)=0; supplied physical-time Planck jets and H_dot error only"}


@cache
def identities():
    t, a = bg.t, bg.a
    q = sp.Function("q")(t)
    M, tau = sp.symbols("M tau", positive=True)
    H = sp.diff(a, t)/a
    qd, qdd = sp.diff(q, t), sp.diff(q, t, 2)
    aE, NE = sp.sqrt(q)*a, sp.sqrt(q)
    HE = sp.simplify(sp.diff(aE, t)/(aE*NE))
    HdE = sp.simplify(sp.diff(HE, t)/NE)
    expected = (sp.diff(H, t)+qdd/(2*q)-3*qd**2/(4*q**2)-H*qd/(2*q))/q
    # Transform the positive Einstein-frame field metric AND its proper-time velocities.
    projected_einstein_kinetic = bg.kinetic/q**2+3*M**2*qd**2/(2*q**3)
    parent_kinetic_solution = -2*M**2*q*sp.diff(H, t)-M**2*qdd+M**2*H*qd
    bridge = sp.simplify((-2*M**2*HdE-projected_einstein_kinetic).subs(bg.kinetic, parent_kinetic_solution))
    q0 = sp.Symbol("q0", positive=True)
    q1, q2 = sp.symbols("q1 q2", real=True)
    # Differentiate first, then evaluate the physical bounce geometry and q jets.
    at_bounce = {a: 1, sp.diff(a, t): 0, sp.diff(a, t, 2): 2/tau**2,
                 q: q0, qd: q1/tau, qdd: q2/tau**2}
    bounce = sp.factor(HdE.subs(at_bounce)*tau**2)
    numerator = 4*q0+q2-3*q1**2/(2*q0)
    eps, eta = sp.symbols("epsilon eta", real=True)
    lower = 4*(1-eps)-eps-3*eps**2/(2*(1-eps))
    factored = (2-eps)*(4-7*eps)/(2*(1-eps))
    vcap, wcap = sp.symbols("epsilon1 epsilon2", nonnegative=True)
    general_margin = (4-2*eta)*q0-wcap-3*vcap**2/(2*q0)
    residuals = {"H_E": sp.simplify(HE-(H+qd/(2*q))/sp.sqrt(q)),
                 "H_E_dot": sp.simplify(HdE-expected), "Einstein_kinetic_bridge": bridge,
                 "bounce_numerator": sp.simplify(2*q0**2*bounce-numerator),
                 "common_error_factorization": sp.factor(lower-factored),
                 "margin_monotonicity": sp.diff(general_margin, q0)-(4-2*eta+3*vcap**2/(2*q0**2)),
                 "threshold_zero": sp.factor(lower.subs(eps, sp.Rational(4, 7)))}
    return {"H_E": HE, "H_E_dot": HdE, "Einstein_time_kinetic": projected_einstein_kinetic,
            "bounce_H_E_dot_times_tau_squared": bounce,
            "necessary_bounce_numerator": numerator, "common_error_lower_margin": sp.factor(lower),
            "residuals": residuals}


def derive():
    result = identities()
    if any(sp.simplify(value) != 0 for value in result["residuals"].values()):
        raise ValueError("Conformal parent identity failed")
    return {**{key: str(value) for key, value in result.items() if key != "residuals"},
            "common_C2_error_necessary_lower_bound": "4/7",
            "general_jet_lower_bound": "(4-2*eta)*(1-epsilon0)-epsilon2-3*epsilon1^2/(2*(1-epsilon0))",
            "examples_not_parent_constructions": {
                "ten_percent_C2": error_bound(sp.Rational(1, 10), sp.Rational(1, 10), sp.Rational(1, 10)),
                "threshold_not_sufficient": error_bound(*([sp.Rational(4, 7)]*3)),
                "small_amplitude_fast_variation": error_bound(sp.Rational(1, 100), 2, sp.Rational(1, 100)),
                "ten_percent_geometry_error": error_bound(*([sp.Rational(1, 10)]*3), geometry_error=sp.Rational(1, 10))},
            "exact_residuals": {key: "0" for key in result["residuals"]},
            "scope": "Classical nonminimal two-derivative scalar parent with healthy Einstein field metric; not all DHOST or UV theories"}
