"""Scale restoration and finite-frequency chart diagnostics; NOT a cutoff."""

from functools import cache

import sympy as sp
from p8 import jets as j
from p8 import quadratic
from p8.matter import ia_completion
from p8.rational_candidates import X

from . import nonlinear_d as model


@cache
def scaled_background_checks():
    """Substitute actual dimensionful functions in the old covariant equations."""
    M, tau = sp.symbols("M tau", positive=True)
    f = model.functions()
    time_map = {model.u: j.t/tau}
    F = M**2/tau**2*f["F"].subs(time_map)
    A3 = M**2*f["a3"].subs(time_map)
    mapping = quadratic.covariant_N_jets(F, sp.Integer(0), -M**2/2, sp.Integer(0), A3, X)
    background = {j.H: model.H.subs(time_map)/tau, j.a: model.d.subs(time_map)}
    def evaluate(expr):
        return sp.factor(quadratic.substitute_functions(
            quadratic.substitute_functions(expr, mapping), background))
    q = quadratic.derive()
    out = {key: evaluate(quadratic.background()[key]) for key in ("EN", "Ea")}
    for key, target in (("GT", M**2), ("FT", M**2),
                        ("Theta", M**2/tau*f["Theta"].subs(time_map)),
                        ("Lambda", M**2*f["Lambda"].subs(time_map))):
        out[key] = sp.cancel(evaluate(q[key])-target)
    A2, A4, A5 = ia_completion(-M**2/2, 0, 0, A3, X)
    out.update({"A2": A2, "A4_homogeneity": sp.cancel(A4-M**2*f["a4"].subs(time_map)),
                "A5_homogeneity": sp.cancel(A5-M**2*f["a5"].subs(time_map))})
    sigma = evaluate(q["Sigma"])
    out["J"] = sp.cancel(sigma+3*(M**2/tau*f["Theta"].subs(time_map))**2/M**2
                          -M**2/tau**2*f["J"].subs(time_map))
    return out


@cache
def derive():
    M, tau = sp.symbols("M tau", positive=True)
    x = sp.Symbol("x", nonnegative=True)
    q = sp.Symbol("q", positive=True)
    f = model.functions()
    u, d = model.u, model.d
    # Curvature reference: E_ref=1/(tau sqrt(d)); never use |H| alone at H=0.
    inner = sp.cancel((4*x+2*(1-x))/(1+x))
    outer = sp.cancel((4*x+2*(x-1))/(1+x))
    kinetic_exact = q*f["J"]/(q*f["Lambda"]**2-f["J"])
    kinetic_principal = f["J"]/f["Lambda"]**2
    rel = sp.cancel(kinetic_exact/kinetic_principal-1)
    return {
        "action_prefactor": (M*tau)**2,
        "physical_functions": {"F": "M^2/tau^2 f(phi/tau,X)",
                               "F2": "-M^2/2", "Ai": "M^2 ai(phi/tau,X)",
                               "clock": "phi=t; X=1", "canonical_tail_clock": "2 M asinh(phi/tau)"},
        "E_curvature_squared_over_E_ref_squared": (inner, outer),
        "curvature_bounds_residuals": {"inner_equals_two": inner-2,
                                       "outer_minus_two": sp.factor(outer-2),
                                       "six_minus_outer": sp.factor(6-outer)},
        "gamma_kinetic_exact": kinetic_exact,
        "gamma_kinetic_relative_error": rel,
        "residuals": {"relative_error": sp.cancel(rel-f["J"]/(q*f["Lambda"]**2-f["J"])),
                      "physical_H_definition": sp.cancel(sp.diff(d, u)/(tau*d)-2*u/(tau*d)),
                      "A3_variation_squared_bound": sp.factor(36/d-(sp.diff(f["a3"], u)/f["a3"])**2)},
        "q": q,
    }
