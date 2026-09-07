"""Literal metric and clock variations before the proportional restriction.

P8(b): +---, EH=-M² R_B/2, interaction=-2 sum beta_n e_n. The beta_n
have mass dimension four; the constant canonical clock has dimension one.
The positive vacuum root r is not the lapse parameter of the bounce action.
"""

from functools import cache

import sympy as sp


@cache
def derive():
    ag, af, ng, nf, r, t = sp.symbols("a_g a_f N_g N_f r t", positive=True)
    phi = sp.Symbol("phi", real=True)
    betas = tuple(sp.Function(f"beta{j}")(phi) for j in range(5))
    generating = sp.Poly((1+t*nf/ng)*(1+t*af/ag)**3, t)
    density = sp.expand(-2*ng*ag**3*sum(betas[j]*generating.nth(j) for j in range(5)))
    at = {af: r*ag, nf: r*ng}
    rho_g = sp.factor((-sp.diff(density, ng)/ag**3).subs(at))
    rho_f = sp.factor((-sp.diff(density, nf)/af**3).subs(at))
    pressure_g = sp.factor((sp.diff(density, ag)/(3*ng*ag**2)).subs(at))
    pressure_f = sp.factor((sp.diff(density, af)/(3*nf*af**2)).subs(at))
    clock = sp.factor((-sp.diff(density, phi)/(ng*ag**3)).subs(at))
    U = betas[0]+3*betas[1]*r+3*betas[2]*r**2+betas[3]*r**3
    V = betas[1]+3*betas[2]*r+3*betas[3]*r**2+betas[4]*r**3
    W = sum(sp.binomial(4, j)*betas[j]*r**j for j in range(5))
    return {"a_g": ag, "a_f": af, "N_g": ng, "N_f": nf, "r": r, "phi": phi,
            "betas": betas, "density": density, "rho_g": rho_g, "rho_f": rho_f,
            "pressure_g": pressure_g, "pressure_f": pressure_f,
            "clock_equation": clock, "U": U, "V": V, "W": W}


def checks():
    d = derive()
    r, phi, U, V = (d[key] for key in ("r", "phi", "U", "V"))
    b0, b1, b4 = sp.symbols("b0 b1 b4", real=True)
    candidate_ratio = -b0/(3*b1)
    return {
        "literal_g_density": sp.factor(d["rho_g"]-2*U),
        "literal_f_density": sp.factor(d["rho_f"]-2*V/r**3),
        "literal_g_spatial_equation": sp.factor(d["pressure_g"]+d["rho_g"]),
        "literal_f_spatial_equation": sp.factor(d["pressure_f"]+d["rho_f"]),
        "literal_clock_equation": sp.factor(d["clock_equation"]-2*sp.diff(d["W"], phi)),
        "potential_metric_decomposition": sp.expand(d["W"]-U-r*V),
        "radial_potential_variation": sp.expand(sp.diff(d["W"], r)-4*V),
        "metric_eliminant_identity_for_nonzero_b1": sp.factor(
            27*b1**3*(b1+b4*candidate_ratio**3)-(27*b1**4-b0**3*b4)),
    }


def controls():
    """Coefficient-point controls, not claims about the specified action."""
    d = derive()
    r, phi, betas = d["r"], d["phi"], d["betas"]
    values = {r: 2, betas[0]: -6, betas[1]: 1, betas[2]: 0,
              betas[3]: 0, betas[4]: -sp.Rational(1, 8)}
    derivatives = {sp.diff(beta, phi): 0 for beta in betas}
    generic = {name: sp.factor(d[name].subs(derivatives).subs(values))
               for name in ("rho_g", "rho_f", "clock_equation")}
    if set(generic.values()) != {0}:
        raise ValueError("The separately named constant-coefficient vacuum control failed")
    sourced = dict(derivatives)
    sourced[sp.diff(betas[0], phi)] = 1
    omitted_clock = sp.factor(d["clock_equation"].subs(sourced).subs(values))
    if omitted_clock != 2:
        raise ValueError("Metric equations alone did not detect the omitted clock equation")
    return {"separate_constant_coefficient_point_all_equations": generic,
            "same_metric_values_with_nonzero_clock_slope": omitted_clock,
            "scope": "Controls of the equations, not the frozen local action or a UV/health claim"}
