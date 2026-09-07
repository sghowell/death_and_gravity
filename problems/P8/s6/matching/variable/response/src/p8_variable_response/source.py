"""A conserved physical-g TT probe and an explicit g-only loading operator.

These are linear external sources. Source amplitude is differentiated at zero
before taking delta limits. No nonlinear stress realization, energy condition,
low temporal band or finite-amplitude backreaction is inferred.
"""

from functools import cache

import sympy as sp

from . import operator


def coefficients():
    d = operator.derive()
    u, c, den = d["u"], d["c"], d["d"]
    return {"u": u, "c": c,
            "j_light": sp.sqrt(c)*den**9/sp.sqrt(c*den**12+8),
            "j_heavy": -2*sp.sqrt(2)*den**3/sp.sqrt(c*den**12+8),
            "physical_J_over_sigma": den**6/2,
            "sigma_dictionary": "sigma=tau^2*Pi/M^2; S_probe=1/2 integral dT a^3 Pi gamma_g"}


def loading(gamma_f):
    """Exact g-only loading on a punctured interval with positive spring.

    For the construction choose gamma_f=zeta*gamma_f,target, where the target
    solves BOTH homogeneous delta=0 equations on [-2/100,-1/100] and zeta is
    smooth, flat zero near the left end and flat one near the right end.
    The returned source then has compact support and gives the target data.
    This routine supplies the differential operator, not a sampled target.
    """
    d = operator.derive()
    u, c, a, b = d["u"], d["c"], d["a"], d["b"]
    kg, kf = d["K1"], d["K2"]
    gg, gf = a/8, c*b/8
    # U_TT from the literal unit-polarization action, not m^2 itself.
    spring = d["mass"]*d["K_relative"]
    k = d["kbar"]
    gamma_g = gamma_f+(sp.diff(kf*sp.diff(gamma_f, u), u)+gf*k**2*gamma_f)/spring
    sigma = 4/a**3*(sp.diff(kg*sp.diff(gamma_g, u), u)+gg*k**2*gamma_g
                    +spring*(gamma_g-gamma_f))
    return {"gamma_g": gamma_g, "sigma": sigma, "spring": spring,
            "K_g": kg, "K_f": kf, "G_g": gg, "G_f": gf}


def g_observability():
    """Source-off physical data to (gamma_g,gamma_g,u,gamma_g,uu,gamma_g,uuu).

    This finite-jet observable is not a low temporal-band observable. At a
    fixed punctured slice its determinant is (U_TT/K_g)^2>0, including the
    auxiliary delta=0 operator; a claim about only (g,g') would be weaker.
    """
    d = operator.derive()
    u = d["u"]
    spring = d["mass"]*d["K_relative"]
    j = spring/d["K1"]
    aa = -(d["a"]*d["kbar"]**2/8+spring)/d["K1"]
    bb = -sp.diff(d["K1"], u)/d["K1"]
    matrix = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [aa, bb, j, 0],
                        [sp.diff(aa, u)+bb*aa, aa+sp.diff(bb, u)+bb**2,
                         sp.diff(j, u)+bb*j, j]])
    return {"matrix": matrix, "relative_coefficient": sp.factor(j),
            "determinant": sp.factor(j**2),
            "proper_jet_clock": "multiply rows by (1,tau^-1,tau^-2,tau^-3)"}


@cache
def checks():
    d, s = operator.derive(), coefficients()
    u = d["u"]
    out = {
        "source_light_physical_projection": sp.simplify(s["j_light"]-d["a"]**3/(2*d["f_sum"])),
        "source_heavy_physical_projection": sp.simplify(s["j_heavy"]+d["a"]**3*d["w2"]/(2*d["f_relative"])),
        "center_light": sp.simplify(s["j_light"].subs({u: 0, d["c"]: 2})-1/sp.sqrt(5)),
        "center_heavy": sp.simplify(s["j_heavy"].subs({u: 0, d["c"]: 2})+2/sp.sqrt(5)),
    }
    # Generic independent Euler identities do not divide by a TT solution.
    kg, kf, gg, gf, spring, a = (sp.Function(name)(u) for name in ("Kg", "Kf", "Gg", "Gf", "Ut", "a"))
    field = sp.Function("f")(u)
    k = d["kbar"]
    g = field+(sp.diff(kf*sp.diff(field, u), u)+gf*k**2*field)/spring
    sigma = 4/a**3*(sp.diff(kg*sp.diff(g, u), u)+gg*k**2*g+spring*(g-field))
    ef = sp.diff(kf*sp.diff(field, u), u)+gf*k**2*field+spring*(field-g)
    eg = 2*(sp.diff(kg*sp.diff(g, u), u)+gg*k**2*g+spring*(g-field))-a**3*sigma/2
    out["actual_unsourced_f_equation"] = sp.expand(ef)
    out["actual_g_source_normalization"] = sp.expand(eg)
    # The exact physical source is transverse with no 0 components. Its
    # energy divergence is H times its spatial trace, and its spatial
    # divergence is k_j e_ij=0 for k along z and diag(1,-1,0)/sqrt(2).
    trace, divergence = sp.trace(sp.diag(1, -1, 0)), sp.diag(1, -1, 0)*sp.Matrix([0, 0, k])
    out["TT_energy_balance_trace"] = trace
    out["TT_spatial_balance"] = sum(v**2 for v in divergence)
    obs = g_observability()
    out["g_only_source_off_jet_observability"] = sp.factor(obs["matrix"].det()-obs["determinant"])
    return out
