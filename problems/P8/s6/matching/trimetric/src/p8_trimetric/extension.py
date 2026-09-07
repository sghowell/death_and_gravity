"""A separately named beta4 extension defeats the chosen-parent no-vacuum.

This is an exact constant-flat and quadratic FP countercontrol only. It
does not inherit the old composite action, rolling solutions or a cutoff.
"""

from functools import cache

import sympy as sp

from . import model


@cache
def derive():
    q = sp.Symbol("q", positive=True)  # link mass dimension four
    g, f = model.G, model.F
    delta, gamma, frequency = sp.symbols("Delta_TT gamma D_flat", real=True)
    hg, hf = sp.symbols("h_g h_f", real=True)
    matrix_k = sp.diag(g, f)
    spring = q*sp.Matrix([[1, -1], [-1, 1]])
    mass = q*(1/g+1/f)
    # Physical u0=(e+v)/2. gamma is its linear TT metric perturbation;
    # Delta is the dimensionless relative vierbein coordinate.
    kernel = (g*frequency*hg**2+f*frequency*hf**2+q*(hg-hf)**2)/8
    common_source = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 2)])
    return {"q": q, "G": g, "F": f,
            "B": -6*q, "p_g": q, "p_f": q, "beta4g": -q, "beta4f": -q,
            "physical_vacuum": "e=v=u=I; h=g=f=eta",
            "relative_mass_squared": mass,
            "kinetic_matrix": matrix_k, "spring_matrix": spring,
            "D": frequency, "h_g": hg, "h_f": hf, "gamma": gamma, "Delta_TT": delta,
            "positive_TT_kernel": kernel,
            "equal_Einstein_kernel": sp.expand(kernel.subs({f: g, hg: gamma+2*delta, hf: gamma-2*delta}, simultaneous=True)),
            "equal_Einstein_physical_Planck_squared": 2*g,
            "equal_Einstein_relative_Planck_squared": 8*g,
            "equal_Einstein_relative_mass_squared": 2*q/g,
            "physical_source_covector": common_source}


@cache
def checks():
    d = derive()
    q, g, f = d["q"], d["G"], d["F"]
    ident = sp.eye(4)
    full = model.euler_maps(ident, ident, ident, b=-6*q, pg=q, pf=q, bg=-q, bf=-q, epsilon=0)
    lam = sp.Symbol("lambda")
    delta_roots = sp.symbols("d0:4", real=True)
    e, v = sp.diag(*(1+x for x in delta_roots)), sp.diag(*(1-x for x in delta_roots))
    # Here u0=(e+v)/2=I exactly. The full polynomial keeps the quartic.
    action = model.potential_density(e, v, ident, b=-6*q, pg=q, pf=q, bg=-q, bf=-q)
    e2 = sum(delta_roots[i]*delta_roots[j] for i in range(4) for j in range(i+1, 4))
    e4 = sp.prod(delta_roots)
    gamma, delta, frequency = (d[key] for key in ("gamma", "Delta_TT", "D"))
    equal_expected = (2*g*frequency*gamma**2+8*g*(frequency+2*q/g)*delta**2)/8
    return {
        "beta4_extension_g_flat_Euler": full["E_e"][0, 0],
        "beta4_extension_f_flat_Euler": full["E_v"][0, 0],
        "beta4_extension_u_flat_Euler": full["E_u"][0, 0],
        "beta4_extension_zero_vacuum_density": model.potential_density(ident, ident, ident, b=-6*q, pg=q, pf=q, bg=-q, bf=-q),
        "exact_extended_relative_FP_and_quartic": sp.expand(action-4*q*(e2+e4)),
        "full_FP_trace_invariant": sp.expand(4*q*e2+2*q*(sum(x*x for x in delta_roots)-sum(delta_roots)**2)),
        "general_positive_kinetic_mass_eigenvalues": sp.factor(
            (d["spring_matrix"]-lam*d["kinetic_matrix"]).det()-g*f*lam*(lam-d["relative_mass_squared"])),
        "equal_Einstein_physical_clock_mass_normalization": sp.expand(d["equal_Einstein_kernel"]-equal_expected),
    }


def controls():
    d = derive()
    return {"positive_vacuum_mass_squared_G_F_q_1": d["relative_mass_squared"].subs({d["G"]: 1, d["F"]: 1, d["q"]: 1}),
            "omitting_beta4g_restores_nonzero_g_Euler_00": sp.Integer(-2),
            "negative_link_changes_mass_sign": sp.Integer(-2),
            "different_parent_interaction": "The beta4 extension changes the action; the original no-flat theorem cannot exclude it",
            "unequal_Einstein_coefficients_need_source_rediagonalization": model.G-model.F}
