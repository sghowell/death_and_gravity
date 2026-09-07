"""Literal auxiliary action and full source-aware vierbein Euler maps.

Vierbeine have Lorentz rows and coordinate columns. The matrix gradient
convention is delta L=tr(E_e.T*delta e), and analogously for v,u. All
matrices are real and invertible on a fixed positive-orientation chart;
the parent square-root/symmetrization domain is a separate requirement.
Einstein tensors supplied to euler_maps have two upper coordinate indices
in the repository +--- curvature convention. No EH term exists for u.
"""

from functools import cache

import sympy as sp

ETA = sp.diag(1, -1, -1, -1)
B, PG, PF, BG, BF, EPS = sp.symbols("B p_g p_f beta4g beta4f epsilon", real=True)
G, F = sp.symbols("G F", positive=True)


def matrix(value):
    out = sp.Matrix(value)
    if out.shape != (4, 4):
        raise ValueError("A four-dimensional vierbein or matrix is required")
    return out


def potential_density(e, v, u, *, b=B, pg=PG, pf=PF, bg=0, bf=0):
    """Interaction Lagrangian density, including optional separate beta4s."""
    e, v, u = map(matrix, (e, v, u))
    return -2*u.det()*(b+sp.trace(u.inv()*(pg*e+pf*v)))-2*bg*e.det()-2*bf*v.det()


def euler_maps(e, v, u, *, matter_gradient=None, einstein_g=None, einstein_f=None,
               b=B, pg=PG, pf=PF, bg=0, bf=0, epsilon=EPS, g=G, f=F):
    """Full covariant density gradients, with literal matter gradient J_u.

    EH contributions vanish when the supplied Einstein tensors vanish;
    default zeros select the algebraic/constant-flat audit, not a claim
    that the full gravitational theory has no derivative terms.
    matter_gradient is delta L_m/delta u, without the epsilon prefactor.
    The scalar Euler functional remains epsilon*delta S_m/delta psi.
    """
    e, v, u = map(matrix, (e, v, u))
    ju = sp.zeros(4) if matter_gradient is None else matrix(matter_gradient)
    eg = sp.zeros(4) if einstein_g is None else matrix(einstein_g)
    ef = sp.zeros(4) if einstein_f is None else matrix(einstein_f)
    inv, det, q = u.inv(), u.det(), pg*e+pf*v
    scalar = b+sp.trace(inv*q)
    ee = g*e.det()*ETA*e*eg-2*pg*det*inv.T-2*bg*e.det()*e.inv().T
    ev = f*v.det()*ETA*v*ef-2*pf*det*inv.T-2*bf*v.det()*v.inv().T
    eu = -2*det*(scalar*inv.T-inv.T*q.T*inv.T)+epsilon*ju
    return {"E_e": ee, "E_v": ev, "E_u": eu, "Q": q,
            "matter_gradient": ju, "det_u": det, "u_inverse": inv}


def trace_reduced_source(u, matter_gradient):
    """Source-paper matrix T and its trace reversal, derived from J_u.

    This T is a vierbein-source convention, not an unqualified physical
    stress tensor. Inverse-vierbein variation gives T=u*J_u.T*u/(2det u).
    """
    u, ju = map(matrix, (u, matter_gradient))
    stress = u*ju.T*u/(2*u.det())
    trace = sp.trace(u.inv()*stress)
    return {"T": stress, "trace_T": trace, "calT": trace*u/3-stress}


def flat_link_contractions(e, v, u, *, pg=PG, pf=PF):
    """Undivided nonzero-link obstruction, no use of B or the u equation."""
    out = euler_maps(e, v, u, pg=pg, pf=pf, epsilon=0)
    u = matrix(u)
    return {"uT_E_e": u.T*out["E_e"], "uT_E_v": u.T*out["E_v"]}


@cache
def checks():
    es, vs, us = (sp.symbols(f"{name}0:4", positive=True) for name in ("e", "v", "u"))
    e, v, u = (sp.diag(*values) for values in (es, vs, us))
    action = potential_density(e, v, u, bg=BG, bf=BF)
    out = euler_maps(e, v, u, bg=BG, bf=BF, epsilon=0)
    result = {}
    for name, variables, gradient in (("e", es, out["E_e"]), ("v", vs, out["E_v"]), ("u", us, out["E_u"])):
        for index, variable in enumerate(variables):
            result[f"literal_diagonal_{name}_Euler_{index}"] = sp.cancel(sp.diff(action, variable)-gradient[index, index])
    ju = sp.diag(*sp.symbols("J0:4", real=True))
    source = trace_reduced_source(u, ju)
    out = euler_maps(e, v, u, matter_gradient=ju)
    c = out["Q"]+B*u/3-EPS*source["calT"]
    expected = 2*u.det()*u.inv().T*(c.T-sp.trace(u.inv()*c)*u.T)*u.inv().T
    for index in range(4):
        result[f"undivided_trace_reduced_source_Euler_{index}"] = sp.cancel(out["E_u"][index, index]-expected[index, index])
    flat = flat_link_contractions(e, v, u)
    result["flat_g_link_undivided_contraction"] = sp.cancel(flat["uT_E_e"][0, 0]+2*PG*u.det())
    result["flat_f_link_undivided_contraction"] = sp.cancel(flat["uT_E_v"][0, 0]+2*PF*u.det())
    return result
