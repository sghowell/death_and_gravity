"""Full source-aware coframe equations restricted AFTER variation to FLRW.

The matter density rho and pressure p are measured in the actual h=u^T eta u
frame, before multiplication by the positive source-strength epsilon. The
actual null stress is n_h=epsilon*(rho+p). All six lapses/scales are positive.
"""

from functools import cache

import sympy as sp
from p8_trimetric import matter, model


@cache
def derive():
    ne, ae, nv, av, nu, au = sp.symbols("n_e a_e n_v a_v n_u a_u", positive=True)
    pg, pf = sp.symbols("p_g p_f", positive=True)
    rho, pressure = sp.symbols("rho p", real=True)
    eps = sp.Symbol("epsilon", positive=True)
    b = sp.Symbol("B", real=True)
    pe, pv = pg*ae/au, pf*av/au
    ce, cv = au*ne/(nu*ae), au*nv/(nu*av)
    lapse = b+3*(pe+pv)+eps*rho/2
    space = b+pe*(ce+2)+pv*(cv+2)-eps*pressure/2
    return {"n_e": ne, "a_e": ae, "n_v": nv, "a_v": av, "n_u": nu, "a_u": au,
            "p_g": pg, "p_f": pf, "rho": rho, "p": pressure, "epsilon": eps, "B": b,
            "P_g": pe, "P_f": pv, "c_e": ce, "c_f": cv,
            "lapse_residual": lapse, "spatial_residual": space,
            "actual_null_stress": eps*(rho+pressure),
            "weighted_cone_residual": pe*(ce-1)+pv*(cv-1)-eps*(rho+pressure)/2}


@cache
def checks():
    d = derive()
    ne, ae, nv, av, nu, au = (d[key] for key in ("n_e", "a_e", "n_v", "a_v", "n_u", "a_u"))
    e, v, u = sp.diag(ne, ae, ae, ae), sp.diag(nv, av, av, av), sp.diag(nu, au, au, au)
    # Independent scalar gradient first; rewrite k0² only after differentiating.
    k0 = sp.Symbol("k0", real=True)
    source = matter.canonical_source(u, sp.Matrix([k0, 0, 0, 0]), (d["rho"]-d["p"])/2)
    eulers = model.euler_maps(e, v, u, b=d["B"], pg=d["p_g"], pf=d["p_f"],
                              matter_gradient=source["J_u"], epsilon=d["epsilon"])
    kinetic = {k0**2: nu**2*(d["rho"]+d["p"])}
    lapse = sp.expand(eulers["E_u"][0, 0]/(-2*au**3)).subs(kinetic)
    space = sp.expand(eulers["E_u"][1, 1]/(-2*nu*au**2)).subs(kinetic)
    q, n, a, capital_n, capital_a = sp.symbols("q n a N A", positive=True)
    symmetric = {ne: n, nv: n, ae: a, av: a, nu: capital_n*n, au: capital_a*a,
                 d["p_g"]: q, d["p_f"]: q}
    return {
        "full_covariant_u_lapse": sp.cancel(lapse-d["lapse_residual"]),
        "full_covariant_u_spatial": sp.cancel(space-d["spatial_residual"]),
        "undivided_weighted_null_cone_identity": sp.cancel(d["spatial_residual"]-d["lapse_residual"]-d["weighted_cone_residual"]),
        "symmetric_lapse": sp.cancel(d["lapse_residual"].subs(symmetric)-d["B"]-6*q/capital_a-d["epsilon"]*d["rho"]/2),
        "symmetric_space": sp.cancel(d["spatial_residual"].subs(symmetric)-d["B"]-2*q/capital_n-4*q/capital_a+d["epsilon"]*d["p"]/2),
        "symmetric_speed_gap": sp.cancel(d["weighted_cone_residual"].subs(symmetric)
                                         -2*q/capital_a*(capital_a/capital_n-1)+d["actual_null_stress"]/2),
    }


def controls():
    """Mixed-sign point is only an auxiliary solution; vacuum is full flat."""
    d = derive()
    algebraic = {d["n_e"]: sp.Rational(1, 2), d["n_v"]: 1, d["n_u"]: 1,
                 d["a_e"]: 1, d["a_v"]: 1, d["a_u"]: 1,
                 d["p_g"]: -2, d["p_f"]: 1, d["B"]: sp.Rational(5, 2),
                 d["rho"]: 1, d["p"]: 1, d["epsilon"]: 1}
    flat = model.euler_maps(sp.eye(4), sp.eye(4), sp.eye(4), b=3, pg=-2, pf=1, bg=2, bf=-1, epsilon=0)
    return {"mixed_sign_auxiliary_lapse": d["lapse_residual"].subs(algebraic),
            "mixed_sign_auxiliary_space": d["spatial_residual"].subs(algebraic),
            "mixed_sign_auxiliary_c_e": d["c_e"].subs(algebraic),
            "mixed_sign_auxiliary_c_f": d["c_f"].subs(algebraic),
            "mixed_sign_actual_null_stress": d["actual_null_stress"].subs(algebraic),
            "mixed_sign_full_flat_e": flat["E_e"][0, 0],
            "mixed_sign_full_flat_v": flat["E_v"][0, 0],
            "mixed_sign_full_flat_u": flat["E_u"][0, 0],
            "mixed_sign_full_flat_relative_spring": sp.Rational(2)*(-2)*1/(-2+1),
            "nonrolling_symmetric_gap": sp.Integer(0)}
