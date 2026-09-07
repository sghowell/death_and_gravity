"""Literal -2 beta(phi_g,phi_f) edge, canonical clocks and undivided balance.

P8(b): +---, R_B=-6(DH+2H²), EH=-G R_B/2. The canonical clock stress
is not separately conserved: its scalar equation contains the edge force.
Other minimally coupled free fields remain separately conserved.
"""

from functools import cache

import sympy as sp


@cache
def derive():
    ni, nj, ai, aj, y, c = sp.symbols("N_g N_f a_g a_f y c", positive=True)
    pg, pf = sp.symbols("phi_g phi_f", real=True)
    hg, hf, qg, qf = sp.symbols("H_g H_f q_g q_f", real=True)
    betas = tuple(sp.Function(f"beta{i}")(pg, pf) for i in range(5))
    lag = -2*(ni*ai**3*betas[0]+betas[1]*(nj*ai**3+3*ni*ai**2*aj)
              +3*betas[2]*(nj*ai**2*aj+ni*ai*aj**2)
              +betas[3]*(3*nj*ai*aj**2+ni*aj**3)+betas[4]*nj*aj**3)
    point = {nj: c*ni, aj: y*ai}
    rho_g = sp.cancel((-sp.diff(lag, ni)/ai**3).subs(point))
    rho_f = sp.cancel((-sp.diff(lag, nj)/aj**3).subs(point))
    pressure_g = sp.cancel((sp.diff(lag, ai)/(3*ni*ai**2)).subs(point))
    pressure_f = sp.cancel((sp.diff(lag, aj)/(3*nj*aj**2)).subs(point))
    U = betas[0]+3*betas[1]*y+3*betas[2]*y**2+betas[3]*y**3
    V = betas[1]+3*betas[2]*y+3*betas[3]*y**2+betas[4]*y**3
    P = 2*(betas[1]+2*betas[2]*y+betas[3]*y**2)
    qprime_force_g = sp.cancel((sp.diff(lag, pg)/(ni*ai**3)).subs(point))
    qprime_force_f = sp.cancel((sp.diff(lag, pf)/(nj*aj**3)).subs(point))
    dy = y*(c*hf-hg)
    Cg = sp.diff(rho_g, y)*dy+sp.diff(rho_g, pg)*qg+sp.diff(rho_g, pf)*c*qf+3*hg*(rho_g+pressure_g)
    Cf = sp.diff(rho_f, y)*dy/c+sp.diff(rho_f, pg)*qg/c+sp.diff(rho_f, pf)*qf+3*hf*(rho_f+pressure_f)
    combined_g = sp.expand(Cg+qg*qprime_force_g)
    combined_f = sp.expand(Cf+qf*qprime_force_f)
    B = 3*P*(y*hf-hg)+2*(sp.diff(U, pf)*qf-sp.diff(V, pg)*qg)
    return {"N_g": ni, "N_f": nj, "a_g": ai, "a_f": aj, "y": y, "c": c,
            "phi_g": pg, "phi_f": pf, "H_g": hg, "H_f": hf, "q_g": qg, "q_f": qf,
            "betas": betas, "L_edge": lag, "U": U, "V": V, "P": P,
            "rho_g": rho_g, "rho_f": rho_f, "pressure_g": pressure_g, "pressure_f": pressure_f,
            "scalar_edge_force_g": qprime_force_g, "scalar_edge_force_f": qprime_force_f,
            "interaction_balance_g": Cg, "interaction_balance_f": Cf,
            "combined_balance_g": combined_g, "combined_balance_f": combined_f, "B": B}


@cache
def checks():
    d = derive()
    pg, pf, y, c, U, V, P = (d[key] for key in ("phi_g", "phi_f", "y", "c", "U", "V", "P"))
    return {
        "literal_g_lapse_density": sp.expand(d["rho_g"]-2*U),
        "literal_f_lapse_density": sp.expand(d["rho_f"]-2*V/y**3),
        "literal_g_null": sp.expand(d["rho_g"]+d["pressure_g"]-(y-c)*P),
        "literal_f_null": sp.expand(d["rho_f"]+d["pressure_f"]-(c-y)*P/(c*y**3)),
        "literal_g_scalar_force": sp.expand(d["scalar_edge_force_g"]+2*(sp.diff(U, pg)+c*sp.diff(V, pg))),
        "literal_f_scalar_force": sp.expand(d["scalar_edge_force_f"]+2*(sp.diff(U, pf)+c*sp.diff(V, pf))/(c*y**3)),
        "source_aware_g_Bianchi": sp.expand(d["combined_balance_g"]-c*d["B"]),
        "source_aware_f_Bianchi": sp.expand(d["combined_balance_f"]+d["B"]/(c*y**3)),
        "combined_two_lapse_reciprocity": sp.expand(d["combined_balance_g"]+c**2*y**3*d["combined_balance_f"]),
    }


def source_checks():
    """Canonical scalar balance and optional conformal-matter trace cancellation.

The certified fixture uses conformal factor one. These generic identities
only explain why a nontrivial physical Jordan metric cannot be ignored.
"""
    q, qp, h, potential_prime, force, alpha, trace = sp.symbols("q qprime H Wprime edge_force alpha trace_E", real=True)
    clock_balance = q*(qp+3*h*q+potential_prime)
    on_shell = {qp: -3*h*q-potential_prime+force-alpha*trace}
    return {"canonical_clock_is_sourced": sp.expand(clock_balance.subs(on_shell)-q*force+alpha*q*trace),
            "conformal_matter_exchange_cancels_locally": sp.expand(clock_balance.subs(on_shell)+alpha*q*trace-q*force)}
