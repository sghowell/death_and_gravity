"""Conditional same-action flat calibration and mixed positive-FP branch.

The global background theorem does NOT require a vacuum. A positive
algebraic cap or a calibrated Einstein equation does not prove the full
u/scalar vacuum equations. The mixed-sign corollary assumes a real full
vacuum with a finite positive FP mass, then proves the TT inverse stays away
from zero on all regular same-action flat-FLRW points.
"""

from functools import cache

import sympy as sp
from p8_trimetric import model


@cache
def derive():
    rg, rf, g, f = sp.symbols("R_g0 R_f0 G F", positive=True)
    pg, pf = sp.symbols("p_g p_f", real=True)
    pe, pv = pg/rg, pf/rf
    total = pe+pv
    kg, kf = g/rg**2, f/rf**2
    spring = 2*pe*pv/total
    return {"R_g0": rg, "R_f0": rf, "G": g, "F": f, "p_g": pg, "p_f": pf,
            "P_g0": pe, "P_f0": pv, "S0": total,
            "b_g": -pg*rg**3, "b_f": -pf*rf**3, "B_effective": -3*total,
            "K_g0": kg, "K_f0": kf, "q_eff0": spring,
            "m_FP0_squared": spring*(1/kg+1/kf),
            "quadratic_FP_invariant_coefficient": pe*pv/(4*total)}


@cache
def checks():
    d = derive()
    rg, rf, g, f, pg, pf = (d[key] for key in ("R_g0", "R_f0", "G", "F", "p_g", "p_f"))
    u = sp.eye(4)
    e, v = u/rg, u/rf
    full = model.euler_maps(e, v, u, b=d["B_effective"], pg=pg, pf=pf,
                            bg=d["b_g"], bf=d["b_f"], epsilon=0)
    result = {f"calibrated_flat_{key}": sp.cancel(full[key][0, 0]) for key in ("E_e", "E_v", "E_u")}
    result["calibrated_total_vacuum_density"] = sp.cancel(model.potential_density(e, v, u, b=d["B_effective"],
                                       pg=pg, pf=pf, bg=d["b_g"], bf=d["b_f"]))
    # Literal full determinant second variation, with general traceful jets.
    # A diagonal test with independent eigenvalues certifies the symmetric
    # quadratic trace polynomial, not only its TT restriction.
    hs, js = sp.symbols("h0:4", real=True), sp.symbols("j0:4", real=True)
    pe, pv = d["P_g0"], d["P_f0"]
    total = d["S0"]
    def e2(values):
        return sum(values[i]*values[j] for i in range(4) for j in range(i+1, 4))
    det_coefficient = (pe*e2(hs)+pv*e2(js)-e2([pe*h+pv*j for h, j in zip(hs, js, strict=True)])/total)/2
    difference = [h-j for h, j in zip(hs, js, strict=True)]
    fp = pe*pv/(4*total)*(sum(difference)**2-sum(x*x for x in difference))
    result["literal_full_FP_trace_invariant"] = sp.cancel(det_coefficient-fp)
    eigenvalue = sp.Symbol("lambda", real=True)
    kinetic = sp.diag(d["K_g0"], d["K_f0"])
    spring = d["q_eff0"]*sp.Matrix([[1, -1], [-1, 1]])
    result["actual_physical_vacuum_mass_normalization"] = sp.factor((spring-eigenvalue*kinetic).det()
                                     -d["K_g0"]*d["K_f0"]*eigenvalue*(eigenvalue-d["m_FP0_squared"]))
    result["general_ratio_EH_kinetic_g"] = sp.cancel(g/rg**4*rg**2-d["K_g0"])
    result["general_ratio_EH_kinetic_f"] = sp.cancel(f/rf**4*rf**2-d["K_f0"])
    # Exact weight-loss identity. A sum of positive denominators fixes its
    # sign on a same-action calibrated Friedmann solution for EITHER p sign.
    ratio, ratio0, p, h = sp.symbols("R R0 p H", real=True)
    loss = p/ratio-p/ratio0
    result["same_action_weight_loss_factorization"] = sp.cancel(loss+p*(ratio**3-ratio0**3)/(ratio*ratio0*(ratio**2+ratio*ratio0+ratio0**2)))
    denominator = ratio*ratio0*(ratio**2+ratio*ratio0+ratio0**2)
    calibrated_lapse = 3*g*h**2-2*p*(ratio**3-ratio0**3)
    result["same_action_weight_loss_calibrated_lapse"] = sp.cancel(2*denominator*loss+3*g*h**2-calibrated_lapse)
    return result


def controls():
    d = derive()
    # Non-unit ratios test the physical-clock kinetic normalization.
    # Choose P_g0=1,P_f0=-2 so S0=-1 and qeff0=4.
    point = {d["R_g0"]: 2, d["R_f0"]: 3, d["p_g"]: 2, d["p_f"]: -6,
             d["G"]: 8, d["F"]: 27}
    return {key: sp.cancel(d[key].subs(point)) for key in
            ("P_g0", "P_f0", "S0", "K_g0", "K_f0", "q_eff0", "m_FP0_squared")}
