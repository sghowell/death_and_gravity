"""Literal HR plus one composite-coupled canonical scalar, in P8 +---.

Only the gravitational symbols/conventions are reused from the frozen HR
input. Its separate matter hypothesis and physical-g verdict are not reused.
"""

from functools import cache

import sympy as sp
from p8_bimetric_general import background as hr

t, a, b, Ng, Nf = hr.t, hr.a, hr.b, hr.Ng, hr.Nf
G, F, M4, BETAS = hr.MG2, hr.MF2, hr.M4, hr.BETAS
ALPHA, BETA = sp.symbols("alpha beta", positive=True)
chi = sp.Function("chi")(t)
potential = sp.Function("V")(chi)


@cache
def equations():
    y, c = b/a, Nf/Ng
    r, s = ALPHA+BETA*y, ALPHA+BETA*c
    ae, ne = ALPHA*a+BETA*b, ALPHA*Ng+BETA*Nf
    hg, hf, he = sp.diff(a, t)/(Ng*a), sp.diff(b, t)/(Nf*b), sp.diff(ae, t)/(ne*ae)
    z = sp.Symbol("elementary_z")
    elementary = sp.Poly((1+c*z)*(1+y*z)**3, z)
    interaction = -M4*Ng*a**3*sum(BETAS[index]*elementary.nth(index) for index in range(5))
    gravity = -3*G*a*sp.diff(a, t)**2/Ng-3*F*b*sp.diff(b, t)**2/Nf+interaction
    matter = ae**3*sp.diff(chi, t)**2/(2*ne)-ne*ae**3*potential
    lagrangian = gravity+matter
    rho = sp.diff(chi, t)**2/(2*ne**2)+potential
    pressure = sp.diff(chi, t)**2/(2*ne**2)-potential
    polynomial = M4*(BETAS[1]+2*BETAS[2]*y+BETAS[3]*y**2)
    ug = M4*(BETAS[0]+3*BETAS[1]*y+3*BETAS[2]*y**2+BETAS[3]*y**3)
    uf = M4*(BETAS[4]+3*BETAS[3]/y+3*BETAS[2]/y**2+BETAS[1]/y**3)
    pg, pf = -ug+(y-c)*polynomial, -uf+(c-y)*polynomial/(c*y**3)
    rg, rf = ALPHA*r**3*rho, BETA*r**3*rho/y**3
    mgp, mfp = ALPHA*s*r**2*pressure, BETA*s*r**2*pressure/(c*y**2)

    def el(field):
        return sp.diff(lagrangian, field)-sp.diff(sp.diff(lagrangian, sp.diff(field, t)), t)

    expected = {"EL_Ng": 3*G*hg**2-ug-rg, "EL_Nf": 3*F*hf**2-uf-rf,
                "EL_a": G*(2*sp.diff(hg, t)/Ng+3*hg**2)+pg+mgp,
                "EL_b": F*(2*sp.diff(hf, t)/Nf+3*hf**2)+pf+mfp}
    scalar = sp.diff(sp.diff(chi, t)/ne, t)/ne+3*he*sp.diff(chi, t)/ne+sp.diff(potential, chi)
    conservation = sp.diff(rho, t)/ne+3*he*(rho+pressure)
    return {"L": lagrangian, "gravity": gravity, "interaction": interaction, "matter": matter,
            "y": y, "c": c, "r": r, "s": s, "a_eff": ae, "N_eff": ne,
            "H_g": hg, "H_f": hf, "H_eff": he,
            "rho": rho, "p": pressure, "null": rho+pressure,
            "U_g": ug, "U_f": uf, "P_g": pg, "P_f": pf,
            "rho_g_matter": rg, "p_g_matter": mgp, "rho_f_matter": rf, "p_f_matter": mfp,
            "P": polynomial, "Q": polynomial-ALPHA*BETA*r**2*pressure,
            "B": Ng*sp.diff(b, t)-Nf*sp.diff(a, t),
            "EL_Ng": sp.factor(el(Ng)/a**3), "EL_Nf": sp.factor(el(Nf)/b**3),
            "EL_a": sp.factor(el(a)/(3*Ng*a**2)), "EL_b": sp.factor(el(b)/(3*Nf*b**2)),
            "EL_chi": el(chi), "scalar_equation": scalar, "conservation": conservation,
            "expected": expected}


def variation_checks():
    data = equations()
    checks = {key: sp.factor(data[key]-value) for key, value in data["expected"].items()}
    checks["scalar_lapse_retaining_variation"] = sp.factor(
        data["EL_chi"]+data["N_eff"]*data["a_eff"]**3*data["scalar_equation"])
    checks["scalar_energy_identity_without_dividing_chidot"] = sp.factor(
        data["conservation"]-sp.diff(chi, t)*data["scalar_equation"]/data["N_eff"])
    checks["unchanged_gravitational_interaction_only"] = sp.factor(data["interaction"]-hr.equations()["interaction"])
    checks["positive_kinetic_null"] = sp.factor(data["null"]-sp.diff(chi, t)**2/data["N_eff"]**2)
    return checks


def bianchi_checks():
    data = equations()
    y, c, r, s, q, branch = (data[key] for key in ("y", "c", "r", "s", "Q", "B"))
    checks = {}
    for label, lapse, h, el_a, el_n, weight, target in (
        ("g", Ng, data["H_g"], "EL_a", "EL_Ng", ALPHA*r**3*s, 3*q*branch/(Ng**2*a)),
        ("f", Nf, data["H_f"], "EL_b", "EL_Nf", BETA*r**3*s/(c*y**3),
         -3*q*branch/(Ng*Nf*a*c*y**3)),
    ):
        rho = data["U_"+label]+data["rho_"+label+"_matter"]
        p = data["P_"+label]+data["p_"+label+"_matter"]
        balance = sp.diff(rho, t)/lapse+3*h*(rho+p)
        checks[label+"_undivided_source_aware_Bianchi"] = sp.factor(balance-weight*data["conservation"]-target)
        checks[label+"_full_Noether_identity"] = sp.factor(
            sp.diff(data[el_n], t)/lapse-3*h*(data[el_a]-data[el_n])+balance)
    rg, pg, rf, pf = (data[key] for key in ("rho_g_matter", "p_g_matter", "rho_f_matter", "p_f_matter"))
    checks["positive_weighted_shared_matter_null"] = sp.factor(rg+pg+c*y**3*(rf+pf)-s*r**3*data["null"])
    checks["g_pressure_branch_null_rewrite"] = sp.factor(rg+pg+(y-c)*data["P"]-ALPHA*r**3*data["null"]-(y-c)*q)
    checks["f_pressure_branch_null_rewrite"] = sp.factor(
        rf+pf+(c-y)*data["P"]/(c*y**3)-BETA*r**3*data["null"]/y**3-(c-y)*q/(c*y**3))
    return checks
