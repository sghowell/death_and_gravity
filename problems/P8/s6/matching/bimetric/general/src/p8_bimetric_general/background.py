"""Two lapse equations for arbitrary-beta HR and separate matter sectors.

P8 +--- curvature and positive-root branch. Each matter action couples to
its own metric only; no common field is coupled to both metrics.
"""

from functools import cache

import sympy as sp

t = sp.Symbol("t", real=True)
a, b, Ng, Nf = (sp.Function(name, positive=True)(t) for name in ("a", "b", "N_g", "N_f"))
KG, KF, VG, VF = (sp.Function(name)(t) for name in ("K_g", "K_f", "V_g", "V_f"))
MG2, MF2, M4 = sp.symbols("M_g_squared M_f_squared m_fourth", positive=True)
BETAS = sp.symbols("beta0:5", real=True)


def ricci(scale, lapse):
    return -6*(sp.diff(scale, t, 2)/(scale*lapse**2)
               +sp.diff(scale, t)**2/(scale**2*lapse**2)
               -sp.diff(scale, t)*sp.diff(lapse, t)/(scale*lapse**3))


def interaction_polynomial(y):
    return M4*(BETAS[1]+2*BETAS[2]*y+BETAS[3]*y**2)


@cache
def equations():
    beta0, beta1, beta2, beta3, beta4 = BETAS
    y, c = b/a, Nf/Ng
    hg, hf = sp.diff(a, t)/(Ng*a), sp.diff(b, t)/(Nf*b)
    interaction = -M4*(beta0*Ng*a**3+beta1*(Nf*a**3+3*Ng*a**2*b)
                        +3*beta2*(Nf*a**2*b+Ng*a*b**2)
                        +beta3*(3*Nf*a*b**2+Ng*b**3)+beta4*Nf*b**3)
    matter_g, matter_f = a**3*KG/(2*Ng)-Ng*a**3*VG, b**3*KF/(2*Nf)-Nf*b**3*VF
    lagrangian = -3*MG2*a*sp.diff(a, t)**2/Ng-3*MF2*b*sp.diff(b, t)**2/Nf+interaction+matter_g+matter_f
    rhog, pg, rhof, pf = KG/(2*Ng**2)+VG, KG/(2*Ng**2)-VG, KF/(2*Nf**2)+VF, KF/(2*Nf**2)-VF
    rg = M4*(beta0+3*beta1*y+3*beta2*y**2+beta3*y**3)
    pv_g = -M4*(beta0+beta1*(2*y+c)+beta2*(y**2+2*c*y)+beta3*c*y**2)
    rf = M4*(beta4+3*beta3/y+3*beta2/y**2+beta1/y**3)
    pv_f = -M4*(beta4+beta3*(2/y+1/c)+beta2*(1/y**2+2/(c*y))+beta1/(c*y**2))

    def el(field):
        return sp.diff(lagrangian, field)-sp.diff(sp.diff(lagrangian, sp.diff(field, t)), t)

    expected = {"EL_Ng": 3*MG2*hg**2-rhog-rg, "EL_Nf": 3*MF2*hf**2-rhof-rf,
                "EL_a": MG2*(2*sp.diff(hg, t)/Ng+3*hg**2)+pg+pv_g,
                "EL_b": MF2*(2*sp.diff(hf, t)/Nf+3*hf**2)+pf+pv_f}
    return {"L": lagrangian, "interaction": interaction, "y": y, "c": c, "H_g": hg, "H_f": hf,
            "rho_g": rhog, "p_g": pg, "rho_f": rhof, "p_f": pf,
            "rho_g_interaction": rg, "p_g_interaction": pv_g,
            "rho_f_interaction": rf, "p_f_interaction": pv_f,
            "EL_Ng": sp.factor(el(Ng)/a**3), "EL_Nf": sp.factor(el(Nf)/b**3),
            "EL_a": sp.factor(el(a)/(3*Ng*a**2)), "EL_b": sp.factor(el(b)/(3*Nf*b**2)),
            "expected": expected}


def variation_checks():
    data = equations()
    checks = {key: sp.simplify(data[key]-value) for key, value in data["expected"].items()}
    z = sp.Symbol("z")
    generating = sp.Poly((1+z*data["c"])*(1+z*data["y"])**3, z)
    checks["full_elementary_symmetric_potential"] = sp.expand(
        data["interaction"]+M4*Ng*a**3*sum(beta*generating.nth(index) for index, beta in enumerate(BETAS)))
    for label, scale, lapse, mass in (("g", a, Ng, MG2), ("f", b, Nf, MF2)):
        raw = -mass*lapse*scale**3*ricci(scale, lapse)/2
        boundary = 3*mass*scale**2*sp.diff(scale, t)/lapse
        checks[label+"_Einstein_boundary"] = sp.simplify(raw-sp.diff(boundary, t)+3*mass*scale*sp.diff(scale, t)**2/lapse)
        checks[label+"_matter_null_stress"] = sp.simplify(
            data["rho_"+label]+data["p_"+label]-(KG if label == "g" else KF)/lapse**2)
    return checks


def bianchi_checks():
    data = equations()
    y, c, hg, hf = (data[key] for key in ("y", "c", "H_g", "H_f"))
    polynomial = interaction_polynomial(y)
    numerator = Ng*sp.diff(b, t)-Nf*sp.diff(a, t)
    checks = {}
    for label, lapse, hubble, el_scale, el_lapse in (("g", Ng, hg, "EL_a", "EL_Ng"), ("f", Nf, hf, "EL_b", "EL_Nf")):
        matter = sp.diff(data["rho_"+label], t)/lapse+3*hubble*(data["rho_"+label]+data["p_"+label])
        balance = sp.diff(data["rho_"+label+"_interaction"], t)/lapse+3*hubble*(data["rho_"+label+"_interaction"]+data["p_"+label+"_interaction"])
        checks[label+"_Noether_matter_balance"] = sp.simplify(
            sp.diff(data[el_lapse], t)/lapse-3*hubble*(data[el_scale]-data[el_lapse])+matter+balance)
        target = 3*polynomial*numerator/(Ng**2*a) if label == "g" else -3*polynomial*numerator/(Ng*Nf*a*c*y**3)
        checks[label+"_undivided_Bianchi_factor"] = sp.factor(balance-target)
    checks["g_interaction_null_factor"] = sp.factor(data["rho_g_interaction"]+data["p_g_interaction"]-(y-c)*polynomial)
    checks["f_interaction_null_factor"] = sp.factor(data["rho_f_interaction"]+data["p_f_interaction"]-(c-y)*polynomial/(c*y**3))
    checks["positive_weighted_interaction_cancellation"] = sp.factor(
        data["rho_g_interaction"]+data["p_g_interaction"]+c*y**3*(data["rho_f_interaction"]+data["p_f_interaction"]))
    checks["dynamic_Hubble_identity_without_dividing_H"] = sp.factor(hf-hg/y-numerator/(Ng*Nf*b))
    return checks
