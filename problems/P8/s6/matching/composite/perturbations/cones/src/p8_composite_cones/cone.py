"""Physical tensor cones from the literal quadratic action, not mass eigenvalues."""

from functools import cache

from p8_composite_modes import model as m
from p8_composite_modes import tensor, vector


@cache
def derive():
    tt, vv = tensor.derive(), vector.derive()
    kg, kf = tt["K_g"]*m.Ne, tt["K_f"]*m.Ne
    ug, uf = tt["U_g_gradient"]/m.Ne, tt["U_f_gradient"]/m.Ne
    cg2, cf2 = m.cancel(ug*m.Ae**2/kg), m.cancel(uf*m.Ae**2/kf)
    cg, cf = m.r/m.s, m.c*m.r/(m.y*m.s)
    branch = m.P-m.alpha*m.beta*m.r**2*m.p
    return {"K_g_T": kg, "K_f_T": kf, "gradient_g_T": ug, "gradient_f_T": uf,
            "c_g_squared": cg2, "c_f_squared": cf2, "c_g": cg, "c_f": cf,
            "weight_g": m.alpha/m.r, "weight_f": m.beta*m.y/m.r,
            "mu": tt["mu"], "Q": branch,
            "c_vector_squared": vv["physical_principal_speed_squared"],
            "Xi": vv["Xi"], "vector_mu": vv["mu"]}


def checks():
    d = derive()
    shear = m.M4*(m.betas[2]+m.betas[3]*m.y)-m.alpha*m.beta**2*m.r*m.p
    return {"literal_g_physical_principal_cone": m.cancel(d["c_g_squared"]-d["c_g"]**2),
            "literal_f_physical_principal_cone": m.cancel(d["c_f_squared"]-d["c_f"]**2),
            "positive_tensor_cone_weights_sum_one": m.cancel(d["weight_g"]+d["weight_f"]-1),
            "weighted_tensor_speeds_equal_matter_speed":
                m.cancel(d["weight_g"]*d["c_g"]+d["weight_f"]*d["c_f"]-1),
            "g_cone_order": m.cancel(d["c_g_squared"]-1
                                     -m.beta*(m.y-m.c)*(m.r+m.s)/m.s**2),
            "f_cone_order": m.cancel(d["c_f_squared"]-1
                                     -m.alpha*(m.c-m.y)*(m.c*m.r+m.y*m.s)/(m.y*m.s)**2),
            "general_pressure_stiffness_split": m.cancel(d["mu"]-m.y*d["Q"]-m.y*(m.c-m.y)*shear),
            "coincident_cone_stiffness": m.cancel(d["mu"].subs(m.c, m.y)-m.y*d["Q"]),
            "coincident_vector_speed": m.cancel(d["c_vector_squared"].subs(m.c, m.y)
                                                -d["vector_mu"]/d["Xi"]),
            "pressure_branch_coincident_vector_zero": m.cancel(
                d["mu"].subs(m.c, m.y).subs(m.p, m.P/(m.alpha*m.beta*m.r**2)))}


def controls():
    d = derive()
    # Algebraic fixtures only, not claimed background solutions.
    faster_g = {m.alpha: 1, m.beta: 1, m.y: 2, m.c: 1}
    faster_f = {m.alpha: 1, m.beta: 1, m.y: 1, m.c: 2}
    omit_pressure = m.y*m.M4*(m.betas[1]+2*m.betas[2]*m.y+m.betas[3]*m.y**2)
    return {"g_can_exceed_matter_cone": d["c_g_squared"].subs(faster_g)-1,
            "f_can_exceed_matter_cone": d["c_f_squared"].subs(faster_f)-1,
            "single_metric_beta_zero_does_not_force_c_equal_y":
                m.cancel(d["c_g_squared"].subs(m.beta, 0)-1),
            "omitted_pressure_falsely_keeps_relative_stiffness": omit_pressure,
            "squared_speed_not_velocity_average": m.cancel(
                d["weight_g"]*d["c_g_squared"]+d["weight_f"]*d["c_f_squared"]-1)}
