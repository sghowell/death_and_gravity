"""Full lapse/scale equations before their sign and global-tail analysis.

B-sector convention: +---, R=-6(DH+2H²), Einstein_00=+3H² in proper
time, EH=-G R/2. There is no direct e/v matter source; it acts on u only.
"""

from functools import cache

import sympy as sp
from p8_trimetric import model


@cache
def derive():
    g, n, a, nu, au = sp.symbols("G_i n_i a_i n_u a_u", positive=True)
    p, b, h, dh = sp.symbols("p_i b_i H_i D_i_H_i", real=True)
    ratio, cone = au/a, au*n/(nu*a)
    lapse = 3*g*h**2-2*p*ratio**3-2*b
    space = g*(2*dh+3*h**2)-2*p*ratio**3/cone-2*b
    null = dh-p*ratio**3*(1/cone-1)/g
    return {"G": g, "n_i": n, "a_i": a, "n_u": nu, "a_u": au,
            "p": p, "b": b, "H": h, "D_H": dh, "R": ratio, "c": cone,
            "lapse_residual": lapse, "spatial_residual": space, "null_residual": null,
            "proper_clock_ratio": n/nu,
            "physical_tensor_kinetic": g/(cone*ratio**2),
            "physical_tensor_gradient": g*cone/ratio**2}


@cache
def checks():
    d = derive()
    g, n, a, nu, au, p, b, h, dh = (d[key] for key in ("G", "n_i", "a_i", "n_u", "a_u", "p", "b", "H", "D_H"))
    e, u = sp.diag(n, a, a, a), sp.diag(nu, au, au, au)
    einstein = sp.diag(3*h**2/n**2, *([-2*dh-3*h**2]*3))
    for index in (1, 2, 3):
        einstein[index, index] /= a**2
    full = model.euler_maps(e, sp.eye(4), u, pg=p, pf=0, bg=b,
                            g=g, einstein_g=einstein, epsilon=0)["E_e"]
    # A genuinely separate lapse-retaining EH minisuperspace variation.
    adot, addot, ndot = sp.symbols("a_dot a_ddot n_dot", real=True)
    mini = -3*g*a*adot**2/n-2*p*(n*au**3+3*nu*au**2*a)-2*b*n*a**3
    momentum = sp.diff(mini, adot)
    momentum_dot = sp.diff(momentum, a)*adot+sp.diff(momentum, adot)*addot+sp.diff(momentum, n)*ndot
    scale_euler = sp.diff(mini, a)-momentum_dot
    clock = {h: adot/(n*a), dh: addot/(n**2*a)-adot**2/(n**2*a**2)-adot*ndot/(n**3*a)}
    result = {"literal_full_covariant_lapse": sp.cancel(full[0, 0]/a**3-d["lapse_residual"]),
              "lapse_retaining_action_variation": sp.cancel(sp.diff(mini, n)/a**3-d["lapse_residual"].subs(clock, simultaneous=True)),
              "scale_retaining_action_variation": sp.cancel(scale_euler/(3*n*a**2)-d["spatial_residual"].subs(clock, simultaneous=True)),
              "undivided_null_identity": sp.cancel(d["spatial_residual"]-d["lapse_residual"]-2*g*d["null_residual"]),
              "ordered_clock_dictionary": sp.cancel(d["proper_clock_ratio"]-d["c"]/d["R"]),
              "physical_TT_characteristic_ratio": sp.cancel(d["physical_tensor_gradient"]/d["physical_tensor_kinetic"]-d["c"]**2)}
    for index in (1, 2, 3):
        result[f"literal_full_covariant_spatial_{index}"] = sp.cancel(full[index, index]/(n*a**2)-d["spatial_residual"])
    # Recover both physical coefficients from actual volume and derivatives.
    measure = n*a**3/(nu*au**3)
    result["physical_TT_kinetic_pullback"] = sp.cancel(g*measure*(nu/n)**2-d["physical_tensor_kinetic"])
    result["physical_TT_gradient_pullback"] = sp.cancel(g*measure*(au/a)**2-d["physical_tensor_gradient"])
    return result
