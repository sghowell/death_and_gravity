"""Literal own-metric-source TT Hessian and relative ADM-shift coefficient.

Coefficients K and U below are per unit TT tensor contraction. The literal
diagonal exponential audit uses diag(1,-1,0), whose squared norm is two.
No composite-matter pressure term or printed foreign normalization is used.
"""

from functools import cache

import sympy as sp


@cache
def derive():
    a, y, c, n, mass2 = sp.symbols("a y c N_g M_squared", positive=True)
    betas = sp.symbols("beta0:5", real=True)
    delta, z, shift = sp.symbols("relative_TT z relative_ADM_shift", real=True)
    generating = (1+c*z)*(1+y*z)*(1+2*y*sp.cosh(delta/2)*z+y**2*z**2)
    e = sp.Poly(generating, z)
    potential = -2*n*a**3*sum(b*e.nth(i) for i, b in enumerate(betas))
    quadratic = sp.diff(potential, delta, 2).subs(delta, 0)/2
    mu = 2*y*(betas[1]+betas[2]*(c+y)+betas[3]*c*y)
    P = 2*(betas[1]+2*betas[2]*y+betas[3]*y**2)
    block_trace = sp.sqrt((c+y)**2-a**2*y**2*shift**2/n**2)
    shift_generating = (1+y*z)**2*(1+block_trace*z+c*y*z**2)
    ep = sp.Poly(shift_generating, z)
    shifted_potential = -2*n*a**3*sum(b*ep.nth(i) for i, b in enumerate(betas))
    shift_coefficient = sp.diff(shifted_potential, shift, 2).subs(shift, 0)/2
    # The regular positive-root chart has c+y>0.
    shift_coefficient = shift_coefficient.subs(sp.sqrt(c**2+2*c*y+y**2), c+y)
    hg, hf, gd, fd = sp.symbols("H_g H_f gamma_g_dot gamma_f_dot", real=True)
    kg = sp.diag(hg+gd/(2*n), hg-gd/(2*n), hg)
    kf = sp.diag(hf+fd/(2*c*n), hf-fd/(2*c*n), hf)
    eg = mass2*n*a**3*(sp.trace(kg*kg)-sp.trace(kg)**2)/2
    ef = mass2*c*n*a**3*y**3*(sp.trace(kf*kf)-sp.trace(kf)**2)/2
    return {"a": a, "y": y, "c": c, "N_g": n, "M_squared": mass2, "betas": betas,
            "delta": delta, "shift": shift, "generating": generating, "potential": potential,
            "literal_TT_quadratic": quadratic, "mu": mu, "P": P,
            "block_trace": block_trace, "literal_shift_coefficient": shift_coefficient,
            "shift_coefficient": a**5*y**2*P/(2*n*(c+y)),
            "K_g": mass2*a**3/(8*n), "K_f": mass2*a**3*y**3/(8*c*n),
            "gradient_g": mass2*n*a/8, "gradient_f": mass2*c*n*a*y/8,
            "relative_potential": n*a**3*mu/8,
            "literal_g_inertia_norm_two": sp.diff(eg, gd, 2)/2,
            "literal_f_inertia_norm_two": sp.diff(ef, fd, 2)/2}


@cache
def checks():
    d = derive()
    k1, k2 = d["K_g"], d["K_f"]
    zero_extra = {d["betas"][2]: 0, d["betas"][3]: 0}
    return {
        "literal_TT_potential_norm_two": sp.factor(d["literal_TT_quadratic"]+2*d["relative_potential"]),
        "literal_g_EH_inertia": sp.factor(d["literal_g_inertia_norm_two"]-2*k1),
        "literal_f_EH_inertia": sp.factor(d["literal_f_inertia_norm_two"]-2*k2),
        "literal_two_by_two_shift": sp.factor(d["literal_shift_coefficient"]-d["shift_coefficient"]),
        "beta1_actual_spring": sp.factor((d["mu"]-d["y"]*d["P"]).subs(zero_extra)),
        "relative_algebraic_mass": sp.factor(d["relative_potential"]/(k1*k2/(k1+k2))
                                               -d["N_g"]**2*d["mu"]/d["M_squared"]*(1+d["c"]/d["y"]**3)),
        "own_metric_scalar_TT_volume": sp.diff(sp.exp(d["delta"])*sp.exp(-d["delta"]), d["delta"]),
    }


@cache
def spatial_curvature_check():
    z = sp.Symbol("z", real=True)
    h = sp.Function("h")(z)
    metric = sp.diag(sp.exp(h), sp.exp(-h), 1)
    inverse = metric.inv()
    def derivative(value, index):
        return sp.diff(value, z) if index == 2 else sp.S.Zero
    connection = [[[sum(inverse[i, ell]*(derivative(metric[ell, k], j)
                    +derivative(metric[ell, j], k)-derivative(metric[j, k], ell))
                    for ell in range(3))/2 for k in range(3)]
                   for j in range(3)] for i in range(3)]
    ricci = sp.zeros(3)
    for i in range(3):
        for j in range(3):
            ricci[i, j] = sum(derivative(connection[k][i][j], k)
                              -derivative(connection[k][i][k], j)
                              +sum(connection[k][k][ell]*connection[ell][i][j]
                                   -connection[k][j][ell]*connection[ell][i][k]
                                   for ell in range(3)) for k in range(3))
    return {"literal_spatial_TT_gradient": sp.simplify(sp.trace(inverse*ricci)+sp.diff(h, z)**2/2)}


def vector_checks():
    """Exact two-shift Schur and high-k sign; not scalar/vector full health."""
    A, B, k2 = sp.symbols("A B k_squared", positive=True)
    C = sp.Symbol("C_shift", real=True)
    u, v, e, f = sp.symbols("u v Eprime Fprime", real=True)
    lag = A*k2*(e-u)**2+B*k2*(f-v)**2+C*(u-v)**2
    solution = sp.solve((sp.diff(lag, u), sp.diff(lag, v)), (u, v))
    reduced = sp.factor(lag.subs(solution))
    coefficient = A*B*C*k2/(A*B*k2+(A+B)*C)
    return {"literal_two_shift_reduction": sp.factor(reduced-coefficient*(e-f)**2),
            "high_k_relative_kinetic_sign": sp.limit(coefficient, k2, sp.oo)-C}


def controls():
    from . import background
    values = {}
    for c in (1, 4):
        d = background.evaluate(0, c)
        C = d["a"]**5*d["y"]**2*d["P"]/(2*(c+d["y"]))
        values[f"c{c}_center_shift_coefficient"] = sp.factor(C)
        values[f"c{c}_center_f_tensor_speed_squared"] = sp.factor(c**2/d["y"]**2)
    if values["c1_center_shift_coefficient"] != -sp.Rational(128, 3) or values["c4_center_shift_coefficient"] != sp.Rational(8, 3):
        raise ValueError("The actual-family opposite vector-sign controls failed")
    return values
