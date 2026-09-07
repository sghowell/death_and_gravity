"""Exact two-TT canonical map, including moving weights and time boundary.

All derivative coefficients here use u=T/tau and M=tau=1 for the algebra;
restore theta,omega with tau^-1 and N,mass_squared with tau^-2. The physical
TT amplitudes are unchanged by that notation. Mass is not assumed positive.
"""

from functools import cache

import sympy as sp

from . import background


@cache
def derive():
    d = background.derive()
    u, c, a, b, y = (d[key] for key in ("u", "c", "a", "b", "y"))
    k1, k2 = a**3/8, b**3/(8*c)
    total, relative = k1+k2, k1*k2/(k1+k2)
    w1, w2 = k1/total, k2/total
    sc = sp.sqrt(8*c)*d["d"]**6/(c*d["d"]**12+8)
    ts = 6*u*(c*d["d"]**12-8)/(d["d"]*(c*d["d"]**12+8))
    tr = -ts
    om = -24*sp.sqrt(2*c)*u*d["d"]**5/(c*d["d"]**12+8)
    ns, nr = sp.diff(ts, u)+ts**2, sp.diff(tr, u)+tr**2
    cg2, cf2 = sp.S.One, c**2/y**2
    cl, ch = w1+cf2*w2, w2+cf2*w1
    cross = sc*(cf2-1)
    mass = 2*y*d["h_u"]*(y**3+c)/(c*(c-y))
    momentum = sp.Symbol("kbar", nonnegative=True)
    q2 = momentum**2/a**2
    return {"u": u, "c": c, "kbar": momentum, "q_squared": q2,
            "K1": k1, "K2": k2, "K_sum": total, "K_relative": relative,
            "f_sum_squared": 2*total, "f_relative_squared": 2*relative,
            "w1": w1, "w2": w2, "sqrt_weights": sc,
            "theta_sum": ts, "theta_relative": tr, "omega": om,
            "N_sum": ns, "N_relative": nr, "cg2": cg2, "cf2": cf2,
            "c_light_squared": cl, "c_heavy_squared": ch, "D": cross,
            "mass_squared": mass, "V_LL": q2*cl-ns,
            "V_HH": mass+q2*ch-nr-4*om**2}


def B(value, time, omega, theta_sum, cross):
    return 2*omega*(sp.diff(value, time)-theta_sum*value)+cross*value


def B_adjoint(value, time, omega, theta_sum, cross):
    return -2*omega*sp.diff(value, time)+(cross-2*sp.diff(omega, time)-2*omega*theta_sum)*value


@cache
def action_checks():
    t = sp.Symbol("T", real=True)
    ll, hh = sp.Function("l")(t), sp.Function("Q")(t)
    ts, tr, om, gl, gh, cross, mass = (sp.Function(name)(t) for name in ("thetaS", "thetaR", "omega", "gradientL", "gradientH", "cross", "mass_squared"))
    ns, nr = sp.diff(ts, t)+ts**2, sp.diff(tr, t)+tr**2
    pre = ((sp.diff(ll, t)-ts*ll-2*om*hh)**2+(sp.diff(hh, t)-tr*hh)**2
           -gl*ll**2-2*cross*ll*hh-(mass+gh)*hh**2)/2
    post = (sp.diff(ll, t)**2+sp.diff(hh, t)**2-(gl-ns)*ll**2-(mass+gh-nr-4*om**2)*hh**2)/2
    post -= 2*om*hh*sp.diff(ll, t)+(cross-2*om*ts)*ll*hh
    boundary = -(ts*ll**2+tr*hh**2)/2
    def euler(field):
        return sp.diff(sp.diff(post, sp.diff(field, t)), t)-sp.diff(post, field)
    return {
        "full_normalization_time_boundary": sp.expand(pre-post-sp.diff(boundary, t)),
        "light_equation_keeps_adjoint_connection": sp.expand(euler(ll)-sp.diff(ll, t, 2)-(gl-ns)*ll-B_adjoint(hh, t, om, ts, cross)),
        "heavy_equation_keeps_rotating_forcing": sp.expand(euler(hh)-sp.diff(hh, t, 2)-(mass+gh-nr-4*om**2)*hh-B(ll, t, om, ts, cross)),
        "formal_adjoint_surface_term": sp.expand(hh*B(ll, t, om, ts, cross)-B_adjoint(hh, t, om, ts, cross)*ll-sp.diff(2*om*ll*hh, t)),
    }


@cache
def checks():
    d = derive()
    u, c = d["u"], d["c"]
    bg = background.derive()
    ks, w, ud, dd, wd, delta = sp.symbols("Ks w ud dd wd Delta", real=True)
    xdot, ydot = ud-w*dd-wd*delta, ud+(1-w)*dd-wd*delta
    values = {
        "actual_sum_normalization": d["theta_sum"]-sp.diff(d["K_sum"], u)/(2*d["K_sum"]),
        "actual_relative_normalization": d["theta_relative"]-sp.diff(d["K_relative"], u)/(2*d["K_relative"]),
        "actual_rotating_weight": d["omega"]-sp.diff(d["w2"], u)/(2*d["sqrt_weights"]),
        "positive_square_root_weights": d["sqrt_weights"]**2-d["w1"]*d["w2"],
        "literal_weighted_kinetic_map": ks*(1-w)*xdot**2+ks*w*ydot**2-ks*(ud-wd*delta)**2-ks*w*(1-w)*dd**2,
        "actual_relative_algebraic_mass": d["mass_squared"]-bg["y"]*bg["P"]*(1+c/bg["y"]**3),
        "locked_light_gradient": d["c_light_squared"]-(1+c*bg["y"])/(1+bg["y"]**3/c),
        "full_high_frequency_gradient_eigenvalues": d["c_light_squared"]*d["c_heavy_squared"]-d["D"]**2-d["cf2"],
        "full_high_frequency_gradient_trace": d["c_light_squared"]+d["c_heavy_squared"]-1-d["cf2"],
    }
    return {key: sp.factor(value) for key, value in values.items()}


@cache
def centers():
    d = derive()
    keys = ("mass_squared", "N_sum", "N_relative", "omega", "c_light_squared", "c_heavy_squared", "D")
    result = {key: sp.factor(d[key].subs(d["u"], 0)) for key in keys}
    result["omega_prime"] = sp.factor(sp.diff(d["omega"], d["u"]).subs(d["u"], 0))
    result["zero_k_heavy_diagonal"] = sp.factor(result["mass_squared"]-result["N_relative"])
    return result


def center_checks():
    c, d = background.C, centers()
    return {"mass_center": sp.factor(d["mass_squared"]-16*(c+8)/(c*(c-2))),
            "Nsum_center": sp.factor(d["N_sum"]-6*(c-8)/(c+8)),
            "Nrelative_center": sp.factor(d["N_relative"]-6*(8-c)/(c+8)),
            "rotation_center_zero": d["omega"],
            "rotation_derivative_not_omitted": sp.factor(d["omega_prime"]+24*sp.sqrt(2*c)/(c+8)),
            "locked_center_cone": sp.factor(d["c_light_squared"]-c*(2*c+1)/(c+8)),
            "c4_positive_diagonal_not_a_rolling_gap": d["zero_k_heavy_diagonal"].subs(c, 4)-22,
            "c4_light_cone_wider": d["c_light_squared"].subs(c, 4)-3}
