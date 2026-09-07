"""Actual regular solutions versus explicitly outside-hypothesis controls."""

from functools import cache

import sympy as sp

from . import edge


def profile_equations(time, Gs, edges, betas, scales, lapses, scalar_profiles):
    """Full homogeneous lapse/scale/current residuals for separate free fields.

    Each supplied scalar trajectory belongs to its own matter action, even
    when two on-shell profiles coincide. This is not a shared coupled field.
    """
    count = len(Gs)
    d = edge.derive()
    rho_v, pressure_v = [sp.Integer(0)]*count, [sp.Integer(0)]*count
    for (i, j), beta in zip(edges, betas, strict=True):
        point = dict(zip(d["beta"], beta, strict=True))
        point.update({d["y"]: scales[j]/scales[i], d["c"]: lapses[j]/lapses[i]})
        for vertex, suffix in ((i, "i"), (j, "j")):
            rho_v[vertex] += d[f"rho_{suffix}"].subs(point)
            pressure_v[vertex] += d[f"pressure_{suffix}"].subs(point)
    h = tuple(sp.diff(a, time)/(n*a) for a, n in zip(scales, lapses, strict=True))
    dh = tuple(sp.diff(hi, time)/n for hi, n in zip(h, lapses, strict=True))
    rho = tuple(sp.diff(phi, time)**2/(2*n**2) for phi, n in zip(scalar_profiles, lapses, strict=True))
    result = {}
    for i in range(count):
        result[f"lapse_{i}"] = sp.simplify(3*Gs[i]*h[i]**2-rho[i]-rho_v[i])
        result[f"spatial_{i}"] = sp.simplify(Gs[i]*(2*dh[i]+3*h[i]**2)+rho[i]+pressure_v[i])
        result[f"separate_scalar_current_{i}"] = sp.simplify(sp.diff(scales[i]**3*sp.diff(scalar_profiles[i], time)/lapses[i], time))
    return {"H": h, "D_H": dh, "matter_rho": rho, "interaction_rho": rho_v,
            "interaction_pressure": pressure_v, "residuals": result}


@cache
def stiff_tree():
    t = sp.Symbol("T", positive=True)
    gs = (sp.Integer(1), sp.Integer(0), sp.Integer(2), sp.Integer(3))
    edges = ((0, 1), (2, 1), (1, 3))
    beta = (-3, 1, 0, 0, -1)
    scale = t**sp.Rational(1, 3)
    phis = tuple(sp.sqrt(2*g/3)*sp.log(t) for g in gs)
    profiles = profile_equations(t, gs, edges, (beta,)*3, (scale,)*4, (sp.Integer(1),)*4, phis)
    return {"T": t, "G": gs, "edges": edges, "beta": beta, "scalar_profiles": phis, **profiles}


@cache
def mixed_algebraic_tree():
    t = sp.Symbol("T", positive=True)
    gs, edges = (1, 1, 1), ((1, 0), (2, 0))
    betas = ((0, 1, -sp.Rational(1, 2), 0, sp.Rational(1, 2)), (-3, 1, 0, 0, -1))
    scale, phi = t**sp.Rational(1, 3), sp.sqrt(sp.Rational(2, 3))*sp.log(t)
    profiles = profile_equations(t, gs, edges, betas, (scale,)*3,
                                (sp.Integer(1), 1/(3*t), sp.Integer(1)), (phi, sp.Integer(0), phi))
    h, hp = profiles["H"][0], profiles["D_H"][0]
    n = 2*(profiles["matter_rho"][0]+profiles["matter_rho"][2])
    return {"T": t, "G": gs, "edges": edges, "betas": betas,
            "P_values": (sp.Integer(0), sp.Integer(2)), "component": (0, 2),
            "A_component": sp.Integer(2), "A_false_all": sp.Integer(3),
            "correct_component_null_residual": sp.simplify(-4*hp-n),
            "false_all_vertex_null_residual": sp.simplify(-6*hp-n), "H_root": h, **profiles}


@cache
def zero_target_exception():
    t = sp.Symbol("T", real=True)
    profiles = profile_equations(t, (0, 1), ((0, 1),), ((0, 0, 0, 0, 0),),
                                (1+t**2, sp.Integer(1)), (sp.Integer(1),)*2, (sp.Integer(0),)*2)
    return {"T": t, "H_root_prime_at_zero": sp.diff(profiles["H"][0], t).subs(t, 0), **profiles}


@cache
def nec_violating_control():
    """Actual GR fluid reconstruction, not a positive-canonical scalar.

    With equal two-vertex metrics and tuned beta1 edges, every interaction
    stress vanishes. These conserved fluid profiles solve both Einstein
    equations; near T=0 they can be realized by a phantom-sign scalar.
    """
    t = sp.Symbol("T", real=True)
    a = 1+t**2
    h = sp.diff(a, t)/a
    dh = sp.diff(h, t)
    rho = 3*h**2
    pressure = -2*dh-rho
    return {"T": t, "H": h, "Hprime": dh, "rho": rho, "pressure": pressure,
            "null_at_bounce": (rho+pressure).subs(t, 0),
            "lapse_residual": sp.simplify(3*h**2-rho),
            "null_residual": sp.simplify(-2*dh-rho-pressure),
            "separate_conservation": sp.simplify(sp.diff(rho, t)+3*h*(rho+pressure))}


@cache
def checks():
    first, mixed, zero, nec = stiff_tree(), mixed_algebraic_tree(), zero_target_exception(), nec_violating_control()
    result = {f"stiff_{key}": value for key, value in first["residuals"].items()}
    result.update({f"mixed_{key}": value for key, value in mixed["residuals"].items()})
    result["actual_mixed_fixed_component_null"] = mixed["correct_component_null_residual"]
    result.update({f"zero_target_{key}": value for key, value in zero["residuals"].items()})
    result.update({f"NEC_violating_{key}": nec[key] for key in ("lapse_residual", "null_residual", "separate_conservation")})
    return result
