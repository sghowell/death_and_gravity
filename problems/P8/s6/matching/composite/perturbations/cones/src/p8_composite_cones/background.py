"""Full-lapse kinematics and null equation after the cone implications."""

import sympy as sp
from p8_composite import background as old
from p8_composite_modes import model as m


def checks():
    d = old.equations()
    y, r = d["y"], d["r"]
    coincident = {old.Nf: old.Ng*y}
    he_difference = d["H_eff"]-d["H_g"]/r
    n, dg_h = sp.symbols("scalar_null Dg_Hg", real=True)
    result = {
        "full_lapse_ratio_derivative_at_coincident_cones": sp.factor(
            (sp.diff(y, old.t)-d["B"]/(old.Ng*old.a)).subs(coincident)),
        "full_lapse_effective_Hubble_on_dynamic_branch": sp.factor(
            he_difference-old.BETA*d["c"]*y*(d["H_f"]-d["H_g"]/y)/(d["s"]*r)),
        "full_lapse_dynamic_error": sp.factor(
            d["H_f"]-d["H_g"]/y-d["B"]/(old.Ng*old.Nf*old.b)),
        "coincident_physical_lapse": sp.factor(d["N_eff"].subs(coincident)-old.Ng*r),
    }
    # y is constant only AFTER c=y and B=0 have both been established.
    # The physical derivative is Dg/r and He=Hg/r on that branch.
    physical_derivative = dg_h/m.r**2
    null_velocity = -m.alpha*m.r**3*n/(2*m.G)
    result["physical_Hubble_monotonicity"] = m.cancel(
        physical_derivative.subs(dg_h, null_velocity)+m.alpha*m.r*n/(2*m.G))
    result["null_coefficient_is_positive"] = m.cancel(m.alpha*m.r/(2*m.G)
                                                      -(m.alpha**2+m.alpha*m.beta*m.y)/(2*m.G))
    # Re-derive the full g null balance from the pinned source terms.
    full_source = d["rho_g_matter"]+d["p_g_matter"]+(d["y"]-d["c"])*d["P"]
    result["literal_full_g_null_source_at_coincident_cones"] = sp.factor(
        (full_source-old.ALPHA*r**3*d["null"]).subs(coincident))
    return result


def coincident_pressure_control():
    """Actual local ODE control if zero vector speed is allowed.

    This constructs a local potential, not an independently frozen free/CD
    solution or a healthy bounce. Its smooth ODE and positive lapse/null
    at the initial point give an ordinary local analytic existence proof.
    """
    from p8_composite import reconstruction

    d = reconstruction.derive()
    y, rho, h = d["y"], d["rho"], d["h"]
    prescribed = (d["X"]+y**2*d["Y"])/(1+y)**2
    velocity = d["yprime"].subs(h, prescribed)
    density_velocity = d["rhoprime"].subs(h, prescribed)
    initial = {y: 1, rho: sp.Rational(1, 2)}
    acceleration = sum(sp.diff(prescribed, z)*v for z, v in
                       ((y, velocity), (rho, density_velocity)))
    return {"h_function": prescribed,
            "coincident_lapse_identity": sp.simplify((d["Nf"]-y*d["Ng"]).subs(h, prescribed)),
            "initial_H": sp.simplify(prescribed.subs(initial)),
            "initial_H_prime": sp.simplify(acceleration.subs(initial)),
            "initial_y_prime": sp.simplify(velocity.subs(initial)),
            "initial_null": d["null"].subs(initial),
            "initial_Ng": sp.simplify(d["Ng"].subs(h, prescribed).subs(initial)),
            "initial_Nf": sp.simplify(d["Nf"].subs(h, prescribed).subs(initial))}
