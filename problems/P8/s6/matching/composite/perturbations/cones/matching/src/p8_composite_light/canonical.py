"""Exact time-dependent two-TT algebra; no unspecified heavy state is discarded."""

from functools import cache

import sympy as sp
from p8_composite_modes import model as m
from p8_composite_modes import tensor


@cache
def physical_coefficients():
    old = tensor.derive()
    k1, k2 = old["K_g"]*m.Ne, old["K_f"]*m.Ne
    total = k1+k2
    relative = k1*k2/total
    w1, w2 = k1/total, k2/total
    cg2, cf2 = (m.r/m.s)**2, (m.c*m.r/(m.y*m.s))**2
    return {"K1": k1, "K2": k2, "K_sum": total, "K_relative": relative,
            "f_sum": sp.sqrt(2*total), "f_relative": sp.sqrt(2*relative),
            "w1": w1, "w2": w2, "cg2": cg2, "cf2": cf2,
            "c_light_squared": m.cancel(w1*cg2+w2*cf2),
            "c_heavy_squared": m.cancel(w2*cg2+w1*cf2),
            "D": sp.sqrt(w1*w2)*(cf2-cg2),
            "mass_squared": m.cancel(old["mu"]/m.s**2*(1/m.G+m.c/(m.F*m.y**3))),
            "relative_stiffness_T": old["U_relative"]/m.Ne}


def B(value, time, omega, theta_sum, cross):
    """Heavy forcing B=2omega(partial_T-theta_sum)+cross; cross=q²D."""
    return 2*omega*(sp.diff(value, time)-theta_sum*value)+cross*value


def B_adjoint(value, time, omega, theta_sum, cross):
    """Formal integration-by-parts adjoint; it is not a retarded inverse."""
    return -2*omega*sp.diff(value, time)+(cross-2*sp.diff(omega, time)-2*omega*theta_sum)*value


@cache
def generic_action():
    t = sp.Symbol("T", real=True)
    light, heavy = sp.Function("l")(t), sp.Function("H")(t)
    ts, tr, om = (sp.Function(name)(t) for name in ("theta_sum", "theta_relative", "omega"))
    light_gradient, heavy_gradient, cross = (sp.Function(name)(t) for name in ("q2_c_light2", "q2_c_heavy2", "q2_D"))
    mass = sp.Function("mass_squared", real=True)(t)
    ns, nr = sp.diff(ts, t)+ts**2, sp.diff(tr, t)+tr**2
    vll, vhh = light_gradient-ns, mass+heavy_gradient-nr-4*om**2
    pre = ((sp.diff(light, t)-ts*light-2*om*heavy)**2+(sp.diff(heavy, t)-tr*heavy)**2)/2
    pre -= (light_gradient*light**2+2*cross*light*heavy+(mass+heavy_gradient)*heavy**2)/2
    post = (sp.diff(light, t)**2+sp.diff(heavy, t)**2-vll*light**2-vhh*heavy**2)/2
    post -= 2*om*heavy*sp.diff(light, t)+(cross-2*om*ts)*light*heavy
    boundary = -(ts*light**2+tr*heavy**2)/2
    return {"t": t, "l": light, "H": heavy, "theta_sum": ts, "theta_relative": tr,
            "omega": om, "cross": cross, "mass_squared": mass,
            "light_gradient": light_gradient, "heavy_gradient": heavy_gradient,
            "N_sum": ns, "N_relative": nr, "V_LL": vll, "V_HH": vhh,
            "pre_boundary_action": pre, "post_boundary_action": post,
            "normalization_boundary": boundary}


def generic_action_checks():
    d = generic_action()
    t, ll, hh, om, ts, cross = (d[key] for key in ("t", "l", "H", "omega", "theta_sum", "cross"))

    def euler(value):
        return sp.diff(sp.diff(d["post_boundary_action"], sp.diff(value, t)), t)-sp.diff(d["post_boundary_action"], value)

    light_equation = sp.diff(ll, t, 2)+d["V_LL"]*ll+B_adjoint(hh, t, om, ts, cross)
    heavy_equation = sp.diff(hh, t, 2)+d["V_HH"]*hh+B(ll, t, om, ts, cross)
    pl, ph = sp.symbols("p_l p_H", real=True)
    velocities = {sp.diff(ll, t): pl+2*om*hh, sp.diff(hh, t): ph}
    hamiltonian = (pl*sp.diff(ll, t)+ph*sp.diff(hh, t)-d["post_boundary_action"]).subs(velocities, simultaneous=True)
    expected_hamiltonian = (pl**2+ph**2+d["V_LL"]*ll**2
                            +(d["mass_squared"]+d["heavy_gradient"]-d["N_relative"])*hh**2)/2
    expected_hamiltonian += 2*om*pl*hh+(cross-2*om*ts)*ll*hh
    f, g = sp.Function("test_f")(t), sp.Function("test_g")(t)
    return {"full_normalization_time_boundary": sp.expand(d["pre_boundary_action"]-d["post_boundary_action"]-sp.diff(d["normalization_boundary"], t)),
            "light_Euler_with_adjoint": sp.expand(euler(ll)-light_equation),
            "heavy_Euler_with_forcing": sp.expand(euler(hh)-heavy_equation),
            "formal_adjoint_with_surface_term": sp.expand(f*B(g, t, om, ts, cross)-B_adjoint(f, t, om, ts, cross)*g-sp.diff(2*om*f*g, t)),
            "full_canonical_Hamiltonian": sp.expand(hamiltonian-expected_hamiltonian)}


def mass_led_pair(light, time, omega, theta_sum, cross, mass_squared):
    """Formal first inverse-mass representative. Requires nonzero mass_squared.

    This does not impose a heavy homogeneous state or prove a Green bound.
    """
    mass_squared = sp.sympify(mass_squared)
    if mass_squared.is_zero is True:
        raise ValueError("The inverse-mass representative requires nonzero mass_squared")
    return -B(light, time, omega, theta_sum, cross)/mass_squared


def operator_checks():
    d = generic_action()
    t, ll, om, ts, cross, mass = (d[key] for key in ("t", "l", "omega", "theta_sum", "cross", "mass_squared"))
    h0 = mass_led_pair(ll, t, om, ts, cross, mass)
    forced = B(ll, t, om, ts, cross)/mass
    remainder = d["V_HH"]-mass
    exact_residual = -sp.diff(forced, t, 2)-remainder*forced
    local_correction = B(ll, t, om, ts, cross)**2/(2*mass)
    correction_euler = sp.diff(sp.diff(local_correction, sp.diff(ll, t)), t)-sp.diff(local_correction, ll)
    hom, response = sp.Function("H_hom")(t), sp.Function("retarded_response")(t)
    coupled_light = sp.diff(ll, t, 2)+d["V_LL"]*ll+B_adjoint(hom-response, t, om, ts, cross)
    retarded_light = (sp.diff(ll, t, 2)+d["V_LL"]*ll
                      -B_adjoint(response, t, om, ts, cross)+B_adjoint(hom, t, om, ts, cross))
    return {"exact_mass_led_heavy_residual": sp.expand(sp.diff(h0, t, 2)+d["V_HH"]*h0+B(ll, t, om, ts, cross)-exact_residual),
            "formal_local_action_first_inverse_mass": sp.expand(correction_euler+B_adjoint(forced, t, om, ts, cross)),
            "retarded_equation_keeps_homogeneous_data": sp.expand(coupled_light-retarded_light)}


def routh_checks():
    t = sp.Symbol("T", real=True)
    u, delta, heavy = (sp.Function(name)(t) for name in ("weighted_u", "relative_delta", "Routh_H"))
    ks, fr, w = (sp.Function(name, positive=True)(t) for name in ("K_sum", "f_relative", "w2"))
    mass = sp.Function("mass_squared", real=True)(t)
    relative = fr**2/2
    stiffness = mass*relative
    action = ks*(sp.diff(u, t)-sp.diff(w, t)*delta)**2+relative*sp.diff(delta, t)**2-stiffness*delta**2
    charge = 2*ks*(sp.diff(u, t)-sp.diff(w, t)*delta)
    conserved = sp.Symbol("J", real=True)
    heavy_euler = sp.diff(sp.diff(action, sp.diff(delta, t)), t)-sp.diff(action, delta)
    fixed_charge = heavy_euler.subs(sp.diff(u, t), conserved/(2*ks)+sp.diff(w, t)*delta)
    normalized = fixed_charge.subs(delta, heavy/fr).doit()/fr
    expected = sp.diff(heavy, t, 2)+(mass-sp.diff(fr, t, 2)/fr)*heavy+sp.diff(w, t)*conserved/fr
    return {"exact_zero_k_Noether_charge": sp.expand(sp.diff(sp.diff(action, sp.diff(u, t)), t)-sp.diff(charge, t)),
            "fixed_charge_Routh_equation": sp.simplify(normalized-expected)}


def physical_checks():
    d = physical_coefficients()
    ks, w = sp.symbols("K_sum w2", positive=True)
    ud, dd, wd, delta = sp.symbols("udot deltadot wdot delta", real=True)
    hdot, gdot = ud-w*dd-wd*delta, ud+(1-w)*dd-wd*delta
    return {"weighted_kinetic_with_rotating_eigenvector": sp.expand(ks*(1-w)*hdot**2+ks*w*gdot**2
                                                                 -ks*(ud-wd*delta)**2-ks*w*(1-w)*dd**2),
            "locked_principal_speed": m.cancel(d["c_light_squared"]-(m.r/m.s)**2*(m.G+m.F*m.c*m.y)/(m.G+m.F*m.y**3/m.c)),
            "frozen_relative_algebraic_mass": m.cancel(d["mass_squared"]-d["relative_stiffness_T"]/d["K_relative"]),
            "weight_sum": m.cancel(d["w1"]+d["w2"]-1)}


def controls():
    """Exact omission witnesses, not additional background or spectral claims."""
    t = sp.Symbol("T", real=True)
    omitted_rotation = B(t, t, 1, 0, 0)-B(t, t, 0, 0, 0)
    omitted_hom = B_adjoint(sp.sin(t), t, 1, 0, 0).subs(t, 0)
    h0 = mass_led_pair(t, t, 1, 0, 0, t**2)
    lost_mass_derivatives = sp.diff(h0, t, 2).subs(t, 1)
    ll, ld, ts = sp.symbols("l lprime theta_sum", real=True)
    omitted_boundary_momentum = sp.diff((ld-ts*ll)**2/2-ld**2/2, ld)
    return {"omitted_rotating_weight_heavy_forcing": omitted_rotation,
            "omitted_homogeneous_heavy_data_light_forcing": omitted_hom,
            "omitted_nonconstant_mass_derivatives_residual": lost_mass_derivatives,
            "omitted_time_boundary_changes_momentum": omitted_boundary_momentum.subs({ll: 1, ts: 1})}
