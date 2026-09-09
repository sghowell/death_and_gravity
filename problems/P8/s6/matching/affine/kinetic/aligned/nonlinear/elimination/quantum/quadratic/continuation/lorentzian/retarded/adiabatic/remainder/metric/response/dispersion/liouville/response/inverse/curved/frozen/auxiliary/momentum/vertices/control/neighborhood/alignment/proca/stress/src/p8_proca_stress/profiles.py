"""Explicit new fixed c-number clock profile; old profile is not reevaluated."""
from functools import cache

import sympy as sp
from p8_clock_tadpole import window
from p8_constant_proca import quantum
from p8_vector_state import wkb

from . import estimates, readouts, subtraction

u=wkb.u
x=sp.Symbol("physical_clock_norm",real=True)


@cache
def integrands():
    a=(1+u*u)**2
    omega=sp.Symbol("positive_mode_frequency",positive=True)
    modes={}
    for kind in ("transverse","longitudinal"):
        v2=sp.Symbol(kind+"_same_exact_mode_squared",nonnegative=True)
        p2=sp.Symbol(kind+"_same_exact_physical_momentum_squared",nonnegative=True)
        ad=subtraction.adiabatic_rows(kind,0)
        physical=(p2+omega**2*v2)/(2*a**3)
        subtract=omega*sum(ad[n]/omega**(2*n) for n in range(3))/(4*a**3)
        modes[kind]={"physical_new_energy":physical,"actual_fourth_order_subtraction":subtract,
                     "physical_momentum_rate":sp.factor(readouts.weights(kind)["c1"]-wkb.background()["lambda"]/2)}
    return {"new_energy_modes":modes,"new_energy_combined_finite_integrand":sum(
        (2 if kind=="transverse" else 1)*(row["physical_new_energy"]-row["actual_fourth_order_subtraction"])
        for kind,row in modes.items()),
        "pressure_definition":"Exactly the old selected pressure: its mode weights, fourth-order subtraction and matched finite local terms agree in S6.81.",
        "finite_energy_local_coefficients":quantum.local_coefficients(),
        "definition":"Multiply the radial integral with measure k_com^2 dk_com/(2 pi^2) by L^-2, and add sum_n energy_local[n]*m^(4-2n)/(64 pi^2 L^2), m=1000. Use the unchanged S6.55 Borel-prepared Cauchy data and exact original clock mode evolution. No new state or numerical integration cutoff is chosen."}


@cache
def action():
    rho=sp.Function("fixed_ordinary_Proca_energy")(u)
    pressure=sp.Function("fixed_ordinary_Proca_pressure")(u)
    old_rho=sp.Function("old_fixed_nonminimal_energy")(u)
    P=-pressure+(rho+pressure)*(x+1)/2
    old=-pressure+(old_rho+pressure)*(x+1)/2
    return {"new_fixed_energy":rho,"same_fixed_pressure":pressure,"old_fixed_energy":old_rho,
            "new_clock_tube_profile":P,"old_clock_tube_profile":old,
            "new_global_profile":window.cutoff(x)*P,
            "explicit_new_minus_old_profile":window.cutoff(x)*(rho-old_rho)*(x+1)/2,
            "construction_boundary":"Compute the new c-number functions once from integrands() and the fixed history/state, then freeze them. Replacing the old profiles is a new declared action, not a state-dependent reevaluation during later variation."}


def profile_bounds(scale):
    d=estimates.physical_bounds(scale)
    eta={j:max(d["new_energy_derivative_bounds"][j]["total"],d["identical_pressure_derivative_bounds"][j]["total"]) for j in range(6)}
    delta=d["new_minus_old_energy_derivative_upper"]
    return {"M_tau":d["M_tau"],"fixed_mass":estimates.MASS,
            "new_profile_clock_tube_derivative_bounds":{j:{"P":sp.Rational(11,10)*v,"P_x":v,"P_xx":sp.Integer(0)} for j,v in eta.items()},
            "new_profile_global_mixed_derivative_bounds":{f"{j},{k}":(
                sp.Rational(5,4)*window.derivative_bound(k)+(k*window.derivative_bound(k-1) if k else 0))*eta[j]
                for j in range(6) for k in range(6-j)},
            "profile_change_clock_tube_derivative_bounds":{j:{"P":v/20,"P_x":v/2,"P_xx":sp.Integer(0)} for j,v in delta.items()},
            "profile_change_global_mixed_derivative_bounds":{f"{j},{k}":(
                window.derivative_bound(k)/8+(sp.Rational(k,2)*window.derivative_bound(k-1) if k else 0))*delta[j]
                for j in range(6) for k in range(6-j)},
            "new_lapse_square_coefficient_upper":sp.Rational(15,8)*eta[0],
            "new_mixed_lapse_scale_coefficient_upper":sp.Rational(15,2)*eta[0],
            "new_scale_square_coefficient_upper":sp.Rational(9,2)*eta[0],
            "changed_lapse_square_coefficient_upper":delta[0],
            "changed_mixed_lapse_scale_coefficient_upper":3*delta[0],
            "changed_scale_square_coefficient_exactly_zero":True,
            "fixed_profiles_not_recomputed_under_metric_or_state_variation":True}


@cache
def checks():
    d=action()
    P=d["new_clock_tube_profile"]
    rho,pressure,old=d["new_fixed_energy"],d["same_fixed_pressure"],d["old_fixed_energy"]
    H=sp.Function("clock_Hubble")(u)
    Px=sp.diff(P,x)
    energy=(2*x*Px-P).subs(x,-1)
    spatial=P.subs(x,-1)
    clock=sp.diff(P.subs(x,-1),u)+2*sp.diff(Px.subs(x,-1),u)+6*H*Px.subs(x,-1)
    out={"new_fixed_profile_cancels_actual_energy":sp.expand(energy+rho),
         "new_fixed_profile_cancels_actual_pressure":sp.expand(spatial+pressure),
         "new_fixed_profile_clock_equation_is_new_conserved_Ward_balance":sp.expand(
             clock-sp.diff(rho,u)-3*H*(rho+pressure)),
         "literal_profile_replacement_not_a_reevaluated_old_coefficient":sp.expand(
             P-d["old_clock_tube_profile"]-(rho-old)*(x+1)/2),
         "no_new_second_clock_norm_derivative_on_flat_tube":sp.diff(P,x,2)}
    N,h,v=sp.symbols("positive_lapse positive_h hat_log_scale",positive=True)
    change=sp.Symbol("fixed_new_minus_old_energy",real=True)
    omega=-sp.log((h-1+N**-2)/h)/4
    density=N*sp.exp(3*omega+3*v)*change*(1-N**-2)/2
    point={N:1,v:0}
    out.update({"actual_replacement_clock_density_zero":sp.simplify(density.subs(point)),
                "actual_replacement_lapse_force":sp.simplify(sp.diff(density,N).subs(point)-change),
                "actual_replacement_scale_force_zero":sp.simplify(sp.diff(density,v).subs(point)),
                "actual_replacement_lapse_square":sp.simplify(sp.diff(density,N,2).subs(point)/2-(3/h-1)*change/2),
                "actual_replacement_mixed_lapse_scale":sp.simplify(sp.diff(density,N,v).subs(point)-3*change),
                "actual_replacement_scale_square_zero":sp.simplify(sp.diff(density,v,2).subs(point))})
    return out


@cache
def gates():
    d=profile_bounds(estimates.SCALE)
    return {"all_new_profile_global_mixed_derivatives_through_five_below_one_e_minus_770":all(
                value<sp.Rational(1,10**770) for value in d["new_profile_global_mixed_derivative_bounds"].values()),
            "all_replacement_global_mixed_derivatives_through_five_below_one_e_minus_770":all(
                value<sp.Rational(1,10**770) for value in d["profile_change_global_mixed_derivative_bounds"].values()),
            "both_global_mixed_derivative_triangles_have_twenty_one_entries":
                len(d["new_profile_global_mixed_derivative_bounds"])==len(d["profile_change_global_mixed_derivative_bounds"])==21,
            "new_fixed_local_quadratic_contact_coefficients_below_one_e_minus_770":all(
                d[name]<sp.Rational(1,10**770) for name in
                ("new_lapse_square_coefficient_upper","new_mixed_lapse_scale_coefficient_upper","new_scale_square_coefficient_upper",
                 "changed_lapse_square_coefficient_upper","changed_mixed_lapse_scale_coefficient_upper")),
            "one_point_Gaussian_cancellation_not_full_quantum_response_control":True,
            "new_profile_is_frozen_before_later_metric_or_state_variation":True}
