"""Actual prepared block inverse: local time constraint plus scalar L1 range kernel."""
from functools import cache

import sympy as sp
from p8_prepared_volterra import primitive as generic
from p8_proca_rank_one_inverse import kernel

from . import coordinates, tree


def scale(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact positive gravitational mass-time product")
    value=sp.Rational(value)
    if value<=0:
        raise ValueError("Require L>0 at fixed ordinary Proca mass")
    return value


def coupling(value):
    value=scale(value)
    return {"M_tau":value,"fixed_mass_time_product":sp.Integer(1000),
            "gamma":1/(64*sp.pi**2*value**2),
            "profile_policy":"This is a normalization formula; no source-dependent state or fixed profile is reselected."}


@cache
def data():
    L=sp.Symbol("positive_fixed_gravitational_scale",positive=True)
    gamma=1/(64*sp.pi**2*L**2)
    return {"gamma":gamma,"time_channel_instantaneous_coefficient":tree.data()["highest_time_channel_coefficient"],
            "time_channel_inverse":1/tree.data()["highest_time_channel_coefficient"],
            "time_channel_inverse_continuous_absolute_upper":tree.data()["instantaneous_time_channel_inverse_absolute_upper"],
            "prepared_diagonal_reference":"diag(A(u), gamma*F_m(partial_u^2)); A=-6*delta(u)^2",
            "prepared_diagonal_inverse":"diag(A(u)^-1, gamma^-1*K_m convolution); K_m has no instantaneous part and lies in L1.",
            "actual_four_primitive_normal_form":"I4(T+gamma*Q)=diag(A,gamma*F_m)+V; V has a weak-log Volterra kernel and all fixed diagonal derivatives.",
            "first_row_remainder_stronger_class":"The first row is entirely local before I4 and at most order three after removing A*eta''''; its I4 kernel and first output derivative are bounded.",
            "uniform_majorant":"v(r)=C*(1+abs(log(r))), 0<r<=1, C finite, not numerically evaluated.",
            "inverse_composition_majorant":"w(r)=Ainv_upper*v(r)+gamma^-1*(abs(K_m)*v)(r), in L1(0,1).",
            "weighted_norm_choice":"Choose finite lambda>0 with integral_0^1 exp(-lambda*r)*w(r)dr < 1/2.",
            "adapted_C0_estimate":"||(eta,w)|| <= 2*exp(lambda)*max(Ainv_upper,gamma^-1*||K_m||_L1)*||I4*g_adapted||.",
            "actual_clock_H_absolute_upper":sp.Rational(8,5),
            "I4_transformed_force_over_original_C0_upper":sp.Rational(13,30),
            "I3_transformed_time_force_over_original_C0_upper":sp.Rational(47,30),
            "physical_scale_over_adapted_C0_upper":sp.Rational(13,5),
            "physical_reconstruction_C0_estimate":"Differentiating the local first integrated row bounds eta' by Ainv_upper*(||I3*g_eta||+C1*||(eta,w)||), with finite C1 from A' and the bounded first-row kernel; n=eta', v=w+H*eta.",
            "smooth_prepared_domain":"All sources and forcing vanish on an initial neighborhood; same original selected state, fixed profile and zero matter-charge perturbation.",
            "no_arbitrary_initial_jet_or_covariance_selection":True,
            "no_numerical_inverse_norm_or_stability_claim":True}


@cache
def checks():
    out={"generic_four_primitive_"+key:value for key,value in generic.checks().items()}
    A,B=sp.symbols("nonzero_time_coefficient nonzero_scalar_block",nonzero=True)
    gamma=sp.Symbol("positive_quantum_factor",positive=True)
    out["actual_rank_one_loop_completed_by_classical_time_channel"]=sp.diag(A,gamma*B)*sp.diag(1/A,1/(gamma*B))-sp.eye(2)
    out["new_scalar_inverse_instantaneous_part_is_zero"]=kernel.data()["instantaneous_inverse"]
    out["new_time_inverse_continuous_bound"]=data()["time_channel_inverse_continuous_absolute_upper"]-sp.Rational(15625,6144)
    out["coupling_scale_law_fixed_mass"]=coupling(1)["gamma"]-4*coupling(2)["gamma"]
    out["time_constraint_has_no_quantum_fourth_derivative"]=sp.diff(
        coordinates.data()["fixed_profile_plus_vector_local_time_channel"],sp.diff(coordinates.eta,coordinates.u,4))
    out["time_constraint_has_no_quantum_third_scale_derivative"]=sp.diff(
        coordinates.data()["fixed_profile_plus_vector_local_time_channel"],sp.diff(coordinates.scale,coordinates.u,3))
    t,s=generic.t,generic.s
    for order in range(4):
        k=generic.local_kernel(order)
        # A polynomial in the lag remains an ordinary bounded kernel after one output derivative.
        direct=(-1)**order*sp.diff((t-s)**3*sp.Function("local_coefficient")(s)/6,s,order)
        out["physical_lapse_reconstruction_bounded_first_row_kernel_"+str(order)]=sp.simplify(sp.diff(k-direct,t))
    # I4[-gN']=-I3[gN] for zero-past gN, so no derivative of raw forcing is needed.
    out["prepared_force_time_transform_primitive_sign"]=sp.simplify(
        sp.diff((t-s)**3/6,s)+(t-s)**2/2)
    u=tree.u
    Hclock=4*u/(1+u*u)
    out["actual_clock_H_absolute_extrema"]=sp.ImmutableMatrix([
        Hclock.subs(u,sp.Rational(1,2))-sp.Rational(8,5),
        Hclock.subs(u,-sp.Rational(1,2))+sp.Rational(8,5)])
    out["actual_clock_H_increasing_on_fixed_interval"]=sp.factor(sp.diff(Hclock,u)-4*(1-u*u)/(1+u*u)**2)
    out["original_force_I4_continuous_primitive_bound"]=sp.Rational(1,6)+4*sp.Rational(8,5)/24-data()["I4_transformed_force_over_original_C0_upper"]
    out["original_force_I3_continuous_primitive_bound"]=sp.Rational(1,2)+4*sp.Rational(8,5)/6-data()["I3_transformed_time_force_over_original_C0_upper"]
    out["physical_scale_reconstruction_continuous_bound"]=1+sp.Rational(8,5)-data()["physical_scale_over_adapted_C0_upper"]
    return out


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("unfixed"),sp.sqrt(2),[],0,-1)
    count=0
    for value in bad:
        for function in (scale,coupling):
            try:
                function(value)
            except (TypeError,ValueError):
                count+=1
    if count!=2*len(bad):
        raise ValueError("An inexact or nonpositive prepared inverse scale was accepted")
    return {"rejected_inputs":count,"fixed_mass_not_varied_with_L":True,
            "no_parent_action_state_profile_or_library_monkeypatched":True}
