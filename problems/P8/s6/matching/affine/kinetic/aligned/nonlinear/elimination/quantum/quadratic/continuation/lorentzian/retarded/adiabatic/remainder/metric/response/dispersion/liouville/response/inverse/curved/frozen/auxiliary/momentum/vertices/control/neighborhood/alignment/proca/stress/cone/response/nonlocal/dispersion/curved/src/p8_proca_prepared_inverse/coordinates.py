"""Prepared physical time-coordinate variation and exact fixed-profile Ward reduction."""
from functools import cache

import sympy as sp
from p8_constant_proca import quantum
from p8_prepared_volterra import tree
from p8_proca_local_response import local
from p8_proca_nonlocal_response import source, tadpole
from p8_vector_metric_local import variation

u,H=tree.u,tree.H
eta,scale=(sp.Function(name)(u) for name in ("prepared_time_shift","prepared_invariant_logscale"))
rho,pressure=(sp.Function(name)(u) for name in ("fixed_vector_rho","fixed_vector_pressure"))
S=rho+pressure


def ward(value):
    value=sp.expand(value.doit())
    value=value.subs(sp.diff(rho,u,2),-3*sp.diff(H,u)*S-3*H*sp.diff(S,u))
    return sp.simplify(sp.expand(value).subs(sp.diff(rho,u),-3*H*S))


@cache
def data():
    ed,wd=sp.diff(eta,u),sp.diff(scale,u)
    eta_row=S*sp.diff(eta,u,2)+(sp.diff(S,u)+3*H*S)*ed-3*sp.diff(H,u)*S*eta-3*S*wd
    scalar_cross=3*sp.diff(pressure,u)*eta+3*S*ed
    return {"physical_lapse_source":ed,"physical_logscale_source":scale+H*eta,
            "normalized_physical_current_transform":"E_eta=-(partial_u+3H) E_N+H E_Z; E_w=E_Z",
            "fixed_profile_plus_vector_local_time_channel":eta_row,
            "fixed_profile_plus_vector_local_scalar_cross_channel":scalar_cross,
            "remaining_scalar_channel":"E_w=3*delta_pressure_Gamma[n=0,v=w]+local_scalar_cross_channel",
            "local_time_channel_maximum_derivative_order":2,
            "quantum_nonlocal_channel_count":1,
            "prepared_inverse_source_map":"eta(u)=integral_from_initial^u n(v)dv; w=v-H*eta",
            "profile_and_state_reselected":False}


@cache
def mode_checks():
    out={}
    time=source.time
    e=sp.Function("physical_prepared_time_shift")(source.u)
    n,v=sp.diff(e,source.u),source.H*e
    mapping={field:sp.diff(n,source.u,j) for j,field in enumerate(source.n)}
    mapping.update({field:sp.diff(v,source.u,j) for j,field in enumerate(source.v)})
    W=sp.Symbol("positive_frequency_squared",positive=True)
    Dx=lambda value:time(value)+2*source.lam*W*sp.diff(value,W)
    for sector in ("T","L"):
        actual=source.data(sector)
        d=actual["rate"]
        dr=sp.factor(actual["r"].subs(mapping,simultaneous=True))
        dd=sp.factor(actual["delta_d"].subs(mapping,simultaneous=True))
        out[sector+"_pure_time_change_frequency"]=sp.factor(dr-sp.diff(e,source.u)-source.lam*e)
        out[sector+"_pure_time_change_canonical_rate"]=sp.factor(dd-e*Dx(d)-d*sp.diff(e,source.u)+sp.diff(e,source.u,2)/2)
        M=sp.Matrix([[2*d,2,0],[-W,0,1],[0,-2*W,-2*d]])
        dM=sp.Matrix([[2*dd,0,0],[-2*dr*W,0,0],[0,-4*dr*W,-2*dd]])
        J=sp.diag(-1,0,1)
        transform=e*M+sp.diff(e,source.u)*J
        out[sector+"_exact_covariance_time_change"]=sp.ImmutableMatrix(
            (transform.applyfunc(Dx)+transform*M-M*transform-dM).applyfunc(sp.factor))
        for component in ("energy","pressure"):
            readout=source.readout(sector,component)
            row=sp.Matrix([[readout["B"]*W,0,readout["A"]]])
            varied=sp.Matrix([[(readout["delta_B"]+2*actual["r"]*readout["B"])*W,0,readout["delta_A"]]])
            varied=varied.subs(mapping,simultaneous=True)
            direct=varied+row*transform
            expected=e*(row.applyfunc(Dx)-3*source.H*row+row*M)
            out[sector+"_"+component+"_exact_time_change_scalar_readout"]=sp.ImmutableMatrix((direct-expected).applyfunc(sp.factor))
    return out


@cache
def subtraction_and_local_checks():
    out={}
    e=sp.Function("physical_prepared_time_shift")(source.u)
    mapping={field:sp.diff(e,source.u,j+1) for j,field in enumerate(source.n)}
    mapping.update({field:sp.diff(source.H*e,source.u,j) for j,field in enumerate(source.v)})
    for sector in ("T","L"):
        for component,label,factor in (("energy","N",-1),("pressure","Z",sp.Rational(1,3))):
            for order in range(3):
                baseline=factor*source.project(variation.data(sector)["baseline"][label][order])
                expected=e*(source.time(baseline)+(1-2*order)*source.lam*baseline-3*source.H*baseline)
                actual=source.adiabatic(sector,component,order).subs(mapping,simultaneous=True)
                out[sector+"_"+component+"_adiabatic_time_change_"+str(order)]=sp.factor(actual-expected)
    e=sp.Function("physical_prepared_time_shift")(local.u)
    mapping={field:sp.diff(e,local.u,j+1) for j,field in enumerate(local.n)}
    mapping.update({field:sp.diff(local.H*e,local.u,j) for j,field in enumerate(local.v)})
    for order in range(3):
        for component in ("energy","pressure"):
            baseline=quantum.local_coefficients()[order][component]
            actual=local.physical_operator(component,order).subs(mapping,simultaneous=True)
            out[component+"_actual_finite_local_time_change_"+str(order)]=sp.factor(actual-e*sp.diff(baseline,local.u))
        rho0,p0=(quantum.local_coefficients()[order][component] for component in ("energy","pressure"))
        out["actual_finite_local_background_Ward_"+str(order)]=sp.factor(sp.diff(rho0,local.u)+3*local.H*(rho0+p0))
    return out


@cache
def checks():
    out=mode_checks().copy()
    out.update(subtraction_and_local_checks())
    ed,wd=sp.diff(eta,u),sp.diff(scale,u)
    r,q=sp.Function("pure_scale_energy_response")(u),sp.Function("pure_scale_pressure_response")(u)
    total_r=r+sp.diff(rho,u)*eta+S*ed
    total_q=q+sp.diff(pressure,u)*eta+S*ed
    raw=sp.diff(total_r,u)+3*H*(total_r+total_q)
    response_ward={sp.diff(r,u):-3*H*(r+q)-3*S*wd}
    out["full_fixed_profile_time_channel_from_actual_varied_Ward"]=ward(raw.subs(response_ward)-data()["fixed_profile_plus_vector_local_time_channel"])
    out["full_fixed_profile_scalar_cross_from_actual_time_change"]=sp.expand(3*total_q-3*q-data()["fixed_profile_plus_vector_local_scalar_cross_channel"])
    out["local_cross_channels_weighted_adjoint"]=ward(
        3*sp.diff(S*eta,u)+9*H*S*eta-data()["fixed_profile_plus_vector_local_scalar_cross_channel"])
    local=-S*ed**2/2-sp.Rational(3,2)*sp.diff(H,u)*S*eta**2-3*S*eta*wd
    euler=lambda field:sp.diff(local,field)-sp.diff(sp.diff(local,sp.diff(field,u)),u)-3*H*sp.diff(local,sp.diff(field,u))
    out["local_time_profile_density_variation"]=ward(euler(eta)-data()["fixed_profile_plus_vector_local_time_channel"])
    out["local_scalar_profile_density_variation"]=ward(euler(scale)-data()["fixed_profile_plus_vector_local_scalar_cross_channel"])
    profile=tadpole.physical_vertices()
    out["actual_frozen_profile_normalized_stress_contact"]=sp.ImmutableMatrix(
        profile["fixed_profile_physical_stress_jacobian"]-sp.Matrix([[profile["rho"]+profile["pressure"],0],[profile["rho"]+profile["pressure"],0]]))
    EN,EV=sp.Function("original_lapse_residual")(u),sp.Function("original_scale_residual")(u)
    out["force_transform_has_only_first_order_prepared_kernel"]=sp.expand(
        (-(sp.diff(EN,u)+3*H*EN)+H*EV).subs(EV,0)+sp.diff(EN,u)+3*H*EN)
    out["prepared_linear_coordinate_map_roundtrip"]=sp.expand((scale+H*eta)-H*eta-scale)
    return out
