"""Endpoint constants and exact finite moment for the causal inverse."""
from functools import cache

import sympy as sp

from . import block, cut

h,z = block.h,block.z


@cache
def data():
    item = block.data()
    T,M,Q,E = (item[key] for key in ("threshold_matrix","spectral_matrix","chart","infinity_chart_constant"))
    threshold = block.clean(T.inv()*M.subs(z,0)*T.inv()/8)
    vector = block.clean(Q*sp.Matrix([-E[0,1]/E[0,0],1]))
    high = block.clean(vector*vector.T/2)
    return {"density_over_sqrt_z_at_threshold": threshold,
            "density_times_log_tau_over_mass2_squared_at_infinity": high,
            "high_frequency_density_vector": vector,
            "instantaneous_inverse_matrix": item["inverse_infinity"],
            "positive_static_spectral_moment": item["static_spectral_moment"],
            "kernel": "R_infinity*delta(t)-2*theta(t)*integral_(2m)^infinity rho(Omega^2)*sin(Omega*t) dOmega",
            "spectral_representation": "R(p)=R_infinity-integral_(4m^2)^infinity rho(tau)/(p+tau) dtau",
            "regular_kernel_small_time": "O(1/[t*log(1/t)^2]) at fixed positive m",
            "regular_kernel_large_time": "O(t^(-3/2)) at fixed positive m",
            "integrated_kernel": "-integral_(4m^2)^infinity rho(tau)/tau*(1-cos(sqrt(tau)*t)) dtau"}


def moment_bound(value):
    value = block.positive_exact(value)
    item = data()
    moment = item["positive_static_spectral_moment"].subs(h,value)
    infinity = item["instantaneous_inverse_matrix"].subs(h,value)
    # PSD matrices: spectral norm <= trace. Loewner -2M <= J(t) <= 0.
    return {"h":value,"instantaneous_spectral_norm_upper":sp.trace(infinity),
            "primitive_spectral_norm_upper":2*sp.trace(moment),
            "domain": "For f(0)=0, ||Rop f||_C0 <= trace(Rinf)||f||_C0 + 2*trace(M)*T*||fprime||_C0. The L1 proof also gives a C0 endomorphism, but no numerical L1 norm is supplied."}


@cache
def checks():
    item,parent = data(),block.data()
    out = {"threshold_density_constant_by_exact_inverse": block.clean(
        item["density_over_sqrt_z_at_threshold"]-sp.Matrix([
            [sp.Rational(67001424075,34630803488)*h**2,-sp.Rational(154675575,129219416)*h],
            [-sp.Rational(154675575,129219416)*h,sp.Rational(357075,482162)]])),
           "threshold_density_rank_one":sp.factor(item["density_over_sqrt_z_at_threshold"].det()),
           "high_density_constant_by_exact_schur_complement":block.clean(
        item["density_times_log_tau_over_mass2_squared_at_infinity"]-sp.Matrix([
            [sp.Rational(149597361,30482432)*h**2,-sp.Rational(4831245,1905152)*h],
            [-sp.Rational(4831245,1905152)*h,sp.Rational(156025,119072)]]))}
    Q,E = parent["chart"],parent["infinity_chart_constant"]
    L=sp.Symbol("large_log_tau_over_m2",positive=True)
    B=sp.Matrix([[E[0,0],E[0,1]],[E[0,1],E[1,1]-2*L-2*sp.I*sp.pi]])
    R=Q*B.inv()*Q.T
    out["leading_cut_inverse_limit"] = block.clean(
        R.applyfunc(lambda entry:sp.limit(entry,L,sp.oo))-parent["inverse_infinity"])
    rho=R.applyfunc(lambda entry:sp.factor(sp.im(entry)/sp.pi))
    out["leading_cut_density_limit"] = block.clean(
        rho.applyfunc(lambda entry:sp.limit(L**2*entry,L,sp.oo))
        -item["density_times_log_tau_over_mass2_squared_at_infinity"])
    tau,s,omega,t,mass = sp.symbols("tau laplace omega time mass",positive=True)
    out["spectral_variable_change"] = sp.factor(2*omega/(s**2+omega**2)-1/(s**2+tau)*2*omega).subs(tau,omega**2)
    out["causal_sine_Laplace_normalization"] = sp.factor(
        sp.integrate(sp.exp(-s*t)*(-2*sp.sin(omega*t)),(t,0,sp.oo))+2*omega/(s**2+omega**2))
    out["primitive_moment_normalization"] = sp.factor(
        sp.integrate(-2*sp.sin(omega*t),(t,0,sp.Symbol("T",positive=True)))
        +2*(1-sp.cos(omega*sp.Symbol("T",positive=True)))/omega)
    out["mass_scaling_of_threshold_coordinate"] = sp.factor(
        (1-4*mass**2/(mass*omega)**2)-(1-4/omega**2))
    out["finite_trace_moment_bound_is_exact"] = sp.factor(
        moment_bound(1)["primitive_spectral_norm_upper"]
        -2*sp.trace(parent["static_spectral_moment"].subs(h,1)))
    out["cut_endpoint_matrix_unchanged"] = block.clean(cut.data()["real_cut_matrix"].subs({z:0,cut.ell:0})-parent["threshold_matrix"])
    return out
