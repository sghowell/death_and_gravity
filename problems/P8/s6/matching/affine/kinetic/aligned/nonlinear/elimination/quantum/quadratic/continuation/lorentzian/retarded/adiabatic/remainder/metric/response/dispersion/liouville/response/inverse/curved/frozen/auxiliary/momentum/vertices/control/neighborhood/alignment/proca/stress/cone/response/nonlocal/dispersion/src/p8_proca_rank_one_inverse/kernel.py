"""Scalar range inverse: exact spectral moments and causal-kernel normalization."""
from functools import cache

import sympy as sp

from . import cut, spectral


def mass(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact positive fixed mass")
    value=sp.Rational(value)
    if value<=0:
        raise ValueError("Require a strictly positive fixed mass")
    return value


@cache
def data():
    return {"positive_density":cut.data()["positive_range_inverse_cut_density"],
            "inverse_spectral_representation":"R(p)=-integral_(4m^2)^infinity rho(tau)/(p+tau) dtau",
            "range_inverse_kernel":"K(t)=-2 theta(t) integral_(2m)^infinity rho(Omega^2) sin(Omega*t) dOmega",
            "instantaneous_inverse":sp.Integer(0),
            "threshold_density_constant":cut.data()["density_over_sqrt_z_at_threshold"],
            "high_density_log_constant":cut.data()["density_times_log_tau_over_mass_squared_squared_at_infinity"],
            "spectral_static_moment":sp.Rational(1,4),
            "absolutely_convergent_kernel_primitive":"J(t)=-integral_(4m^2)^infinity rho(tau)/tau*(1-cos(sqrt(tau)*t)) dtau",
            "primitive_lower":-sp.Rational(1,2),"primitive_upper":sp.Integer(0),
            "small_time_bound":"O(1/[t*log(1/t)^2]) at fixed m>0",
            "large_time_bound":"O(t^-3/2) at fixed m>0",
            "L1_half_line_kernel":True,"no_numerical_L1_norm_claim":True,
            "short_window_C0_inverse_norm_tends_to_zero":True,
            "same_mass_scaling":"K_m(t)=m K_1(m t); the full half-line L1 norm is mass independent for m>0",
            "scope":"Scalar inverse of the isolated normalized rank-one massive block only. The inactive lapse channel still has no inverse here."}


def moment_bound(value):
    value=mass(value)
    return {"mass":value,"instantaneous_C0_coefficient":sp.Integer(0),
            "primitive_C0_upper":sp.Rational(1,2),
            "static_positive_spectral_moment":sp.Rational(1,4),
            "inverse_squared_frequency_spectral_moment":sp.Rational(9,560)/value**2,
            "C1_zero_initial_bound":"For f(0)=0 on [0,T], ||Rop f||_C0 <= (T/2)||fprime||_C0; this is additional to the separate L1-based C0 endomorphism proof."}


@cache
def checks():
    s,omega,t,T=sp.symbols("positive_Laplace positive_frequency positive_time positive_final_time",positive=True)
    out={"new_causal_sine_Laplace_normalization":sp.factor(
            sp.integrate(sp.exp(-s*t)*(-2*sp.sin(omega*t)),(t,0,sp.oo))+2*omega/(s*s+omega*omega)),
         "new_kernel_primitive_normalization":sp.factor(
            sp.integrate(-2*sp.sin(omega*t),(t,0,T))+2*(1-sp.cos(omega*T))/omega),
         "new_finite_primitive_moment_bound":sp.Rational(1,2)-2*cut.data()["spectral_static_moment"],
         "new_zero_instantaneous_inverse":data()["instantaneous_inverse"],
         "new_fixed_mass_rescaling_of_next_moment":sp.factor(
             moment_bound(2)["inverse_squared_frequency_spectral_moment"]
             -moment_bound(1)["inverse_squared_frequency_spectral_moment"]/4)}
    A,L=sp.symbols("positive_spectral_scale positive_log",positive=True)
    leading=-1/(2*L-sp.Rational(14,15)+2*sp.I*sp.pi)
    out["new_no_hidden_high_frequency_instantaneous_term"]=sp.limit(leading,L,sp.oo)
    out["new_threshold_density_to_frequency_coordinate"]=sp.limit(
        (1-4*spectral.m**2/(2*spectral.m+A)**2)/A,A,0)-1/spectral.m
    return out

