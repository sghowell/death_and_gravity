"""Positive scalar cut density, threshold gap and differentiable endpoints."""
from functools import cache

import sympy as sp

from . import spectral

z,m=spectral.z,spectral.m
ell=sp.Symbol("cut_atanh_sqrt_fraction",real=True)


@cache
def data():
    P=spectral.data()["positive_pair_polynomial"]
    D=sp.Rational(16,15)+z-3*z*z+sp.sqrt(z)*P*ell
    U=sp.sqrt(z)*P/2
    density=U/(D*D+sp.pi**2*U*U)
    lower=sp.Rational(16,15)+z*(3*(z-sp.Rational(5,6))**2+sp.Rational(23,12))
    return {"real_positive_H_on_upper_cut":D,"positive_imaginary_H_over_pi":U,
            "positive_range_inverse_cut_density":density,
            "cut_substitution":"ell=atanh(sqrt(z)), z=1-4m^2/tau, 0<z<1",
            "real_cut_H_continuous_lower":lower,"uniform_first_sheet_H_real_part_lower":sp.Rational(16,15),
            "range_inverse_absolute_upper_on_first_sheet":sp.Rational(15,16),
            "density_over_sqrt_z_at_threshold":sp.Rational(675,512),
            "density_times_log_tau_over_mass_squared_squared_at_infinity":sp.Rational(1,2),
            "spectral_static_moment":sp.Rational(1,4),
            "spectral_inverse_squared_frequency_moment":sp.Rational(9,560)/m**2,
            "no_instantaneous_range_inverse":True,
            "full_two_source_inverse":False}


def fraction(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact rational open-cut fraction")
    value=sp.Rational(value)
    if not 0<value<1:
        raise ValueError("Require 0<z<1 on the open massive cut")
    item=data()
    mapping={z:value,ell:sp.atanh(sp.sqrt(value))}
    return {key:sp.factor(item[key].subs(mapping,simultaneous=True)) for key in
            ("real_positive_H_on_upper_cut","positive_imaginary_H_over_pi","positive_range_inverse_cut_density","real_cut_H_continuous_lower")}


@cache
def checks():
    item=data()
    D,U=item["real_positive_H_on_upper_cut"],item["positive_imaginary_H_over_pi"]
    actual=sp.factor(sp.im(-1/(D+sp.I*sp.pi*U))/sp.pi)
    out={"new_scalar_inverse_cut_jump":sp.factor(actual-item["positive_range_inverse_cut_density"]),
         "strict_real_cut_lower_square":sp.factor(item["real_cut_H_continuous_lower"]
             -(sp.Rational(16,15)+4*z-5*z*z+3*z**3)),
         "cut_lower_remainder_exact":sp.factor(
             D-item["real_cut_H_continuous_lower"]-sp.sqrt(z)*spectral.data()["positive_pair_polynomial"]*(ell-sp.sqrt(z))),
         "positive_atanh_remainder_derivative":sp.simplify(
             sp.diff(sp.atanh(sp.sqrt(z))-sp.sqrt(z),z)-sp.sqrt(z)/(2*(1-z))),
         "new_threshold_gap_matches_radial_endpoint":D.subs({z:0,ell:0})-spectral.data()["H_at_threshold"],
         "new_threshold_density_constant":sp.factor(
             (item["positive_range_inverse_cut_density"]/sp.sqrt(z)).subs({z:0,ell:0})-sp.Rational(675,512)),
         "new_static_spectral_moment":sp.Rational(1,4)-1/spectral.data()["H_at_zero"],
         "new_next_spectral_moment":sp.factor(item["spectral_inverse_squared_frequency_moment"]
             -spectral.data()["H_slope_at_zero"]/spectral.data()["H_at_zero"]**2)}
    L=sp.Symbol("large_log_tau_over_mass_squared",positive=True)
    leading=(2)/((2*L-sp.Rational(14,15))**2+4*sp.pi**2)
    out["new_high_density_log_coefficient"]=sp.limit(L**2*leading,L,sp.oo)-sp.Rational(1,2)
    out["new_inverse_high_frequency_limit"]=sp.limit(-1/(2*L-sp.Rational(14,15)+2*sp.I*sp.pi),L,sp.oo)
    return out

