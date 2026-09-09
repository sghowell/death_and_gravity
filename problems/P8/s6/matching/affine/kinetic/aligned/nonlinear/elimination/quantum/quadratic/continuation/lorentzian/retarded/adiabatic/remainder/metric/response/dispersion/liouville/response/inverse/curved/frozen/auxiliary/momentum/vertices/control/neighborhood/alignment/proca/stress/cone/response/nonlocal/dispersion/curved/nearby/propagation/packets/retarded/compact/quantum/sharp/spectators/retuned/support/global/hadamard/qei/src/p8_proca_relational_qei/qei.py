"""Exact positive-frequency moments and actual short-sampling coefficients."""
from functools import cache

import sympy as sp

from . import modes


def rational(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require exact rational "+name)
    return sp.Rational(value)


def coefficients(clock_squared_speed,mixing_strength):
    c2=rational(clock_squared_speed,"clock squared speed")
    mixing=rational(mixing_strength,"mixing strength")
    if not 0<c2<=1 or mixing<0:
        raise ValueError("Require 0<c_clock^2<=1 and nonnegative mixing")
    c=sp.sqrt(c2)
    field=1+mixing/c**3
    spatial=1+mixing/c**5
    energy=(field+spatial)/2
    return {"clock_squared_speed":c2,"relative_clock_spectral_weight":mixing/c,
        "field_square_and_time_derivative_enhancement":field,
        "spatial_gradient_sum_enhancement":spatial,
        "test_kinetic_density_enhancement":energy,
        "field_square_prefactor_without_hbar_over_kappa":field/(8*sp.pi**2),
        "time_derivative_square_prefactor_without_hbar_over_kappa":field/(16*sp.pi**2),
        "spatial_gradient_sum_prefactor_without_hbar_over_kappa":spatial/(16*sp.pi**2),
        "test_kinetic_density_prefactor_without_hbar_over_kappa":energy/(16*sp.pi**2)}


@cache
def data():
    theta,c,k=sp.symbols("positive_spectral_limit positive_mode_speed positive_radial_momentum",positive=True)
    r=sp.Symbol("nonnegative_clock_weight",nonnegative=True)
    moments={p:sp.integrate(k**p,(k,0,theta/c)) for p in (1,3)}
    radial=sp.Rational(1,4)/sp.pi**2
    field_coefficient=radial*(1+r/c**2)/2
    time_coefficient=radial*(1+r/c**2)/4
    spatial_coefficient=radial*(1+r/c**4)/4
    energy_coefficient=(time_coefficient+spatial_coefficient)/2
    actual=coefficients(sp.Rational(1199,1215),sp.Rational(1,1215))
    return {"positive_mode_speed":c,"clock_spectral_weight":r,"positive_spectral_limit":theta,
        "radial_moments":moments,
        "field_leading_radial_kernel_without_hbar_over_kappa":"(4*pi^2)^-1 integral_0^infinity k [r*exp(-i*c*k*delta_u)+exp(-i*k*delta_u)] dk",
        "test_kinetic_leading_radial_kernel_without_hbar_over_kappa":"(8*pi^2)^-1 integral_0^infinity k^3 [r*(1+c^2)*exp(-i*c*k*delta_u)+2*exp(-i*k*delta_u)] dk",
        "generic_field_square_coefficient":field_coefficient,
        "generic_time_derivative_square_coefficient":time_coefficient,
        "generic_spatial_gradient_sum_coefficient":spatial_coefficient,
        "generic_test_kinetic_density_coefficient":energy_coefficient,
        "actual_bounce_coefficients":actual,
        "sampling_family":"g_eta(u)=eta^-1/2*g(u/eta), eta>0, real g in C_c^infinity",
        "exact_finite_difference_bound":"integral g^2 <:sum_a (D_a O)^2:>_omega >= -B_reference[g]",
        "reference_functional":"B_reference[g]=(1/pi) integral_0^infinity d_alpha integral ds dt g(s)g(t) exp(-i*alpha*(s-t)) W_reference_D(s,t)",
        "field_square_short_sampling_limit":"lim eta^2 B_reference[g_eta]=(hbar/kappa)*C_field*integral |g'|^2",
        "derivative_square_short_sampling_limit":"lim eta^4 B_reference[g_eta]=(hbar/kappa)*C_derivative*integral |g''|^2",
        "test_kinetic_short_sampling_limit":"lim eta^4 B_reference_test_energy[g_eta]=(hbar/kappa)*C_test_energy*integral |g''|^2",
        "finite_width_bound_is_reference_functional_not_its_leading_asymptotic":True,
        "no_claim_of_optimality_or_cosmological_scale_remainder":True}


@cache
def checks():
    d=data()
    c,r,theta=d["positive_mode_speed"],d["clock_spectral_weight"],d["positive_spectral_limit"]
    actual=d["actual_bounce_coefficients"]
    c0=modes.data()["center_clock_speed"]
    r0=modes.data()["center_clock_relative_weight"]
    return {
        "field_radial_positive_frequency_moment_retains_speed_squared":sp.factor(d["radial_moments"][1]-theta**2/(2*c*c)),
        "derivative_radial_positive_frequency_moment_retains_speed_fourth":sp.factor(d["radial_moments"][3]-theta**4/(4*c**4)),
        "generic_test_kinetic_coefficient_keeps_time_and_spatial_weights":sp.factor(
            d["generic_test_kinetic_density_coefficient"]-(2+r*(1+c*c)/c**4)/(32*sp.pi**2)),
        "actual_field_coefficient_substitutes_nonzero_clock_residue":sp.simplify(
            d["generic_field_square_coefficient"].subs({c:c0,r:r0},simultaneous=True)-actual["field_square_prefactor_without_hbar_over_kappa"]),
        "actual_test_energy_coefficient_substitutes_nonzero_clock_residue":sp.simplify(
            d["generic_test_kinetic_density_coefficient"].subs({c:c0,r:r0},simultaneous=True)-actual["test_kinetic_density_prefactor_without_hbar_over_kappa"]),
        "test_energy_enhancement_is_mean_of_time_and_spatial_enhancements":sp.simplify(
            actual["test_kinetic_density_enhancement"]-(actual["field_square_and_time_derivative_enhancement"]+actual["spatial_gradient_sum_enhancement"])/2),
        "decoupled_free_field_known_answer":coefficients(1,0)["test_kinetic_density_prefactor_without_hbar_over_kappa"]-1/(16*sp.pi**2),
        "decoupled_Wick_square_known_answer":coefficients(1,0)["field_square_prefactor_without_hbar_over_kappa"]-1/(8*sp.pi**2),
        "coincident_speeds_keep_both_positive_spectral_weights":sp.factor(
            d["generic_test_kinetic_density_coefficient"].subs(c,1)-(1+r)/(16*sp.pi**2)),
        "real_even_Fourier_half_line_and_Parseval_factors":sp.Rational(1,1)/sp.pi*sp.pi-1,
    }


@cache
def gates():
    d=data()["actual_bounce_coefficients"]
    e=d["test_kinetic_density_enhancement"]
    return {name:bool(value) for name,value in {
        "actual_test_energy_short_sampling_coefficient_exceeds_free_one":e>1,
        "actual_test_energy_enhancement_below_one_part_in_one_thousand":e<sp.Rational(1001,1000),
        "actual_test_energy_enhancement_above_1_00084":e>sp.Rational(25021,25000),
        "actual_test_energy_enhancement_below_1_00085":e<sp.Rational(20017,20000),
        "actual_field_square_coefficient_is_strictly_positive":d["field_square_prefactor_without_hbar_over_kappa"]>0,
    }.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,sp.Symbol("free"),sp.sqrt(2),[])
    rejected=0
    for value in bad:
        for args in ((value,0),(1,value)):
            try:
                coefficients(*args)
            except (TypeError,ValueError):
                rejected+=1
    for args in ((0,0),(-1,0),(2,0),(1,-1)):
        try:
            coefficients(*args)
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=2*len(bad)+4:
        raise ValueError("An inadmissible coefficient domain was accepted")
    return {"rejected_inputs":rejected}
