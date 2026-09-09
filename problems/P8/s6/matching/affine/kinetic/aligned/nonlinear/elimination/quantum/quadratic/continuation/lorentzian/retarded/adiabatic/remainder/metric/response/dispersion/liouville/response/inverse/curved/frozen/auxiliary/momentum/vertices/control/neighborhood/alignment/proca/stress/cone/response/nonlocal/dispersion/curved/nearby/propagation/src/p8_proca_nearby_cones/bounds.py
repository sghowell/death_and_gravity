"""Whole-box exact polynomial enclosures and the actual-solution bridge."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import bounds as old_bounds
from p8_proca_nearby_bounce import taylor
from sympy.polys.domains import QQ

from . import rational

TIME_WINDOW=sp.Rational(1,10**7)
LAPSE_CENTER=sp.Rational(1000001,1000000)
LAPSE_RADIUS=sp.Rational(4,10**7)
TRACE_RADIUS=sp.Rational(3,10**5)
RADII=(TIME_WINDOW,LAPSE_RADIUS,TRACE_RADIUS)
CLOCK_SPEED_EXCESS_LOWER=sp.Rational(25,10**7)
CLOCK_SPEED_EXCESS_UPPER=sp.Rational(8,10**6)


def exact_duration(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact rational characteristic interval")
    value=sp.Rational(value)
    if not 0<value<=TIME_WINDOW:
        raise ValueError("Require a positive half-width at most 10^-7")
    return value


def _qq(value):
    return QQ(int(value.p),int(value.q))


def _sympy(value):
    return sp.Rational(int(value.numerator),int(value.denominator))


def polynomial_box(poly):
    if poly.ring != rational.FIELD.ring:
        raise ValueError("Polynomial belongs to a different exact box ring")
    shifted=poly.shift_list([QQ.zero,_qq(LAPSE_CENTER),QQ.zero])
    center=shifted.get((0,0,0),QQ.zero)
    radii=tuple(_qq(value) for value in RADII)
    deviation=sum(abs(value)*radii[0]**powers[0]*radii[1]**powers[1]*radii[2]**powers[2]
                  for powers,value in shifted.items() if powers!=(0,0,0))
    return {"lower":_sympy(center-deviation),"upper":_sympy(center+deviation),
            "center":_sympy(center),"deviation_upper":_sympy(deviation),
            "unshifted_monomials":len(poly),"shifted_monomials":len(shifted)}


def rational_box(value):
    if value.field != rational.FIELD:
        raise ValueError("Rational function belongs to a different exact box field")
    num,den=polynomial_box(value.numer),polynomial_box(value.denom)
    if den["lower"]*den["upper"]<=0:
        raise ValueError("The whole-box denominator enclosure meets zero")
    endpoints=[a/b for a in (num["lower"],num["upper"])
               for b in (den["lower"],den["upper"])]
    return {"lower":min(endpoints),"upper":max(endpoints),"numerator":num,"denominator":den}


@cache
def enclosures():
    s=rational.system()
    return {name:rational_box(s[name]) for name in (
        "clock_physical_speed_squared_excess","on_constraint_rescaled_lapse_Hessian",
        "constraint_matter_square","M")}


@cache
def solution_bridge():
    parent=taylor.evolution()
    rho=TIME_WINDOW/sp.Rational(1,100)
    b=rho/(1-rho*rho)
    e_upper=1+sp.Rational(1,10**6)
    nlo,nhi=LAPSE_CENTER-LAPSE_RADIUS,LAPSE_CENTER+LAPSE_RADIUS
    scale_radius=sp.Rational(2,10**9)
    Rlo,Rhi=1-scale_radius,1+scale_radius
    pivot=enclosures()["on_constraint_rescaled_lapse_Hessian"]
    M=enclosures()["M"]
    max_M=max(abs(M["lower"]),abs(M["upper"]))
    clock_kinetic=-pivot["upper"]/(4*e_upper**3*max_M**2)
    return {"parent_enlarged_time_window":taylor.TIME_WINDOW,
            "actual_phase_displacement_upper":parent["Picard_image_radius_upper"],
            "primitive_and_trace_linear_complex_time_radius":sp.Rational(1,100),
            "inherited_complex_trace_linear_upper":old_bounds.coefficients()["full_trace_linear_absolute_upper"],
            "odd_time_Cauchy_ratio":rho,"actual_trace_linear_absolute_upper":b,
            "rescaled_trace_absolute_upper":2*(parent["Picard_image_radius_upper"]+b),
            "declared_rescaled_trace_radius":TRACE_RADIUS,
            "actual_lapse_lower":nlo,"actual_lapse_upper":nhi,
            "actual_eomega_lower":sp.Integer(1),"actual_eomega_upper":e_upper,
            "actual_log_hat_scale_absolute_upper":sp.Rational(1,10**9),
            "actual_hat_scale_lower":Rlo,"actual_hat_scale_upper":Rhi,
            "physical_light_comoving_speed_lower":nlo/(e_upper*Rhi),
            "physical_light_comoving_speed_upper":nhi/Rlo,
            "clock_diagonal_kinetic_lower":clock_kinetic,
            "matter_diagonal_kinetic_lower":1/nhi,
            "tensor_principal_coefficient_lower":1/e_upper**4,
            "actual_scalar_and_tensor_gradients_positive":True,
            "sourcefree_ordinary_Proca_classical_positive_three_mode_parent_unchanged":True}


def characteristic_separation(value=TIME_WINDOW):
    value=exact_duration(value)
    return {"half_width":value,"clock_interval":(-value,value),
            "clock_physical_speed_excess_lower":CLOCK_SPEED_EXCESS_LOWER,
            "clock_physical_speed_excess_upper":CLOCK_SPEED_EXCESS_UPPER,
            "physical_light_comoving_speed_lower":sp.Rational(99,100),
            "physical_light_comoving_speed_upper":sp.Rational(101,100),
            "classical_comoving_characteristic_separation_lower":
                2*value*sp.Rational(99,100)*CLOCK_SPEED_EXCESS_LOWER,
            "classical_comoving_characteristic_separation_upper":
                2*value*sp.Rational(101,100)*CLOCK_SPEED_EXCESS_UPPER,
            "same_emission_event_and_same_final_clock_slice":True,
            "finite_frequency_wavepacket_error_or_cutoff_established":False,
            "not_a_UV_exclusion":True}


@cache
def checks():
    s=rational.system()
    c=s["clock_physical_speed_squared_excess"]
    ring=rational.FIELD.ring
    x,n,z=ring.gens
    toy=3+2*x+5*(n-_qq(LAPSE_CENTER))**2+7*x*z
    box=polynomial_box(toy)
    exact_deviation=2*TIME_WINDOW+5*LAPSE_RADIUS**2+7*TIME_WINDOW*TRACE_RADIUS
    speed=characteristic_separation()
    return {
        "whole_box_positive_and_negative_monomial_bound_anchor":box["deviation_upper"]-exact_deviation,
        "whole_box_shift_retains_exact_center":box["center"]-3,
        "uncancelled_rational_fraction_reproduces_speed_excess":
            sp.Integer(0) if c+1==s["clock_physical_speed_squared"] else sp.Integer(1),
        "actual_solution_box_uses_parent_phase_image":solution_bridge()["actual_phase_displacement_upper"]-sp.Rational(1,10**6),
        "characteristic_interval_includes_both_sides_of_bounce":
            speed["classical_comoving_characteristic_separation_lower"]-sp.Rational(99,2*10**14),
        "separation_scales_linearly_with_certified_interval":
            characteristic_separation(TIME_WINDOW/2)["classical_comoving_characteristic_separation_lower"]
            -speed["classical_comoving_characteristic_separation_lower"]/2}


@cache
def gates():
    e=enclosures()
    b=solution_bridge()
    cs=e["clock_physical_speed_squared_excess"]
    pv=e["on_constraint_rescaled_lapse_Hessian"]
    w=e["constraint_matter_square"]
    M=e["M"]
    return {key:bool(value) for key,value in {
        "entire_time_box_inside_actual_enlarged_solution":TIME_WINDOW==taylor.TIME_WINDOW,
        "odd_trace_linear_coefficient_has_complex_majorant_below_one":b["inherited_complex_trace_linear_upper"]<1,
        "actual_shifted_rescaled_trace_fits_box":b["rescaled_trace_absolute_upper"]<TRACE_RADIUS,
        "physical_lapse_strictly_above_one":b["actual_lapse_lower"]>1,
        "declared_eomega_upper_dominates_square_root_lapse":b["actual_eomega_upper"]**2>b["actual_lapse_upper"],
        "entire_box_clock_speed_squared_excess_above_five_e_minus_six":cs["lower"]>sp.Rational(5,10**6),
        "entire_box_clock_speed_squared_excess_below_sixteen_e_minus_six":cs["upper"]<sp.Rational(16,10**6),
        "strict_clock_speed_not_only_speed_squared_lower":cs["lower"]>(1+CLOCK_SPEED_EXCESS_LOWER)**2-1,
        "strict_clock_speed_upper":cs["upper"]<(1+CLOCK_SPEED_EXCESS_UPPER)**2-1,
        "whole_box_lapse_Hessian_rescaled_stays_negative":pv["upper"]<-2 and pv["lower"]>-3,
        "eliminated_positive_matter_square_is_physical":w["lower"]>sp.Rational(9,1000) and w["upper"]<sp.Rational(11,1000),
        "curvature_lapse_derivative_never_zero":M["upper"]<-sp.Rational(24,100) and M["lower"]>-sp.Rational(26,100),
        "clock_kinetic_Schur_complement_above_seven":b["clock_diagonal_kinetic_lower"]>7,
        "matter_kinetic_lower_above_ninety_nine_hundredths":b["matter_diagonal_kinetic_lower"]>sp.Rational(99,100),
        "tensor_principal_lower_above_ninety_nine_hundredths":b["tensor_principal_coefficient_lower"]>sp.Rational(99,100),
        "physical_light_coordinate_speed_lower_above_ninety_nine_hundredths":b["physical_light_comoving_speed_lower"]>sp.Rational(99,100),
        "physical_light_coordinate_speed_upper_below_one_point_zero_one":b["physical_light_comoving_speed_upper"]<sp.Rational(101,100),
        "finite_positive_characteristic_separation_lower":characteristic_separation()["classical_comoving_characteristic_separation_lower"]>sp.Rational(49,10**14),
        "no_all_frequency_Legendre_chart_claim":True,
        "no_quantum_state_profile_or_response_transfer":True,
        "finite_band_signal_cutoff_and_common_parent_UV_matching_remain_open":True}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.sqrt(2),[],0,-1,2*TIME_WINDOW)
    rejected=0
    for value in bad:
        for fn in (exact_duration,characteristic_separation):
            try:
                fn(value)
            except (TypeError,ValueError):
                rejected+=1
    F=rational.FIELD
    for value in (1/F.gens[0],1/(F.gens[1]-_qq(LAPSE_CENTER))):
        try:
            rational_box(value)
        except ValueError:
            rejected+=1
    if rejected!=2*len(bad)+2:
        raise ValueError("A bad interval or zero-crossing denominator was accepted")
    return {"rejected_inputs":rejected,"no_unverified_point_scan_or_floating_inequality":True}
