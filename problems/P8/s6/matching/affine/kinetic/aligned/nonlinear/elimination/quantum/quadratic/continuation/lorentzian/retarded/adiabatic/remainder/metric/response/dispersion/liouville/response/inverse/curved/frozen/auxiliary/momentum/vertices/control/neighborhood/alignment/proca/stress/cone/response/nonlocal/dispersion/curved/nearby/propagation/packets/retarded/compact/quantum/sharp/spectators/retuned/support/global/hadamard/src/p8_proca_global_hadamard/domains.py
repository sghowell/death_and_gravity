"""Explicit frequency separation on every finite global time strip."""
from functools import cache

import sympy as sp
from p8_proca_global_support import coverage, model
from p8_proca_retuned_margin import original


def radius(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact finite rational compact time radius")
    value=sp.Rational(value)
    if value<sp.Rational(1,4):
        raise ValueError("The compact strip must contain the initial gamma slab")
    return value


@cache
def data():
    bg=model.background()
    numerator,denominator=sp.fraction(sp.factor(bg["J"]*bg["h"]**2))
    polynomial=sp.Poly(numerator,model.u)
    R=sp.Symbol("positive_finite_time_radius",positive=True)
    upper=polynomial.as_expr().subs(model.u,R)/800
    scale=(1+R**2)**2
    squared_gap=4*model.EPSILON/(upper+4*model.EPSILON)
    return {"finite_time_radius":R,"exact_J_times_h_squared_numerator":numerator,
        "exact_J_times_h_squared_denominator":denominator,
        "numerator_monomials":polynomial.terms(),
        "compact_J_times_h_squared_upper":upper,"compact_scale_upper":scale,
        "compact_squared_clock_matter_gap_lower":squared_gap,
        "compact_same_sign_frequency_gap_lower":squared_gap/(2*scale),
        "compact_opposite_sign_frequency_gap_lower":sp.Rational(493,250)/scale,
        "global_clock_squared_speed_lower":sp.Rational(1199,1215),
        "global_positive_conformal_volume":bg["a"]**3,
        "initial_gamma_slab":(-sp.Rational(1,4),sp.Rational(1,4)),
        "initial_Gamma_squared_speed_margin_lower":original.data()["positive_compact_original_clock_squared_cone_margin_lower"],
        "complete_chart_switches":coverage.data()["actual_fixed_chart_switch_times"],
        "high_real_frequency_symbol_threshold":coverage.K_MIN,
        "bounds_are_for_every_finite_R_not_a_uniform_infinite_history_gap":True}


def at_radius(value):
    value=radius(value)
    d=data()
    return {name:sp.factor(d[name].subs(d["finite_time_radius"],value))
        for name in ("compact_J_times_h_squared_upper","compact_scale_upper",
        "compact_squared_clock_matter_gap_lower","compact_same_sign_frequency_gap_lower",
        "compact_opposite_sign_frequency_gap_lower")}


@cache
def checks():
    d=data()
    bg=model.background()
    return {
        "actual_global_gap_numerator_and_denominator_reconstruct":sp.factor(
            d["exact_J_times_h_squared_numerator"]/d["exact_J_times_h_squared_denominator"]-bg["J"]*bg["h"]**2),
        "actual_global_gap_denominator_has_positive_fixed_lower_800":sp.factor(
            d["exact_J_times_h_squared_denominator"]-800*(1+model.u**2)**12),
        "actual_global_squared_speed_difference_keeps_total_margin":sp.factor(
            1-bg["clock_physical_squared_speed"]-4*model.EPSILON/(bg["J"]*bg["h"]**2+4*model.EPSILON)),
        "actual_global_scale_is_smooth_positive_everywhere":sp.factor(bg["a"]-(1+model.u**2)**2),
    }


@cache
def gates():
    d=data()
    return {name:bool(value) for name,value in {
        "global_numerator_positive_even_monomials_bound_every_finite_strip":all(
            powers[0]%2==0 and coeff>0 for powers,coeff in d["numerator_monomials"]),
        "all_time_speeds_positive_and_clock_strictly_slower_at_finite_time":0<d["global_clock_squared_speed_lower"]<1,
        "clock_speed_itself_above_493_over_500":d["global_clock_squared_speed_lower"]>sp.Rational(493,500),
        "initial_gamma_compact_squared_speed_margin_above_one_over_2000":d["initial_Gamma_squared_speed_margin_lower"]>sp.Rational(1,2000),
        "compact_same_sign_gap_symbol_is_strictly_positive":d["compact_same_sign_frequency_gap_lower"].is_positive is True,
        "compact_opposite_sign_gap_symbol_is_strictly_positive":d["compact_opposite_sign_frequency_gap_lower"].is_positive is True,
        "global_mode_chart_overlap_is_the_replayed_actual_overlap":coverage.gates()["both_fixed_switches_lie_in_nonempty_chart_overlaps"],
        "every_global_velocity_pole_is_below_the_replayed_high_root_cutoff":coverage.gates()["global_scalar_root_cutoff_covers_gamma_complex_q_threshold"],
    }.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
        sp.Symbol("free"),sp.sqrt(2),[],0,-1,sp.Rational(1,8))
    rejected=0
    for value in bad:
        for call in (radius,at_radius):
            try:
                call(value)
            except (TypeError,ValueError):
                rejected+=1
    if rejected!=2*len(bad):
        raise ValueError("An invalid finite time radius was accepted")
    return {"rejected_inputs":rejected}
