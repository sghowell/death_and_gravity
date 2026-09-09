"""All-time two-chart coverage and explicit complex-frequency pole separation."""
from functools import cache
from itertools import pairwise

import sympy as sp
from p8_proca_retuned_margin import original

from . import model

CUT=sp.Rational(3,16)
K_MIN=sp.Integer(64)


def exact(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require exact finite rational global support endpoints")
    return sp.Rational(value)


def time_pair(start,end):
    start,end=exact(start),exact(end)
    if start>end:
        raise ValueError("Global retarded endpoints must be ordered")
    nodes=[start]+[c for c in (-CUT,CUT) if start<c<end]+[end]
    segments=[]
    for left,right in pairwise(nodes):
        mid=(left+right)/2
        segments.append({"start":left,"end":right,"chart":"gamma" if -CUT<=mid<=CUT else "unitary"})
    primitive=model.background()["physical_matter_null_primitive"]
    return {"source_time":start,"output_time":end,
        "finite_compact_time_radius":max(sp.Rational(1,4),abs(start),abs(end)),
        "ordered_chart_segments":segments,"maximum_required_segments":3,
        "actual_physical_matter_null_radius":primitive.subs(model.u,end)-primitive.subs(model.u,start),
        "support_radius_is_actual_integral_not_a_segmentwise_overestimate":True}


@cache
def data():
    radius=sp.Symbol("finite_compact_time_radius",positive=True)
    old_bounds=original.data()
    J0upper=36+sp.Rational(1,50)+sp.Rational(1,800)
    return {"gamma_regular_time_domain":(-sp.Rational(1,4),sp.Rational(1,4)),
        "unitary_regular_time_domain":"abs(u)>=1/8 on each finite strip abs(u)<=R, R>=1/4",
        "actual_fixed_chart_switch_times":(-CUT,CUT),
        "gamma_absolute_Lambda_lower":sp.Rational(1,5),
        "gamma_absolute_Lambda_upper":sp.Rational(1,2),
        "global_absolute_matter_density_upper":sp.Rational(1,10),
        "gamma_half_matter_mixing_squared_upper":sp.Rational(1,800),
        "gamma_scale_squared_upper":(sp.Rational(17,16))**4,
        "unitary_absolute_Theta_lower_on_finite_strip":3/(8*(1+radius**2)),
        "global_J_lower_on_finite_strip":sp.Rational(1199,800)/(1+radius**2)**6,
        "new_global_clock_squared_speed_lower":old_bounds["positive_all_time_clock_squared_speed_lower"],
        "gamma_J_plus_margin_plus_half_w_squared_upper":J0upper,
        "gamma_complex_q_threshold_strict_upper":50*J0upper,
        "declared_gamma_complex_q_threshold":sp.Integer(2000),
        "declared_global_complex_scalar_root_threshold":K_MIN,
        "gamma_finite_q_pole_definition":"D=q*Lambda^2-J-delta_J-w^2/2",
        "gamma_high_complex_pole_separation":"abs(D)>=abs(q)*Lambda^2/2 when abs(q)>=2000",
        "all_time_scale_lower":sp.Integer(1),
        "complete_global_scalar_multiplier_type":"C_R*(1+norm(zeta))^9*exp(S_m(t,s)*norm(Im(zeta)))",
        "complete_global_transfer_frequency_degree":9,
        "whole_global_clock_conformal_time_length":sp.pi/2,
        "global_result_is_not_global_continuation_of_off_clock_nearby_solution":True}


@cache
def checks():
    d=data()
    o=model.old
    center={symbol:value.subs(model.u,0) for symbol,value in model.substitution().items()}
    D=o.q*o.lam**2-o.J-model.margin.delta_J-o.w**2/2
    pole=sp.Rational(152,25)
    return {"gamma_actual_center_velocity_chart_pole_is_152_over_25":sp.factor(D.subs(center).subs(o.q,pole)),
        "gamma_center_pole_matches_retuned_total_margin_formula":pole-6-16*model.EPSILON,
        "explicit_global_high_complex_threshold_square_is_4096":K_MIN**2-4096,
        "at_most_three_canonical_segment_bounds_give_frequency_degree_nine":3*(2+1)-d["complete_global_transfer_frequency_degree"],
        "global_null_primitive_limits_have_total_pi_over_two":
            sp.limit(model.background()["physical_matter_null_primitive"],model.u,sp.oo)
            -sp.limit(model.background()["physical_matter_null_primitive"],model.u,-sp.oo)-sp.pi/2,
        "global_coincident_time_radius_zero":time_pair(0,0)["actual_physical_matter_null_radius"]}


@cache
def gates():
    d=data()
    rows={"gamma_strip_h_is_below_five_quarters":sp.Rational(17,16)**3<sp.Rational(5,4),
        "gamma_strip_scale_squared_below_two":d["gamma_scale_squared_upper"]<2,
        "gamma_pole_threshold_strictly_below_two_thousand":d["gamma_complex_q_threshold_strict_upper"]<2000,
        "gamma_J_bound_follows_from_existing_larger_compact_interval":original.data()["compact_J_h_squared_upper"]<36,
        "global_scalar_root_cutoff_covers_gamma_complex_q_threshold":K_MIN**2/2>2000,
        "both_fixed_switches_lie_in_nonempty_chart_overlaps":sp.Rational(1,8)<CUT<sp.Rational(1,4),
        "global_clock_squared_speed_has_positive_uniform_lower":0<d["new_global_clock_squared_speed_lower"]<1,
        "whole_global_clock_kinetic_pivot_positive_even_polynomial_proof_replayed":
            original.margin_bounds.polynomial_bounds()["global_positive_even_remainder"],
        "global_reference_time_is_covered_by_gamma_chart":time_pair(0,0)["ordered_chart_segments"][0]["chart"]=="gamma",
        "large_global_interval_uses_no_more_than_three_regular_segments":len(time_pair(-10**6,10**6)["ordered_chart_segments"])==3}
    return {name:bool(value) for name,value in rows.items()}
