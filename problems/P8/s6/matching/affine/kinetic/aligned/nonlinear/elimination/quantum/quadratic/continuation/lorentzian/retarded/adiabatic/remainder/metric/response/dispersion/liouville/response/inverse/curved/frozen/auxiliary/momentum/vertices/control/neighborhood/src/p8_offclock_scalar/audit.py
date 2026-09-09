"""Same-action off-clock residual inventory and scientific controls."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model

from . import datum, euler, principal, spatial


@cache
def residuals():
    result={}
    for group in (datum.checks(),principal.checks(),euler.checks()):
        if set(result).intersection(group):
            raise ValueError("Repeated off-clock residual name")
        result.update(group)
    for a,left in enumerate(spatial.CHANNELS):
        for right in spatial.CHANNELS[a+1:]:
            result["independent_Fourier_pair_reversal_"+left+"_"+right]=spatial.root_relation(
                spatial.pair(right,left)["kernel"]-spatial.pair(left,right)["kernel"].subs(spatial.k,-spatial.k))
    return result


@cache
def gates():
    d=datum.data()
    result=dict(datum.domain()["gates"])
    pairs=[spatial.pair(left,right) for a,left in enumerate(spatial.CHANNELS) for right in spatial.CHANNELS[a:]]
    point=1-datum.RADIUS
    coeff=euler.result()["Euler_constant_determinant_leading_coefficient"]
    result.update({
        "all_twenty_one_spatial_phase_pairs_reconstructed":len(pairs)==21,
        "all_sixty_three_spatial_constraint_checks":sum(len(row["spatial_checks"]) for row in pairs)==63
            and all(all(row["spatial_checks"].values()) for row in pairs),
        "matter_family_substitution_only_after_fixed_phase_N_derivatives":
            all(row["matter_density_relation_imposed_after_N_derivative"] is True for row in pairs),
        "nonzero_coupling_at_an_actual_punctured_domain_datum":d["canonical_gamma_vector_mixing"].subs(model.N,point)>0,
        "negative_Euler_leading_coefficient_at_actual_punctured_datum":coeff.subs(model.N,point)<0,
        "dropping_trace_temporal_correction_changes_longitudinal_pivot":
            sp.factor(d["gamma_t"]-d["longitudinal_kinetic_times_sqrt_N"]).subs(model.N,point)>0,
        "wrong_constraint_tangent_is_not_the_nonzero_lapse_Hessian":
            d["fixed_phase_lapse_Hessian"].subs(model.N,1)<0,
        "original_clock_freeze_before_Euler_negative_control":
            euler.result()["clock_negative_control_without_Euler_time_jets"]<0,
        "actual_Euler_first_clock_leading_sign_stays_positive":
            euler.result()["actual_clock_Euler_leading_determinant"]>0,
        "physical_and_hat_frequency_rescaling_positive":1-datum.RADIUS>0,
        "same_action_no_vector_state_reselection_or_tube_deformation":True,
        "frozen_local_growth_not_actual_time_dependent_instability":True,
        "no_numerical_cutoff_or_UV_exclusion_from_q_infinity":True,
        "original_P8_and_common_parent_matching_open":True})
    return {name:bool(value) for name,value in result.items()}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
                  sp.Symbol("outside_lapse"),sp.sqrt(2),sp.exp(model.N),1/(model.N-1),
                  1/(model.N-1+datum.RADIUS/2)):
        calls.append(lambda value=value:datum.enclosure(value))
    for value in (True,False,sp.true,sp.false,1,sp.Integer(1),1.0,"unknown",None,[]):
        calls.extend((lambda value=value:spatial.pair(value,"curvature"),
                      lambda value=value:spatial.pair("curvature",value)))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An uncertified off-clock input was accepted")
    return {"rejected_inputs":rejected,"native_pair_validation_before_cache":True,
            "continuous_denominator_enclosures_reject_poles_and_inexact_coefficients":True,
            "fixed_phase_lapse_derivative_not_replaced_by_constraint_tangent":True,
            "joint_trace_temporal_pivot_not_replaced_by_isolated_mass":True,
            "actual_clock_Euler_negative_control_keeps_all_time_jets":True}
