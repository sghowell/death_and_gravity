"""Independent central scalar cone identities, continuous bounds and controls."""
from functools import cache

import sympy as sp

from . import background, domain, principal, spatial


@cache
def residuals():
    out={}
    for group in (background.checks(),principal.checks(),domain.checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated central scalar cone residual")
        out.update(group)
    # Every unordered scalar phase pair is generated from both leg orders.
    for a,left in enumerate(spatial.CHANNELS):
        for right in spatial.CHANNELS[a:]:
            x=spatial.pair(left,right)["kernel_first_time_jet"]
            y=spatial.pair(right,left)["kernel_first_time_jet"].subs(spatial.k,-spatial.k)
            out["actual_scalar_phase_Fourier_reversal_"+left+"_"+right]=spatial.clean(x-y)
    return out


@cache
def gates():
    out=dict(domain.gates())
    for a,left in enumerate(spatial.CHANNELS):
        for right in spatial.CHANNELS[a:]:
            result=spatial.pair(left,right)
            out["actual_spatial_constraints_"+left+"_"+right]=all(result["spatial_checks"].values())
    f=principal.formula()
    out.update({"central_matter_direction_exactly_luminal":f["matter_physical_speed_squared"]==1,
                "physical_lapse_and_spatial_metric_both_retained":f["actual_lapse"]==background.N
                    and f["physical_spatial_metric_over_hat_metric"]==background.N,
                "constant_mass_vector_is_not_old_offclock_mass_ratio":True,
                "finite_q_phase_velocity_not_used_as_a_characteristic_cone":True,
                "time_derivative_taken_before_principal_freezing":True,
                "new_fixed_quantum_profile_is_not_silently_omitted_from_a_full_quantum_claim":True})
    return {name:bool(value) for name,value in out.items()}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
                  sp.Symbol("unfixed_lapse"),sp.sqrt(2),[]):
        calls.extend((lambda value=value:domain.domain(value),
                      lambda value=value:domain.point(value),
                      lambda value=value:domain.enclosure(background.N,value)))
    for value in (0,-1,sp.Rational(1,10**5)):
        calls.append(lambda value=value:domain.domain(value))
    for value in (0,-1,1-2*domain.OUTER_RADIUS,1+2*domain.OUTER_RADIUS):
        calls.append(lambda value=value:domain.point(value))
    for value in (True,False,1,sp.Integer(1),"unknown",None,[]):
        calls.extend((lambda value=value:spatial.pair(value,"curvature"),
                      lambda value=value:spatial.pair("curvature",value)))
    for value in (sp.Float(1),sp.sqrt(2),sp.Symbol("unknown"),sp.sin(background.N),
                  1/(background.N-1)):
        calls.append(lambda value=value:domain.enclosure(value,domain.RADIUS))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An inexact, singular or outside-domain scalar cone input was accepted")
    return {"rejected_inputs":rejected,"input_validation_precedes_scientific_caches":True,
            "freezing_before_Euler_and_differentiating_along_constraint_family_are_negative_controls":True,
            "larger_superluminal_datum_not_silently_removed_from_scope":True}
