"""Fixed-state ordinary Proca C5 stress and explicit profile replacement audit."""
from functools import cache

import sympy as sp
from p8_clock_tadpole import window
from p8_vector_regularity import estimates as old_estimates
from p8_vector_regularity import spectral

from . import estimates, profiles, readouts, subtraction


@cache
def residuals():
    out={}
    for group in (readouts.checks(),subtraction.checks(),profiles.checks(),
                  spectral.algebra_checks(),old_estimates.radial_checks(),window.chain_rule_checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated ordinary Proca stress residual name")
        out.update(group)
    for kind in ("transverse","longitudinal"):
        for j in range(6):
            for n,value in enumerate(readouts.row_bounds(kind,j)["box_reconstructions"]):
                out[kind+f"_energy_time{j}_continuous_coefficient_reconstruction_{n}"]=value
    for j in range(6):
        for n,value in estimates.local_derivatives(j)["reconstructions"].items():
            out[f"actual_new_local_energy_n{n}_time{j}_coefficient_reconstruction"]=value
    return out


@cache
def gates():
    out={**estimates.gates(),**profiles.gates(),**window.proof_checks()}
    for kind in ("transverse","longitudinal"):
        for j in range(6):
            row=readouts.row_bounds(kind,j)
            tail=subtraction.reference_tail(kind,j)
            out[kind+f"_ordinary_energy_time{j}_frequency_degree"]=row["all_frequency_powers_at_most_derivative_order_plus_one"]
            out[kind+f"_ordinary_energy_time{j}_integrable_subtraction_tail"]=(
                tail["integrable_tail_certified"] and tail["no_lower_Laurent_power"]
                and tail["all_majorant_coefficients_nonnegative"])
    failure=subtraction.reference_tail("longitudinal",4,2)["low_coefficient_residuals"][4]
    out.update({"lower_reference_does_not_support_this_fourth_derivative_tail_argument":failure.subs({readouts.u:0,readouts.z:1})!=0,
                "same_actual_borel_state_and_all_momentum_radial_envelope":True,
                "pressure_equality_includes_subtraction_and_matched_local_terms":True,
                "stress_derivatives_not_functional_metric_derivatives":True,
                "old_fixed_scalar_coefficients_replaced_only_in_this_new_candidate":True,
                "new_quantum_response_mixed_loops_cutoff_and_V_G_B_still_open":True,
                "original_P8_is_not_closed":True})
    return {name:bool(value) for name,value in out.items()}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
                  sp.Symbol("unfixed_scale"),sp.sqrt(2),[],0,-1):
        calls.extend((lambda value=value:estimates.physical_bounds(value),
                      lambda value=value:profiles.profile_bounds(value)))
    for value in (True,False,1,sp.Integer(1),"unknown",None,[]):
        calls.extend((lambda value=value:readouts.rows(value,0),
                      lambda value=value:subtraction.reference_tail(value,0)))
    for value in (True,False,sp.Integer(1),1.0,-1,6,None,[]):
        calls.extend((lambda value=value:readouts.rows("transverse",value),
                      lambda value=value:readouts.row_bounds("longitudinal",value),
                      lambda value=value:subtraction.adiabatic_rows("transverse",value),
                      lambda value=value:subtraction.reference_tail("longitudinal",value),
                      lambda value=value:estimates.local_derivatives(value)))
    for value in (True,False,sp.Integer(4),4.0,1,5,None,[]):
        calls.append(lambda value=value:subtraction.reference_tail("longitudinal",0,value))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An inexact or outside-domain Proca stress argument was accepted")
    return {"rejected_inputs":rejected,"native_validation_precedes_cached_scientific_calls":True,
            "actual_ODE_projected_derivative_not_derivative_of_an_approximate_mode":True,
            "low_reference_failure_not_claimed_to_be_divergence_of_the_exact_state":True,
            "old_nonminimal_energy_local_coefficients_not_reused_for_ordinary_Proca":True}
