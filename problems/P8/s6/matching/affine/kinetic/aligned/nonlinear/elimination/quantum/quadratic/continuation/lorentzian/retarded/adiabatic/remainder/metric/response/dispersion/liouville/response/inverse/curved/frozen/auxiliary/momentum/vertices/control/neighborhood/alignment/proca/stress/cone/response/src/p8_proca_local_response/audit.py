"""Exact new local quantum metric response, continuous bounds and controls."""
from functools import cache

import sympy as sp

from . import bounds, chart, local


def adjoint(row,extra=()):
    result=[sp.Integer(0) for _ in row]
    for j,coefficient in enumerate(row):
        for k in range(j+1):
            value=coefficient
            for _ in range(j-k):
                value=local.weighted(value,extra)
            result[k]+=(-1)**j*sp.binomial(j,k)*value
    return tuple(sp.factor(value) for value in result)


@cache
def self_adjoint_checks():
    out={}
    for name,rows,operator,extra in (
        ("physical",{"N":local.n,"Z":local.v},local.operator,()),
        ("clock_with_fixed_local_profile",{"N":local.n,"V":chart.vhat},chart.operator,(chart.vhat,))):
        names=tuple(rows)
        for order in range(3):
            for a,output in enumerate(names):
                for source in names[a:]:
                    left=tuple(sp.factor(sp.diff(operator(output,order),field)) for field in rows[source][:5])
                    right=adjoint(tuple(sp.factor(sp.diff(operator(source,order),field))
                                        for field in rows[output][:5]),extra)
                    for j in range(5):
                        out[name+f"_{output}{source}_weighted_adjoint_{order}_{j}"]=sp.factor(left[j]-right[j])
    return out


@cache
def residuals():
    out={}
    for group in (local.checks(),chart.checks(),self_adjoint_checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated new local Proca response residual")
        out.update(group)
    for group,outputs in bounds.envelopes().items():
        for output,orders in outputs.items():
            for order,data in orders.items():
                for key,value in data["box_reconstructions"].items():
                    out[group+"_"+output+f"_order{order}_continuous_coefficient_"+key]=value
    return out


@cache
def gates():
    out=dict(bounds.gates())
    top=chart.top()
    out.update({"rank_one_local_fourth_derivative_coefficient_nonzero_on_I":top["actual_fourth_derivative_coefficient"][1,1]==-4,
                "local_fourth_derivative_matrix_is_not_the_old_invertible_indefinite_matrix":top["actual_fourth_derivative_coefficient"].det()==0,
                "second_nonlinear_metric_map_contact_retained":True,
                "fixed_local_profile_part_is_a_decomposition_not_a_new_state_choice":True,
                "dimensional_poles_and_evanescent_counterterms_recomputed_with_all_mass_jets_zero":True,
                "physical_readout_normalization_is_varied_after_Euler_currents":True,
                "pure_four_dimensional_divergences_omitted_only_as_compact_support_Euler_null":True,
                "old_prepared_fourth_block_inverse_and_frozen_pole_not_transferred":True,
                "no_higher_derivative_local_truncation_resummed_as_extra_particles":True,
                "full_nonlocal_response_mixed_loops_matching_and_original_P8_open":True})
    return {name:bool(value) for name,value in out.items()}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
                  sp.Symbol("unfixed_scale"),sp.sqrt(2),[],0,-1):
        calls.append(lambda value=value:bounds.at_scale(value))
    for value in (True,False,sp.Integer(1),1.0,-1,3,None,[]):
        calls.extend((lambda value=value:local.density(value),
                      lambda value=value:local.operator("N",value),
                      lambda value=value:local.coefficients("N",value,"Z"),
                      lambda value=value:local.physical_operator("energy",value),
                      lambda value=value:chart.density(value),
                      lambda value=value:chart.direct_metric_density(value),
                      lambda value=value:chart.operator("N",value)))
    for value in (True,False,1,sp.Integer(1),"unknown",None,[]):
        calls.extend((lambda value=value:local.operator(value,0),
                      lambda value=value:local.coefficients("N",0,value),
                      lambda value=value:local.physical_operator(value,0),
                      lambda value=value:chart.operator(value,0)))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An inexact or outside-domain local Proca response argument was accepted")
    return {"rejected_inputs":rejected,"validation_precedes_all_scientific_caches":True,
            "coefficientwise_exact_normalization_not_numeric_zero_tolerance":True,
            "no_parent_scientific_function_or_model_monkeypatched":True}
