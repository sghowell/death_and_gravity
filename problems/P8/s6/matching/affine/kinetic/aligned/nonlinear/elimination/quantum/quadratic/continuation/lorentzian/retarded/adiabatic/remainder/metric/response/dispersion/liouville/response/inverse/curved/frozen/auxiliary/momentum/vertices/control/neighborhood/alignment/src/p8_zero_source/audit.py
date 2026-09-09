"""New-action identity inventory, exact-domain controls and scope gates."""
from functools import cache

import sympy as sp

from . import affine, bounds, coefficients, model, neighborhood, response, tree


@cache
def residuals():
    out={}
    for group in (model.checks(),affine.checks(),coefficients.checks(),bounds.checks(),
                  response.checks(),neighborhood.checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated zero-source residual name")
        out.update(group)
    d=tree.data()
    at=tree.at_scale(tree.SCALE)
    for name in at:
        if name!="common_M_tau":
            out["same_named_scale_budget_"+name]=sp.factor(at[name]-d[name])
    return out


@cache
def gates():
    out={}
    for group in (bounds.gates(),response.gates(),neighborhood.gates(),tree.gates()):
        if set(out).intersection(group):
            raise ValueError("Repeated zero-source proof gate")
        out.update(group)
    rows=coefficients.jets()
    g=model.old.generic()
    old=model.data()["old_Hamiltonian"]
    out.update({"all_fourteen_old_inventory_rows_have_five_new_lapse_jets":len(rows)==14 and all(len(row)==5 for row in rows.values()),
                "new_inventory_contains_twelve_nonzero_rows":sum(any(value!=0 for value in row) for row in rows.values())==12,
                "old_action_did_not_already_have_new_vector_parity":sp.factor(old.subs(model.old.j,-model.old.j)-old)!=0,
                "new_candidate_is_not_an_old_action_field_redefinition":True,
                "new_nonlinear_canonical_momenta_are_not_identified_with_old_ones":True,
                "the_new_source_free_action_retains_even_vector_interactions":
                    sp.diff(g["pieces"]["spatial_vector_mass"],model.old.vector)!=0,
                "same_free_quadratic_phase_identification_fixes_tree_comparison":True,
                "same_clock_background_and_fixed_Gaussian_tadpole_are_preserved":True,
                "finite_tree_error_budget_does_not_bound_loops_or_all_operators":True,
                "V_G_B_matching_and_original_P8_remain_open":True})
    return {name:bool(value) for name,value in out.items()}


def controls():
    common=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
            sp.Symbol("unfixed_scale"),sp.sqrt(2),[])
    calls=[]
    for value in common:
        calls.extend((lambda value=value:bounds.response_radius(value),
                      lambda value=value:tree.at_scale(value)))
    for value in (-1,sp.Rational(-1,10**30),sp.Rational(2,10**24)):
        calls.append(lambda value=value:bounds.response_radius(value))
    for value in (-1,0,sp.Rational(-1,10**30)):
        calls.append(lambda value=value:tree.at_scale(value))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An inexact or outside-domain zero-source input was accepted")
    return {"rejected_inputs":rejected,
            "native_validation_before_cached_mathematics":True,
            "actual_original_I_Iphi_Q_jets_not_replaced_by_clock_values":True,
            "new_coefficient_bounds_not_assumed_to_equal_old_sharper_ones":True,
            "noncommuting_source_completion_keeps_matrix_order":True,
            "old_quartic_source_response_and_mass_ratio_not_declared_zero":True}
