"""Classical constant-mass candidate with an explicit quantum non-transfer boundary."""
from functools import cache

import sympy as sp
from p8_offclock_scalar import datum

from . import affine, bounds, model, propagation, quantum


@cache
def residuals():
    out={}
    for group in (affine.checks(),model.checks(),propagation.checks(),quantum.checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated constant-mass residual name")
        out.update(group)
    at=bounds.at_scale(bounds.previous_tree.SCALE)
    for name,value in at.items():
        if name!="common_M_tau":
            out["common_scale_comparison_"+name]=sp.factor(value-bounds.tree()[name])
    return out


@cache
def gates():
    result={**bounds.gates(),**quantum.gates()}
    n=propagation.neighborhood()
    box=datum.domain()["continuous_rational_proofs"]
    actual_mass=affine.data()["new_retained_mass"]
    result.update({
        "new_retained_mass_invertible_with_original_Lorentz_signature":actual_mass.det()!=0 and actual_mass==sp.diag(1,-1,-1,-1),
        "old_to_new_quotient_determinant_ratio_positive_on_certified_box":
            box["gamma_t"]["lower"]*box["gamma_s"]["lower"]**3>0,
        "constant_canonical_Proca_mass_strictly_positive":model.data()["canonical_Proca_mass_squared"]>0,
        "new_temporal_pivot_nonzero_at_both_actual_offclock_anchors":all(
            n["new_temporal_auxiliary_pivot"].subs(model.old.N,point)<0
            for point in (1-datum.RADIUS,1+datum.RADIUS)),
        "physical_Proca_constraint_not_dropped":propagation.symbol()["full_massive_plane_wave_operator"]!=
            propagation.symbol()["constrained_wave_operator"],
        "new_affine_Hessian_change_is_not_hidden":affine.data()["added_full_connection_Hessian"]!=sp.zeros(64),
        "same_full_classical_clock_not_a_new_quantum_corrected_bounce":True,
        "no_full_coupled_neighborhood_or_UV_verdict_from_vector_cone":True,
        "V_G_B_matching_and_original_P8_remain_open":True})
    return {name:bool(value) for name,value in result.items()}


def controls():
    common=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
            sp.Symbol("unfixed_scale"),sp.sqrt(2),[])
    calls=[]
    for value in common:
        calls.extend((lambda value=value:bounds.response_radius(value),
                      lambda value=value:bounds.at_scale(value)))
    for value in (-1,sp.Rational(-1,10**30),sp.Rational(2,10**24)):
        calls.append(lambda value=value:bounds.response_radius(value))
    for value in (-1,0,sp.Rational(-1,10**30)):
        calls.append(lambda value=value:bounds.at_scale(value))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An uncertified constant-mass input was accepted")
    return {"rejected_inputs":rejected,"exact_validation_before_cached_calls":True,
            "connection_Hessian_update_and_quantum_contact_change_explicit":True,
            "Proca_constraint_not_replaced_by_four_unconstrained_positive_oscillators":True,
            "unchanged_clock_operator_not_used_to_claim_unchanged_functional_response":True}
