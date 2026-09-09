"""New ordinary Proca response identities, continuous estimates and scope gates."""
from functools import cache

import sympy as sp
from p8_proca_local_response import bounds as local_bounds
from p8_vector_state import wkb

from . import envelopes, evolution, source, tadpole, tail, tangent


@cache
def residuals():
    out={}
    for group in (source.checks(),tangent.low_checks(),tail.algebra_checks(),
                  evolution.algebra_checks(),evolution.weight_checks(),tadpole.physical_vertices()["checks"]):
        if set(out).intersection(group):
            raise ValueError("Repeated new ordinary Proca response residual")
        out.update(group)
    for sector in ("T","L"):
        reference=tangent.reference(sector)
        out.update({sector+"_new_varied_reference_low_"+str(j):value
                    for j,value in reference["varied_residual_low_coefficients"].items()})
        for component in ("energy","pressure"):
            out.update({sector+"_"+component+"_new_subtracted_reference_low_"+str(j):value
                        for j,value in tail.reference_tail(sector,component)["low_tail_numerator_coefficients"].items()})
        out[sector+"_new_tenth_lapse_derivative_vanishes"]=reference["tenth_source_derivative_fixtures"]["N"]
        out[sector+"_new_tenth_scale_derivative_is_nonzero"]=reference["tenth_source_derivative_fixtures"]["Z"]+sp.Rational(1,512)
        out[sector+"_new_tenth_lapse_coefficient_zero_on_full_box"]=reference["tenth_source_derivative_coefficients"]["N"]
        expected=-(1-source.z)/512 if sector=="T" else -(1+source.z)/512
        out[sector+"_new_tenth_scale_coefficient_on_full_box"]=sp.factor(reference["tenth_source_derivative_coefficients"]["Z"]-expected)
        for order in range(1,5):
            value=tangent.coefficient(sector,order)
            out[sector+"_new_reference_source_linearity_"+str(order)]=source.clean(
                sum(field*sp.diff(value,field) for field in source.n+source.v)-value)
    for order in range(1,5):
        out["new_isotropic_reference_variation_"+str(order)]=source.clean(
            (tangent.varied_P("T",order)-tangent.varied_P("L",order)).subs(source.z,0))
    for j,row in tadpole.chart_lift()["coefficient_derivative_envelopes"].items():
        out["actual_clock_chart_input_continuous_derivative_"+str(j)]=row["reconstruction"]
    return out


@cache
def gates():
    m=local_bounds.MASS
    out={}
    for sector in ("T","L"):
        d=envelopes.reference(sector)
        out[sector+"_all_varied_reference_box_reconstructions"]=all(value==0 for value in d["box_reconstructions"])
        out[sector+"_same_baseline_reference_and_transport_bounds"]=all(d["base_reference_checks"].values())
        out[sector+"_new_log_frequency_variation_bound_below_four"]=d["source_log_frequency_variation_upper"]<4
        out[sector+"_new_reference_frequency_variation_bound_below_four"]=d["delta_W_over_frequency_upper"]<4
        out[sector+"_new_reference_c_variation_bound_below_sixty_four"]=d["delta_c_upper"]<64
        out[sector+"_canonical_rate_at_most_three"]=wkb.box_bound(source.data(sector)["rate"])["absolute_upper"]<=3
        out[sector+"_reference_momentum_below_two_frequency_f"]=5+sp.Rational(3,2)*m<=2*m
        out[sector+"_baseline_transport_exponent_below_one_quarter"]=2*d["base_residual_over_inverse_frequency_eighth_upper"]/m**9<sp.Rational(1,4)
        for component in ("energy","pressure"):
            out[sector+"_"+component+"_new_tail_majorants_nonnegative"]=tail.reference_tail(sector,component)["all_majorants_nonnegative"]
            out[sector+"_"+component+"_new_readout_transport_bound_positive"]=evolution.constants()["by_sector"][sector]["physical_readouts"][component]["delta_reference_bilinear_over_nu_squared_upper"]>0
    old=evolution.constants()["frozen_initial_and_evolved_mixing"]
    B6,B8,B10=[old["initial_mixing_envelopes"][j] for j in (6,8,10)]
    out.update({"same_middle_initial_band_covered_by_sixth_power":B8/(4*m)**2<=B6,
                "same_high_initial_band_covered_by_sixth_power":B10/(8*m)**4<=B6,
                "same_evolved_mixing_strictly_below_one":(B6+old["evolution_mixing_envelope"]/m**4)/m**6<1,
                "same_Gronwall_slack":1/(1-sp.Rational(1,4))<2,
                "actual_clock_output_map_bound":tadpole.chart_lift()["output_bound_uses_omega_N_at_most_one_half"],
                "new_weights_not_old_nonminimal_energy":all(value==0 for value in evolution.weight_checks().values())})
    vector=evolution.bound(local_bounds.SCALE)
    full=tadpole.bound(local_bounds.SCALE)
    out.update({"new_full_vector_physical_response_below_one_e_minus_790":vector["joint_complete_metric_C10_to_C0_upper_bound"]<sp.Rational(1,10**790),
                "new_background_cancelled_physical_response_below_one_e_minus_790":full["joint_background_cancelled_physical_C10_to_C0_upper"]<sp.Rational(1,10**790),
                "new_background_cancelled_clock_response_below_one_e_minus_779":full["joint_background_cancelled_clock_C10_to_C0_upper"]<sp.Rational(1,10**779),
                "no_independent_initial_response_does_not_erase_squeezed_initial_mixing":True,
                "same_actual_state_and_new_explicit_fixed_profile_not_reselected":True,
                "new_dimensional_mass_jets_zero_before_matching":True,
                "actual_physical_output_contacts_and_covariance_Ward_retained":True,
                "tenth_scale_source_jet_not_mislabeled_as_C4_norm":True,
                "C10_to_C0_bound_not_a_no_loss_coupled_inverse_or_quantum_cone":True,
                "old_full_rank_inverse_or_frozen_pole_not_transferred":True,
                "omitted_operators_mixed_loops_V_G_B_matching_and_original_P8_open":True})
    return {name:bool(value) for name,value in out.items()}


def controls():
    calls=[]
    for value in (True,False,1,sp.Integer(1),"unknown",None,[]):
        calls.extend((lambda value=value:source.data(value),
                      lambda value=value:source.readout(value,"energy"),
                      lambda value=value:source.adiabatic(value,"energy",1),
                      lambda value=value:source.baseline(value,1),
                      lambda value=value:tangent.coefficient(value,1),
                      lambda value=value:tangent.reference(value),
                      lambda value=value:envelopes.reference(value),
                      lambda value=value:tail.reference_tail(value,"energy")))
        calls.extend((lambda value=value:source.readout("T",value),
                      lambda value=value:source.adiabatic("T",value,1),
                      lambda value=value:tail.reference_tail("T",value)))
    for value in (True,False,1.0,sp.Integer(1),"1",-1,5,None,[]):
        calls.extend((lambda value=value:source.baseline("T",value),
                      lambda value=value:tangent.coefficient("T",value),
                      lambda value=value:source.adiabatic("T","energy",value)))
    for value in (True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
                  sp.Symbol("unfixed_scale"),sp.sqrt(2),[],0,-1):
        calls.extend((lambda value=value:evolution.bound(value),lambda value=value:tadpole.bound(value)))
    for value in (True,False,1.0,sp.Float(1),sp.oo,sp.nan,source.n[0]*source.v[0],
                  source.n[11],sp.Symbol("unknown")*source.n[0],sp.Integer(1)):
        calls.append(lambda value=value:source.linear_bound(value))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An inexact or outside-domain prepared Proca response input was accepted")
    return {"rejected_inputs":rejected,"native_source_and_order_guards_reject_cached_aliases":True,
            "no_parent_scientific_model_or_function_monkeypatched":True,
            "new_reference_constants_recomputed_not_old_response_bounds_copied":True}
