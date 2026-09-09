"""Exact complete-frequency scalar source bridge on the new constraint solution."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import principal
from p8_proca_nearby_packets import finite
from p8_proca_nearby_retarded import source
from p8_proca_retuned_margin import model


@cache
def data():
    s=model.system()
    e,R,w=sp.symbols("new_conformal_factor new_hat_scale new_rescaled_matter_density",positive=True)
    p=principal
    expressions={p.a:-3*model.N*e/4,p.m:model.N/(2*e**3),p.g:model.N*e/2,
        p.r:-model.N/(2*e**3),p.rn:-e*s["M"].as_expr(),
        p.h:s["on_constraint_rescaled_lapse_Hessian"].as_expr()/e,
        p.ell:w/e,p.alpha:((2*s["A"]*model.zf+s["B"])/3).as_expr(),
        p.beta:2*w*s["M"].as_expr(),p.H:s["hat_Hubble"].as_expr()}
    d=source.data()
    return {"literal_new_coefficient_substitution":expressions,
        "new_constraint_matter_square":s["constraint_matter_square"].as_expr(),
        "new_conformal_fourth_power":(model.Nf**2*s["h"]/s["D"]).as_expr(),
        "full_original_weighted_phase_generator":d["weighted_first_order_generator"],
        "full_original_polynomial_q_coefficient_matrices":d["polynomial_q_coefficients"],
        "generic_exact_old_chart_to_packet_map":d["old_to_packet_canonical_map"],
        "generic_complete_packet_Laurent_generator":finite.k*finite.data()["J"]+finite.data()["L0"]
            +finite.data()["L1"]/finite.k+finite.data()["L2"]/finite.k**2+finite.data()["L3"]/finite.k**3,
        "new_source_action_per_coordinate_four_volume":R**3*model.N*e**3*sp.Symbol("J")*p.chi,
        "new_original_chart_source_vector":sp.Matrix([0,0,0,model.N*e**3]),
        "observable_is_new_relational_matter_linearization":p.chi,
        "new_coordinate_matter_frequency":model.N/(e*R),
        "complex_q_definition":"q=(zeta_1^2+zeta_2^2+zeta_3^2)/R(u)^2",
        "all_symbolic_source_and_Laurent_identities_are_universal_algebra_only":True,
        "new_actual_coefficients_and_state_not_inherited_from_fast_clock_parent":True}


@cache
def checks():
    rows={"new_source_bridge_"+name:value for name,value in source.checks().items()}
    rows.update({"new_exact_Laurent_"+name:value for name,value in finite.checks().items()})
    s=model.system()
    B=s["rescaled_principal_B"]
    F,G=s["lapse_flow"],s["rescaled_trace_flow"]
    moving=B.diff(model.uf)+B.diff(model.Nf)*F+B.diff(model.zf)*G-s["conformal_log_flow"]*B
    Gbar=moving-s["hat_Hubble"]*B-s["m0"]*B*B-model.Nf*s["constraint_matter_square"]
    rows["new_complete_moving_gradient_matches_new_rational_principal_formula"]=(
        Gbar-s["gradient_fixed_lapse_part"]-s["gradient_lapse_velocity_part"]*F).as_expr()
    rows["new_actual_clock_frequency_uses_full_moving_gradient"]=(
        -4*s["h"]*s["M"]**2*Gbar/(s["D"]*s["on_constraint_rescaled_lapse_Hessian"])
        -s["clock_physical_speed_squared"]).as_expr()
    e,N,M,P,kc,c,R=sp.symbols("e N M P K_clock c R",nonzero=True)
    rows["new_clock_kinetic_retains_fixed_phase_pivot_and_conformal_cube"]=sp.factor(
        -(P/e)/(4*(e*M)**2)+P/(4*e**3*M**2))
    rows["new_clock_gradient_over_scale_equals_kinetic_times_coordinate_frequency_squared"]=sp.factor(
        (N*N/e**2)*kc*c*c/R**2-kc*(N*c/(e*R))**2)
    rows["new_matter_gradient_over_scale_equals_kinetic_times_matter_frequency_squared"]=sp.factor(
        N*e/R**2-(e**3/N)*(N/(e*R))**2)
    return rows
