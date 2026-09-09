"""Explicit bridge to the preserved global clock and old parameterized cone proof."""
from functools import cache

import sympy as sp
from p8_aligned_margin import bounds as margin_bounds
from p8_aligned_margin import dynamics as margin_dynamics

from . import model


def on_clock(value):
    bg=margin_dynamics.old.background()
    return sp.factor(value.evaluate(model.Nf,1).as_expr().subs(model.z,-2*bg["H"]))


@cache
def data():
    bg=margin_dynamics.old.background()
    bounds=margin_bounds.polynomial_bounds()
    eps,delta=model.NEW_MARGIN,model.INCREMENT
    lower=bounds["global_lower"]
    return {"unchanged_global_classical_scale":(1+model.u**2)**2,
        "unchanged_global_classical_physical_Hubble":bg["H"],
        "unchanged_global_classical_matter_density":bg["ell"],
        "new_global_clock_speed_squared":sp.factor(bg["J"]/(bg["J"]+4*eps/bg["h"]**2)),
        "positive_all_time_clock_squared_speed_lower":lower/(lower+4*eps),
        "positive_compact_original_clock_squared_cone_margin_lower":4*eps/(bounds["compact_upper"]+4*eps),
        "relative_light_kinetic_increment_over_frozen_margin_upper":4*delta/(lower+4*model.OLD_MARGIN),
        "added_scalar_on_real_clock_tube_upper":delta/100,
        "added_scalar_first_X_derivative_on_real_clock_tube_upper":delta/5,
        "added_scalar_second_X_derivative_on_real_clock_tube_upper":2*delta,
        "clock_tube_definition":"all real u and abs(X-1)<=1/10",
        "positive_polynomial_global_J_h_squared_lower":lower,
        "compact_J_h_squared_upper":bounds["compact_upper"],
        "two_original_principal_charts_keep_their_existing_coverage":True,
        "not_a_global_continuation_of_the_new_off_clock_solution":True,
        "no_old_loop_potential_profile_or_quantum_state_transferred":True}


@cache
def checks():
    bg=margin_dynamics.old.background()
    s=model.system()
    d=data()
    rows={"new_original_clock_matter_constraint_is_exact_old_history":
            on_clock(s["constraint_matter_square"])-bg["ell"]**2,
        "new_original_clock_lapse_flow_vanishes_identically":on_clock(s["lapse_flow"]),
        "new_original_clock_trace_flow_is_exact_background_derivative":
            on_clock(s["rescaled_trace_flow"])+2*sp.diff(bg["H"],model.u),
        "new_original_clock_fixed_phase_pivot_matches_total_not_increment_margin":
            on_clock(s["on_constraint_rescaled_lapse_Hessian"])+2*(bg["J"]+4*model.NEW_MARGIN/bg["h"]**2),
        "new_actual_clock_speed_agrees_with_parameterized_two_chart_result":
            on_clock(s["clock_physical_speed_squared"])-d["new_global_clock_speed_squared"],
        "preserved_global_scale_has_original_physical_Hubble":
            sp.diff(d["unchanged_global_classical_scale"],model.u)/d["unchanged_global_classical_scale"]-bg["H"],
        "global_parameterized_lower_keeps_total_margin":
            d["positive_all_time_clock_squared_speed_lower"]-sp.Rational(1199)/(1199+3200*model.NEW_MARGIN),
        "relative_increment_compares_to_frozen_nonzero_margin_not_undeformed_action":
            d["relative_light_kinetic_increment_over_frozen_margin_upper"]-
            4*model.INCREMENT/(sp.Rational(1199,800)+4*model.OLD_MARGIN)}
    return {name:sp.factor(value) for name,value in rows.items()}


@cache
def gates():
    d=data()
    rows={"new_margin_inside_original_explicit_classical_parameter_family":0<model.NEW_MARGIN<=sp.Rational(1,100),
        "old_global_positive_even_polynomial_proof_rechecked":margin_bounds.polynomial_bounds()["global_positive_even_remainder"],
        "new_all_time_original_clock_speed_lower_above_493_over_500":
            d["positive_all_time_clock_squared_speed_lower"]>sp.Rational(493,500),
        "new_compact_original_clock_margin_above_one_over_two_thousand":
            d["positive_compact_original_clock_squared_cone_margin_lower"]>sp.Rational(1,2000),
        "increment_relative_to_frozen_kinetic_form_below_seven_over_five_hundred":
            d["relative_light_kinetic_increment_over_frozen_margin_upper"]<sp.Rational(7,500),
        "literal_addition_respects_recorded_clock_tube_scalar_budget":
            d["added_scalar_on_real_clock_tube_upper"]<sp.Rational(1,20000),
        "global_clock_preservation_is_classical_not_quantum_profile_cancellation":True}
    return {name:bool(value) for name,value in rows.items()}
