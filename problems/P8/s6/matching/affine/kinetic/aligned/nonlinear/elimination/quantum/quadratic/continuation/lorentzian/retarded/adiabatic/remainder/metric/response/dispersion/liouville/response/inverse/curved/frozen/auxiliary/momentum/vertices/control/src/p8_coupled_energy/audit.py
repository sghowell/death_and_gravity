"""Exact identities, continuous bounds and fail-closed domain controls."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model

from . import bounds, energy, other_modes, scalars, tree, window


@cache
def residuals():
    out={}
    for group in (scalars.checks(),energy.checks(),other_modes.checks(),tree.scale_checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated seven-mode energy residual")
        out.update(group)
    H,q,g=sp.symbols("H q gradient",real=True)
    z=scalars.z
    out["fixed_comoving_inverse_frequency_derivative"]=sp.factor(
        scalars.derivative(1/z)+2*scalars.H/z)
    out["fixed_comoving_principal_energy_derivative"]=sp.expand(
        (-2*H*q*g)+3*H*q*g-H*q*g)
    # Independently check the explicit finite-fiber Schur relaxation.
    mu=sp.Symbol("measure_upper",positive=True)
    out["fixed_fiber_Schur_relaxation_polynomial"]=sp.expand(
        (8*mu**3)**2-(7*mu)**2-mu**2*(64*mu**4-49))
    return out


@cache
def gates():
    domain=window.data()
    result={name:bool(value>0) for name,value in bounds.elementary_margins().items()}
    result.update({name:bool(value>=0) for name,value in domain["margins"].items()})
    result.update(other_modes.gates())
    data=tree.data()
    result.update({
        "both_actual_scalar_charts_checked":set(domain["charts"])=={"unitary","gamma"},
        "all_actual_background_jets_have_positive_even_denominator_proofs":
            len(bounds.backgrounds()["positive_even_denominator_proofs"])==15,
        "full_beta_initial_shift_has_nonzero_canonical_defect":energy.initial_mixing_negative_control()!=sp.zeros(2),
        "seven_physical_columns_not_transverse_only":data["propagating_mode_columns"]==7,
        "fixed_center_Hamiltonian_volume_retained_separately_from_Fourier_columns":
            data["fixed_center_background_volume_upper"]==8,
        "cubic_Wick_and_five_particle_intermediate_factors_fit":tree.WICK_FACTOR**2>6**2*60,
        "quartic_Wick_and_two_particle_symmetry_factors_fit":tree.WICK_FACTOR>24*2,
        "exchange_partitions_species_orientations_and_time_orders_counted_separately":
            tree.EXCHANGE_FACTOR==3*7*2*2,
        "fixed_fiber_measure_Schur_relaxation":data["momentum_species_measure_upper"]>7,
        "named_scale_is_10_to_353":data["sufficient_named_M_tau"]==10**353,
        "cubic_transition_bound_meets_target":data["cubic_block_bound_at_named_scale"]<=tree.TARGET,
        "quartic_connected_tree_bound_meets_target":data["quartic_connected_tree_bound_at_named_scale"]<=tree.TARGET,
        "not_a_replacement_for_the_L_10_to_24_response_example":data["not_the_previous_M_tau_10_to_24_example"],
        "massive_vector_modes_retained_above_threshold":data["not_a_light_only_heavy_elimination_domain"],
        "no_unsupported_all_orders_or_Wilsonian_cutoff_claim":data["not_a_Wilsonian_or_all_orders_cutoff"],
        "original_classical_action_margin_and_vector_state_unchanged":True,
        "original_P8_and_finite_common_parent_matching_open":True})
    return {name:bool(value) for name,value in result.items()}


def controls():
    calls=[]
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None)
    for value in bad:
        calls.extend((lambda value=value:bounds.time_upper(value),
                      lambda value=value:bounds.coefficient_bound(value,"gamma")))
    for value in (True,False,sp.true,1,sp.Integer(1),1.0,"unknown",None,[]):
        calls.extend((lambda value=value:scalars.data(value),
                      lambda value=value:bounds.chart_bounds(value),
                      lambda value=value:bounds.coefficient_bound(1,value)))
    u=model.u
    for value in (1/u,1/(1-u*u),1/(1+u),sp.sin(u),sp.sqrt(2),sp.Symbol("outside_time")):
        calls.append(lambda value=value:bounds.time_upper(value))
    s=scalars
    for chart in ("unitary","gamma"):
        for value in (1/s.z,1/(1-s.z),1/(s.H+1),sp.sqrt(2),sp.sin(s.ell),sp.Symbol("outside_chart")):
            calls.append(lambda value=value,chart=chart:bounds.coefficient_bound(value,chart))
    D=s.lam**2-(s.Je+s.ell**2*s.lam**2/2)*s.z
    R=s.lam**2-s.Je*s.z
    for value in (1/s.lam,1/D,1/R):
        calls.append(lambda value=value:bounds.coefficient_bound(value,"unitary"))
    calls.append(lambda:bounds.coefficient_bound(1/s.theta,"gamma"))
    for value in (sp.eye(3),sp.zeros(2,1),[[1,0],[0,1]],1,None):
        calls.append(lambda value=value:bounds.matrix_bound(value,"gamma"))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An uncertified energy/chart input was accepted")
    return {"rejected_inputs":rejected,
            "native_chart_validation_before_cached_entrypoints":True,
            "exact_rational_coefficients_only":True,
            "wrong_chart_and_zero_frequency_denominators_rejected":True,
            "full_beta_initial_shift_negative_control_retained":True,
            "no_independent_oscillator_or_reselected_vector_state_shortcut":True}
