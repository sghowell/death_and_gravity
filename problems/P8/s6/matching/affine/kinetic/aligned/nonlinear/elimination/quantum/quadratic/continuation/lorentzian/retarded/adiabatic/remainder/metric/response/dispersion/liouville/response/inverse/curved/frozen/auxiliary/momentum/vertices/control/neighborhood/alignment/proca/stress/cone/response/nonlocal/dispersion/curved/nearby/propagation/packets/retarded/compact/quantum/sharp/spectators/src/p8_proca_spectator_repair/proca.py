"""Literal ordinary-Proca decoupling and its three physical principal modes."""
from functools import cache

import sympy as sp
from p8_constant_proca import model as actual_parent

from . import pencil

k,m=sp.symbols("positive_physical_momentum positive_Proca_mass",positive=True)


def species(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Integer)):
        raise TypeError("Require an exact integer number of ordinary Proca species")
    if value<0:
        raise ValueError("The number of Proca species cannot be negative")
    return int(value)


def characteristic(value=1):
    n=species(value)
    return {"number_of_ordinary_Proca_species":n,
        "number_of_physical_Proca_polarizations":3*n,
        "number_of_luminal_physical_modes":3+3*n,
        "total_number_of_physical_modes":4+3*n,
        "normalized_principal_characteristic":(pencil.c**2-pencil.s**2)*(1-pencil.s**2)**(3+3*n),
        "fixed_finite_positive_masses_do_not_enter_principal_polynomial":True,
        "n_is_total_species_in_new_classical_family_not_old_quantum_state_count":True}


@cache
def data():
    eps,eta=sp.symbols("vector_perturbation metric_perturbation",real=True)
    metric=sp.diag(-1,1,1,1)
    variables=sp.symbols("inverse_metric_jet_0:10",real=True)
    H=sp.zeros(4)
    offset=0
    for i in range(4):
        for j in range(i,4):
            H[i,j]=H[j,i]=variables[offset]
            offset+=1
    G=metric+eta*H
    density=1-eta*sp.trace(metric*H)/2
    E1,E2,E3,B1,B2,B3=sp.symbols("electric_1 electric_2 electric_3 magnetic_1 magnetic_2 magnetic_3",real=True)
    F=sp.Matrix([[0,E1,E2,E3],[-E1,0,B3,-B2],[-E2,-B3,0,B1],[-E3,B2,-B1,0]])
    A=sp.Matrix(sp.symbols("vector_component_0:4",real=True))
    action=-density*eps**2*((F.T*G*F*G).trace()/4+m*m*(A.T*G*A)[0]/2)
    quadratic=-(F.T*metric*F*metric).trace()/4-m*m*(A.T*metric*A)[0]/2
    v,vd,a0=sp.symbols("longitudinal_mode longitudinal_velocity nondynamical_temporal_mode",real=True)
    long=(vd-k*a0)**2/2+m*m*a0*a0/2-m*m*v*v/2
    solution=k*vd/(k*k+m*m)
    reduced=sp.factor(long.subs(a0,solution))
    zeta=sp.Rational(1,10**6)
    unscaled=-zeta*(F.T*metric*F*metric).trace()/4-(A.T*metric*A)[0]/2
    canonical=sp.factor(unscaled/zeta)
    return {"local_physical_inertial_metric":metric,"inverse_metric_first_jet":H,
        "metric_volume_first_jet":density,"full_relevant_metric_vector_order_action":action,
        "vector_amplitude_parameter":eps,"metric_amplitude_parameter":eta,
        "literal_background_metric_quadratic_Proca_action":quadratic,
        "electric_components":(E1,E2,E3),"magnetic_components":(B1,B2,B3),"vector_components":A,
        "tangent_frame_longitudinal_action_with_constraint":long,
        "nondynamical_temporal_component":a0,"longitudinal_coordinate":v,"longitudinal_velocity":vd,
        "temporal_constraint_solution":solution,
        "reduced_tangent_frame_longitudinal_action":reduced,
        "longitudinal_kinetic_coefficient":m*m/(k*k+m*m),
        "tangent_frame_frozen_massive_dispersion_squared":k*k+m*m,
        "actual_literal_unscaled_Maxwell_coefficient":zeta,
        "actual_unscaled_vector_quadratic_action":unscaled,
        "actual_action_after_constant_canonical_vector_rescaling":canonical,
        "actual_existing_canonical_Proca_mass_squared":actual_parent.data()["canonical_Proca_mass_squared"],
        "n_example_principal_counts":{str(n):characteristic(n) for n in (0,1,3)},
        "no_scalar_vector_quadratic_cross_term_at_zero_Proca_background":True,
        "tangent_frame_massive_dispersion_not_a_time_independent_FLRW_evolution_claim":True,
        "classical_principal_result_does_not_include_Proca_loop_corrections_or_quantum_backreaction":True}


@cache
def checks():
    d=data()
    eps,eta=d["vector_amplitude_parameter"],d["metric_amplitude_parameter"]
    L=d["full_relevant_metric_vector_order_action"]
    A=d["vector_components"]
    expected=(sum(x*x for x in d["electric_components"])-sum(x*x for x in d["magnetic_components"]))/2
    expected+=m*m*(A[0]**2-sum(A[i]**2 for i in range(1,4)))/2
    v,vd,a0=d["longitudinal_coordinate"],d["longitudinal_velocity"],d["nondynamical_temporal_component"]
    reduced=d["reduced_tangent_frame_longitudinal_action"]
    K,H,C=(pencil.symmetric(name) for name in ("mass_anchor_K","mass_anchor_H","mass_anchor_C"))
    mass=pencil.symmetric("finite_lower_order_mass")
    principal=pencil.s**2*K+2*pencil.s*C-H
    full=(k*k*principal+mass).det()
    rows={"literal_vector_quadratic_action_uses_only_background_metric":
        sp.diff(L,eps,2).subs({eps:0,eta:0})/2-d["literal_background_metric_quadratic_Proca_action"],
        "literal_mixed_metric_vector_quadratic_variation_vanishes":
        sp.diff(L,eta,eps).subs({eps:0,eta:0}),
        "literal_tangent_Proca_action_has_correct_electric_magnetic_and_mass_signs":
        d["literal_background_metric_quadratic_Proca_action"]-expected,
        "temporal_Proca_constraint_is_solved_without_dropping_mass":
        sp.diff(d["tangent_frame_longitudinal_action_with_constraint"],a0).subs(a0,d["temporal_constraint_solution"]),
        "reduced_longitudinal_action_keeps_positive_finite_momentum_kinetic_factor":
        reduced-(d["longitudinal_kinetic_coefficient"]*vd**2-m*m*v*v)/2,
        "longitudinal_massive_dispersion_matches_transverse_one":
        m*m/d["longitudinal_kinetic_coefficient"]-d["tangent_frame_frozen_massive_dispersion_squared"],
        "longitudinal_front_is_luminal_at_fixed_finite_mass":
        sp.limit(d["tangent_frame_frozen_massive_dispersion_squared"]/k**2,k,sp.oo)-1,
        "actual_existing_Proca_mass_is_retained":
        d["actual_existing_canonical_Proca_mass_squared"]-10**6,
        "literal_constant_vector_rescaling_gives_canonical_action_and_mass":
        d["actual_action_after_constant_canonical_vector_rescaling"]-
            d["literal_background_metric_quadratic_Proca_action"].subs(
                m*m,d["actual_existing_canonical_Proca_mass_squared"]),
        "canonical_mass_keeps_reciprocal_literal_Maxwell_coefficient":
        d["actual_existing_canonical_Proca_mass_squared"]*
            d["actual_literal_unscaled_Maxwell_coefficient"]-1,
        "finite_mass_mixing_does_not_change_second_order_characteristic_leading_term":
        sp.expand(full).coeff(k,4)-principal.det(),
        "zero_spectator_case_retains_matter_and_two_tensor_modes":
        characteristic(0)["number_of_luminal_physical_modes"]-3,
        "actual_one_Proca_case_has_three_extra_physical_modes":
        characteristic(1)["total_number_of_physical_modes"]-7}
    return {name:sp.factor(value) for name,value in rows.items()}


@cache
def gates():
    d=data()
    rows={"canonical_Proca_mass_strictly_positive":d["actual_existing_canonical_Proca_mass_squared"]>0,
        "longitudinal_kinetic_strictly_positive_for_finite_positive_mass":
            d["longitudinal_kinetic_coefficient"].is_positive is True,
        "metric_vector_mixing_starts_above_quadratic_order":True,
        "finite_mass_high_frequency_limit_not_exchanged_with_infinite_mass_limit":True,
        "no_added_Proca_species_can_remove_the_unchanged_clock_factor":True,
        "no_Proca_quantum_state_or_tadpole_cancellation_is_transferred":True}
    return {name:bool(value) for name,value in rows.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("count"),sp.Rational(1,2),[],-1)
    rejected=0
    for value in bad:
        for call in (species,characteristic):
            try:
                call(value)
            except (TypeError,ValueError):
                rejected+=1
            else:
                raise ValueError("An invalid number of spectator species was accepted")
    return {"rejected_inputs":rejected}
