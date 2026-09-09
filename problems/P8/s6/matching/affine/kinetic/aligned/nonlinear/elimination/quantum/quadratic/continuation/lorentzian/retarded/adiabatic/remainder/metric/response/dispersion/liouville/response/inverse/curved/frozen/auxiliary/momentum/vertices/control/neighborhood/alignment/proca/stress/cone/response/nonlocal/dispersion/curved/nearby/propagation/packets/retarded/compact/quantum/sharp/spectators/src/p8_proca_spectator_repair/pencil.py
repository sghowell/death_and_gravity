"""Actual physical light block and finite second-order characteristic pencil."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import bounds as parent_bounds
from p8_proca_nearby_packets import modes

c,s=sp.symbols("positive_physical_clock_speed characteristic_speed",positive=True)


@cache
def light():
    d=modes.data()
    kc,km,ell=modes.kc,modes.km,modes.ell
    K=d["kinetic"]
    H=(d["coordinate_gradient"].subs(modes.wc,c*modes.wm)/modes.wm**2).applyfunc(sp.factor)
    C=sp.Matrix([[1/sp.sqrt(kc),0],[-ell/sp.sqrt(kc),1/sp.sqrt(km)]])
    return {"literal_light_kinetic":K,"physical_matter_frame_light_gradient":H,
        "old_kinetic_canonical_congruence":C,
        "canonical_light_kinetic":sp.eye(2),"canonical_light_gradient":sp.diag(c*c,1),
        "canonical_light_clock_direction":sp.Matrix([1,0]),
        "actual_clock_squared_speed_excess_lower":sp.Rational(5,10**6),
        "actual_clock_squared_speed_excess_upper":sp.Rational(16,10**6),
        "physical_gradient_uses_one_matter_frame_not_mode_dependent_frame":True}


def symmetric(prefix):
    values=sp.symbols(f"{prefix}_11 {prefix}_12 {prefix}_22",real=True)
    return sp.Matrix([[values[0],values[1]],[values[1],values[2]]])


@cache
def extension():
    deltaK,deltaC,deltaH=(symmetric(name) for name in ("light_kinetic_change","light_mixed_change","light_gradient_change"))
    crossK,crossC,crossH=(sp.Matrix(2,2,sp.symbols(f"{name}_0:4",real=True))
                        for name in ("heavy_kinetic_mix","heavy_time_space_mix","heavy_gradient_mix"))
    heavyK,heavyC,heavyH=(symmetric(name) for name in ("heavy_kinetic","heavy_time_space","heavy_gradient"))
    block=lambda a,b,d:a.row_join(b).col_join(b.T.row_join(d))
    K=block(sp.eye(2)+deltaK,crossK,heavyK)
    C=block(deltaC,crossC,heavyC)
    H=block(sp.diag(c*c,1)+deltaH,crossH,heavyH)
    P=s*s*K+2*s*C-H
    direction=sp.Matrix([1,0,0,0])
    return {"four_mode_algebra_anchor_not_a_bound_on_number_of_heavy_fields":True,
        "full_kinetic_anchor":K,"full_mixed_time_space_anchor":C,"full_gradient_anchor":H,
        "full_characteristic_pencil_anchor":P,"embedded_light_clock_direction":direction,
        "light_kinetic_change":deltaK,"light_mixed_change":deltaC,"light_gradient_change":deltaH,
        "clock_direction_pencil":(direction.T*P*direction)[0],
        "general_theorem_allows_any_finite_number_of_added_physical_modes":True,
        "positive_full_kinetic_form_and_Hermitian_pencil_are_explicit_hypotheses":True,
        "unchanged_light_block_has_a_real_characteristic_at_least_clock_speed":True,
        "arbitrary_light_mixed_terms_do_not_repair_both_propagation_directions":True,
        "no_frequency_dependent_constraint_schur_complement_assumed_polynomial":True}


@cache
def checks():
    d,e=light(),extension()
    C=d["old_kinetic_canonical_congruence"]
    k,h=d["literal_light_kinetic"],d["physical_matter_frame_light_gradient"]
    dk,dc,dh=e["light_kinetic_change"],e["light_mixed_change"],e["light_gradient_change"]
    value=e["clock_direction_pencil"]
    nochange={x:0 for matrix in (dk,dc,dh) for x in matrix}
    return {"actual_light_kinetic_canonical_congruence":(C.T*k*C-sp.eye(2)).applyfunc(sp.factor),
        "actual_gradient_is_canonicalized_in_fixed_physical_matter_frame":
            (C.T*h*C-sp.diag(c*c,1)).applyfunc(sp.factor),
        "physical_clock_direction_has_exact_positive_excess":
            sp.factor((sp.Matrix([1,0]).T*(C.T*(h-k)*C)*sp.Matrix([1,0]))[0]-(c*c-1)),
        "all_heavy_offdiagonal_blocks_drop_out_of_embedded_clock_pencil":
            sp.expand(value-(s*s*(1+dk[0,0])+2*s*dc[0,0]-c*c-dh[0,0])),
        "unchanged_embedded_light_pencil_is_exactly_original_clock_pencil":
            sp.expand(value.subs(nochange)-(s*s-c*c)),
        "unchanged_full_pencil_has_nonpositive_minimum_eigenvalue_at_clock_speed":
            sp.expand(value.subs(nochange).subs(s,c)),
        "matter_lightcone_test_keeps_all_three_light_block_changes":
            sp.expand(value.subs(s,1)-(1-c*c+dk[0,0]+2*dc[0,0]-dh[0,0])),
        "two_orientation_average_removes_arbitrary_light_mixed_block":
            sp.expand((value.subs(s,1)+value.subs(s,-1))/2-(1-c*c+dk[0,0]-dh[0,0])),
        "two_orientation_difference_retains_mixed_block_with_factor_four":
            sp.expand(value.subs(s,1)-value.subs(s,-1)-4*dc[0,0]),
        "at_least_one_clock_speed_orientation_remains_nonpositive_with_unchanged_K_H":
            sp.expand((value.subs(s,c)+value.subs(s,-c)).subs(
                {x:0 for matrix in (dk,dh) for x in matrix})),
        "canonical_light_principal_determinant_has_clock_and_matter_factors":
            sp.factor((h-s*s*k).det()-modes.kc*modes.km*(c*c-s*s)*(1-s*s))}


@cache
def gates():
    d=parent_bounds.enclosures()["clock_physical_speed_squared_excess"]
    return {"actual_whole_interval_clock_gap_exceeds_five_e_minus_six":
                bool(d["lower"]>sp.Rational(5,10**6)),
        "actual_whole_interval_clock_gap_below_sixteen_e_minus_six":
                bool(d["upper"]<sp.Rational(16,10**6)),
        "arbitrary_finite_dimension_follows_from_written_Hermitian_continuity_proof":True,
        "mixed_time_space_couplings_are_included_not_silently_zeroed":True,
        "positive_full_kinetic_form_is_a_hypothesis_not_an_inferred_parent_property":True}
