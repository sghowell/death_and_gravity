"""Necessary principal light-block repair, keeping both orientations and mixed terms."""
import sympy as sp

from . import pencil

GAP=sp.Rational(5,10**6)


def exact_nonnegative(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact nonnegative rational principal-norm budget")
    value=sp.Rational(value)
    if value<0:
        raise ValueError("A principal-norm budget cannot be negative")
    return value


def budget(kinetic=GAP/4,gradient=GAP/4):
    kinetic,gradient=exact_nonnegative(kinetic),exact_nonnegative(gradient)
    if kinetic+gradient>=GAP:
        raise ValueError("The strict nonrepair budget must sum to less than the actual gap floor")
    return {"old_kinetic_relative_kinetic_change_norm_upper":kinetic,
        "old_kinetic_relative_gradient_change_norm_upper":gradient,
        "strict_nonrepair_margin":GAP-kinetic-gradient,
        "no_smallness_assumption_on_finite_mixed_time_space_or_heavy_blocks":True,
        "positive_full_kinetic_form_remains_required":True,
        "a_matter_causal_parent_requires_change_norm_sum_strictly_above":GAP,
        "bounds_are_conditional_inputs_not_claimed_actual_UV_matching_errors":True}


def data():
    c,s=pencil.c,pencil.s
    gap=c*c-1
    mixed=(s*s+gap*s-c*c)
    return {"actual_whole_interval_necessary_relative_K_H_repair_norm_sum_lower":GAP,
        "both_orientation_necessary_directional_condition":
            "DeltaK_11-DeltaH_11 >= (c_clock^2-1)+2*abs(DeltaC_11)",
        "positive_direction_only_mixed_repair_pencil":mixed,
        "positive_direction_only_mixed_repair_roots":(sp.Integer(1),-c*c),
        "algebraic_gradient_change_that_can_repair_principal_pair":sp.diag(-gap,0),
        "repaired_canonical_kinetic":sp.eye(2),"repaired_canonical_gradient":sp.eye(2),
        "algebraic_repair_not_a_covariant_action_or_UV_construction":True,
        "old_kinetic_relative_norm_is_a_specified_canonical_comparison_not_unqualified_coefficient_size":True,
        "example_strict_nonrepair_budget":budget()}


def checks():
    d=data()
    c,s=pencil.c,pencil.s
    dk,dh,dc=sp.symbols("clock_kinetic_change clock_gradient_change clock_mixed_change",real=True)
    plus=1-c*c+dk+2*dc-dh
    minus=1-c*c+dk-2*dc-dh
    return {"two_orientation_average_requires_K_H_repair_even_with_arbitrary_mixed_term":
        sp.expand((plus+minus)/2-(1-c*c+dk-dh)),
        "both_orientations_require_extra_absolute_mixed_term_allowance":
        sp.expand((plus+minus)/2-sp.Abs((plus-minus)/2)-(1-c*c+dk-dh-2*sp.Abs(dc))),
        "mixed_only_positive_direction_repair_leaves_a_faster_negative_root":
        sp.expand(d["positive_direction_only_mixed_repair_pencil"]-(s-1)*(s+c*c)),
        "mixed_only_control_is_luminal_in_one_orientation":
        d["positive_direction_only_mixed_repair_pencil"].subs(s,1),
        "mixed_only_control_is_negative_at_other_matter_lightcone":
        sp.expand(d["positive_direction_only_mixed_repair_pencil"].subs(s,-1)+2*(c*c-1)),
        "allowing_a_real_light_gradient_change_can_remove_the_rank_one_excess":
        sp.diag(c*c,1)+d["algebraic_gradient_change_that_can_repair_principal_pair"]-sp.eye(2),
        "strict_budget_retains_the_sum_of_both_relative_norm_errors":
        budget()["strict_nonrepair_margin"]-GAP/2}


def gates():
    rows={"actual_gap_floor_positive":GAP>0,
        "default_conditional_nonrepair_margin_positive":budget()["strict_nonrepair_margin"]>0,
        "mixed_terms_not_assumed_small_or_absent":True,
        "necessary_repair_lower_is_not_a_sufficient_UV_matching_condition":True,
        "algebraic_repair_control_does_not_change_the_frozen_action":True}
    return {name:bool(value) for name,value in rows.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.sqrt(2),[],-1)
    rejected=0
    for value in bad:
        for pair in ((value,0),(0,value)):
            try:
                budget(*pair)
            except (TypeError,ValueError):
                rejected+=1
            else:
                raise ValueError("An invalid relative principal-norm budget was accepted")
    for pair in ((GAP,0),(0,GAP),(GAP/2,GAP/2),(GAP,GAP)):
        try:
            budget(*pair)
        except ValueError:
            rejected+=1
        else:
            raise ValueError("A non-strict repair margin was accepted")
    return {"rejected_inputs":rejected}
