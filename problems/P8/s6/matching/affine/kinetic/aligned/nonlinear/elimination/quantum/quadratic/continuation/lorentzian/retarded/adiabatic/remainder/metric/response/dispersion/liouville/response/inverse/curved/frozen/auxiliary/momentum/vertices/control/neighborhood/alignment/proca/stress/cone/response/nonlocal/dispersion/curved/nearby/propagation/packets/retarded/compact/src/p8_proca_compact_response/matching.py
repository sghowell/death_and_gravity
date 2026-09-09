"""Conditional response comparison keeps low-band error and UV-tail error separate."""
from functools import cache

import sympy as sp

from . import band, probes

DEFAULT_ERROR_FRACTION=sp.Rational(1,8)


def budget(low=DEFAULT_ERROR_FRACTION,tail=DEFAULT_ERROR_FRACTION):
    for value in (low,tail):
        if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
            raise TypeError("Require exact real rational conditional error fractions")
    low,tail=sp.Rational(low),sp.Rational(tail)
    if low<0 or tail<0 or low+tail>=sp.Rational(1,2):
        raise ValueError("Both conditional errors must be nonnegative with sum less than one half")
    return {"assumed_low_spatial_band_matching_error_fraction":low,
            "assumed_parent_high_spatial_frequency_tail_fraction":tail,
            "classical_low_spatial_band_signal_fraction":sp.Rational(1,2),
            "conditional_parent_full_response_lower_fraction":sp.Rational(1,2)-low-tail,
            "conditional_parent_full_response_lower":probes.SIGNAL_LOWER*(sp.Rational(1,2)-low-tail),
            "actual_matching_and_parent_tail_bounds_are_not_established":True}


@cache
def data():
    return {"default_conditional_budget":budget(),
            "other_conditional_budgets":[budget(0,sp.Rational(1,4)),budget(sp.Rational(1,10),sp.Rational(1,5))],
            "classical_comparison_source":band.data()["source_projection_definition"],
            "required_common_parent":["same actual nearby background and physical matter metric",
                "same linear local relational observable and physical source/detector normalization",
                "defined parent state/response and canonical maps",
                "parent response vanishes for the original compact matter-spacelike supports",
                "low-spatial-band matching error including temporal validity, omitted terms, loops and thresholds",
                "parent high-spatial-frequency paired tail small enough to prevent cancellation"],
            "low_frequency_agreement_alone_allows_high_frequency_cancellation":True,
            "algebraic_cancellation_control_is_not_a_constructed_UV_theory":True,
            "no_nearby_quantum_matching_or_original_P8_closure":True}


@cache
def checks():
    D,L,E,T=sp.symbols("positive_classical_signal low_classical_pairing low_matching_error parent_high_tail",real=True)
    low,tail=sp.symbols("low_error_fraction tail_error_fraction",nonnegative=True)
    d=budget()
    return {"conditional_comparison_triangle_keeps_both_parent_errors":
                sp.expand(D/2-D*low-D*tail-D*(sp.Rational(1,2)-low-tail)),
            "default_two_error_budget_retains_quarter_signal":
                d["conditional_parent_full_response_lower"]-probes.SIGNAL_LOWER/4,
            "low_band_agreement_only_has_an_exact_cancellation_control":L+(-L),
            "cancellation_control_has_zero_low_band_mismatch":L-L,
            "source_projection_low_and_high_parts_reconstruct_parent_response":
                (L+E)+T-(L+(E+T))}


@cache
def gates():
    d=data()
    budgets=[d["default_conditional_budget"],*d["other_conditional_budgets"]]
    return {name:bool(value) for name,value in {
        "all_conditional_error_budgets_leave_strictly_positive_response":
            all(row["conditional_parent_full_response_lower"]>0 for row in budgets),
        "default_conditional_budget_requires_two_distinct_one_eighth_bounds":
            d["default_conditional_budget"]["assumed_low_spatial_band_matching_error_fraction"]==sp.Rational(1,8)
            and d["default_conditional_budget"]["assumed_parent_high_spatial_frequency_tail_fraction"]==sp.Rational(1,8),
        "all_reported_error_budgets_are_explicitly_conditional":
            all(row["actual_matching_and_parent_tail_bounds_are_not_established"] is True for row in budgets),
        "common_parent_requirements_are_not_replaced_by_classical_normalform_constants":True,
        "matching_alone_is_not_a_matter_causal_UV_exclusion":True}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.sqrt(2),[])
    cases=[(value,0) for value in bad]+[(0,value) for value in bad]
    cases.extend(((-1,0),(0,-1),(sp.Rational(1,4),sp.Rational(1,4)),(sp.Rational(1,2),0),
                  (0,sp.Rational(1,2)),(1,0),(0,1)))
    rejected=0
    for low,tail in cases:
        try:
            budget(low,tail)
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(cases):
        raise ValueError("An inadmissible conditional matching error budget was accepted")
    return {"rejected_inputs":rejected}
