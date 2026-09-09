"""Native exact residual aggregation and explicit continuous-domain gates."""
from functools import cache

import sympy as sp

from . import (
    canonical,
    covariance,
    domains,
    green,
    jets,
    recursion,
    riccati,
    transitions,
)


@cache
def residuals():
    out={}
    for module in (canonical,green,riccati,jets,covariance,domains,transitions,recursion):
        values=module.checks()
        if set(out)&set(values):
            raise ValueError("An identity name was reused")
        out.update(values)
    out.update(recursion.mode_checks())
    return out


@cache
def gates():
    out={}
    for module in (domains,jets,transitions,recursion):
        values=module.gates()
        if set(out)&set(values):
            raise ValueError("A proof gate name was reused")
        out.update(values)
    d=riccati.principal()
    out.update({
        "leading_clock_kinetic_parameter_positive":d["clock_kinetic"].is_positive is True,
        "leading_covariance_has_positive_determinant":bool(
            sp.factor(d["positive_graph_imaginary_part"].det()).is_positive is True),
        "all_three_Riccati_frequency_sum_denominators_positive":all(
            value.is_positive is True for value in d["sum_denominators"]),
        "cutoff_endpoints_are_both_accepted":covariance.cutoff_weight(0)==0 and covariance.cutoff_weight(1)==1,
    })
    return {name:bool(value) for name,value in out.items()}


@cache
def controls():
    return {"rejected_inputs":sum(module.controls()["rejected_inputs"]
        for module in (domains,riccati,covariance))}
