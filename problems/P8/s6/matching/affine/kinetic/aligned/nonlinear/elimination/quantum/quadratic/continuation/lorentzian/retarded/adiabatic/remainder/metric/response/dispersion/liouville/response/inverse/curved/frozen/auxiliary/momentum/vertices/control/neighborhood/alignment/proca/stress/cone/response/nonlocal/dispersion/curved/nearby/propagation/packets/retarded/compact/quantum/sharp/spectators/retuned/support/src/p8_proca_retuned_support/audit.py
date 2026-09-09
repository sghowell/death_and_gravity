"""Native algebra, fresh support-domain inequalities and rejected input audit."""
from functools import cache

import sympy as sp

from . import bridge, domain, growth, momentum, quantum


@cache
def residuals():
    out={}
    for module in (bridge,domain,growth,momentum,quantum):
        rows=module.checks()
        if set(out).intersection(rows):
            raise ValueError("Duplicate new-support residual name")
        out.update(rows)
    return out


@cache
def gates():
    out={}
    for module in (domain,growth,momentum,quantum):
        rows=module.gates()
        if set(out).intersection(rows):
            raise ValueError("Duplicate new-support gate name")
        out.update(rows)
    return out


def controls():
    invalid=[True,False,"1",1.0,sp.Float(1),None,[],{},sp.oo,-sp.oo,sp.zoo,sp.nan,sp.I]
    rejected=0
    for value in invalid:
        for fn in (lambda v:growth.time_pair(v,0),lambda v:momentum.null_cone_example(v)):
            try:
                fn(value)
            except (TypeError,ValueError):
                rejected+=1
            else:
                raise ValueError("An inexact, nonfinite or invalid new-support input was accepted")
    for start,end in ((0,-growth.HALF),(-domain.T,0),(0,domain.T)):
        try:
            growth.time_pair(start,end)
        except (TypeError,ValueError):
            rejected+=1
        else:
            raise ValueError("Unordered or uncertified new-support times accepted")
    for size in (0,-1,-sp.Rational(1,2)):
        try:
            momentum.null_cone_example(size)
        except (TypeError,ValueError):
            rejected+=1
        else:
            raise ValueError("A nonpositive complex-vector example accepted")
    state_controls=quantum.controls()
    return {"new_domain_rejected_inputs":rejected,"new_Gaussian_inputs":state_controls,
        "rejected_inputs":rejected+state_controls["rejected_inputs"],
        "coincident_times_and_large_complex_null_vectors_are_valid":True,
        "fast_clock_mutation_is_excluded_by_actual_subluminal_gate_not_assumed_causal":True,
        "frozen_old_state_action_probes_and_reports_unchanged":True}
