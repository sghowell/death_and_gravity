"""Literal new action, fresh actual solution and principal cone audit."""
from functools import cache

from . import central, domain, flow, model, original


@cache
def residuals():
    out={}
    for module in (model,domain,flow,original,central):
        rows=module.checks()
        if set(rows).intersection(out):
            raise ValueError("Duplicate retuned-margin residual")
        out.update(rows)
    return out


@cache
def gates():
    out={}
    for module in (domain,flow,original,central):
        rows=module.gates()
        if set(rows).intersection(out):
            raise ValueError("Duplicate retuned-margin proof gate")
        out.update(rows)
    return out


def controls():
    groups={"new_actual_solution_duration":flow.controls(),"central_family":central.controls()}
    return {"groups":groups,"rejected_inputs":sum(row["rejected_inputs"] for row in groups.values()),
        "no_frozen_action_report_background_state_or_scientific_library_changed":True}
