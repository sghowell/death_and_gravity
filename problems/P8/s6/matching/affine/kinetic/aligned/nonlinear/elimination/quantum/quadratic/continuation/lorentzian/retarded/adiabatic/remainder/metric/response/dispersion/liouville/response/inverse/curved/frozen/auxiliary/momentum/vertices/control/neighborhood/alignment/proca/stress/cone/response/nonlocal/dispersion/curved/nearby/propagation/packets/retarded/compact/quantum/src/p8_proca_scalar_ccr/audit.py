"""Actual reduced quadratic scalar CCR and positive Gaussian linear-field audit."""
from functools import cache

from . import canonical, growth, response, state


@cache
def residuals():
    out={}
    for module in (canonical,growth,state,response):
        rows=module.checks()
        if set(out).intersection(rows):
            raise ValueError("Repeated exact scalar CCR residual")
        out.update(rows)
    return out


@cache
def gates():
    out={}
    for module in (growth,response):
        rows=module.gates()
        if set(out).intersection(rows):
            raise ValueError("Repeated scalar CCR proof gate")
        out.update(rows)
    return out


def controls():
    return {"initial_Gaussian_weight":state.controls(),
            "rejected_inputs":state.controls()["rejected_inputs"],
            "no_parent_source_or_scientific_library_changed":True}
