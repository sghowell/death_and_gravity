"""Exact central jets and explicit analytic local-solution proof gates."""
from functools import cache

from . import existence, jets, taylor


@cache
def residuals():
    out={}
    for group in (jets.checks(),existence.checks(),taylor.checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated nearby-bounce residual")
        out.update(group)
    return out


def gates():
    first,second=existence.gates(),taylor.gates()
    if set(first).intersection(second):
        raise ValueError("Repeated nearby-bounce gate")
    return {**first,**second}


def controls():
    first,second=existence.controls(),taylor.controls()
    return {"rejected_inputs":first["rejected_inputs"]+second["rejected_inputs"],
            "coarse_domain_controls":first,"enlarged_domain_controls":second,
            "no_parent_action_profile_state_or_scientific_function_changed":True}
