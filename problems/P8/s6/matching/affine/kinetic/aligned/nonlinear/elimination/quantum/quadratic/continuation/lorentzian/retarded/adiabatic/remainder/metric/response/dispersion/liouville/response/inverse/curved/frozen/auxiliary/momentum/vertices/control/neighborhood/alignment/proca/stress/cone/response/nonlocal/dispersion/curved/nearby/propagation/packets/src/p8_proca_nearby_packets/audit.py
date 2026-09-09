"""Actual finite-frequency packet audit; no interacting cutoff or UV transfer."""
from functools import cache

from . import domains, finite, modes, normal_form, observable, packets


@cache
def residuals():
    out={}
    for module in (finite,modes,domains,normal_form,observable,packets):
        group=module.checks()
        if set(out).intersection(group):
            raise ValueError("Repeated actual nearby packet residual")
        out.update(group)
    return out


@cache
def gates():
    out={}
    for module in (domains,normal_form,observable,packets):
        group=module.gates()
        if set(out).intersection(group):
            raise ValueError("Repeated actual nearby packet gate")
        out.update(group)
    return out


def controls():
    groups={"complex_domain":domains.controls(),"normal_form":normal_form.controls(),"packets":packets.controls()}
    return {"rejected_inputs":sum(group["rejected_inputs"] for group in groups.values()),
            "groups":groups,"no_parent_action_state_profile_or_scientific_library_changed":True}
