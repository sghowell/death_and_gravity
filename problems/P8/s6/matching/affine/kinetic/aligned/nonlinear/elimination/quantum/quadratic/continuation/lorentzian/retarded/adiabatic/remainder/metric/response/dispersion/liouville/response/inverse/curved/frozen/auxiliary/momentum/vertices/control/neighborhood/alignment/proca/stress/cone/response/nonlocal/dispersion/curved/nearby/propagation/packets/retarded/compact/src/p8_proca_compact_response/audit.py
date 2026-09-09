"""Quantitative compact classical response with explicit conditional UV boundaries."""
from functools import cache

from . import band, bumps, matching, probes


@cache
def residuals():
    out={}
    for module in (bumps,probes,band,matching):
        rows=module.checks()
        if set(out).intersection(rows):
            raise ValueError("Repeated compact classical response residual")
        out.update(rows)
    return out


@cache
def gates():
    out={}
    for module in (probes,band,matching):
        rows=module.gates()
        if set(out).intersection(rows):
            raise ValueError("Repeated compact classical response gate")
        out.update(rows)
    return out


def controls():
    groups={"compact_bump":bumps.controls(),"conditional_matching":matching.controls()}
    return {"rejected_inputs":sum(group["rejected_inputs"] for group in groups.values()),
            "groups":groups,"no_parent_source_or_scientific_library_changed":True}
