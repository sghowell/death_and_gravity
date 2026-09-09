"""Sharper actual complex jets, complete transfer and newly named compact probes."""
from functools import cache

from . import basis, jets, probes, transfer


@cache
def residuals():
    out={}
    for module in (jets,basis,transfer,probes):
        rows=module.checks()
        if set(out).intersection(rows):
            raise ValueError("Repeated sharper response identity")
        out.update(rows)
    return out


@cache
def gates():
    out={}
    for module in (jets,basis,transfer,probes):
        rows=module.gates()
        if set(out).intersection(rows):
            raise ValueError("Repeated sharper response gate")
        out.update(rows)
    return out


def controls():
    groups={"complex_jets":jets.controls(),"frequency_domain":transfer.controls()}
    return {"groups":groups,"rejected_inputs":sum(v["rejected_inputs"] for v in groups.values()),
            "no_parent_source_scientific_library_or_old_probe_changed":True}
