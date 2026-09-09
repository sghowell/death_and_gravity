"""Finite principal-extension and ordinary-Proca spectator screening audit."""
from functools import cache

from . import pencil, proca, repair


@cache
def residuals():
    out={}
    for module in (pencil,repair,proca):
        rows=module.checks()
        if set(out).intersection(rows):
            raise ValueError("Repeated spectator-repair identity")
        out.update(rows)
    return out


@cache
def gates():
    out={}
    for module in (pencil,repair,proca):
        rows=module.gates()
        if set(out).intersection(rows):
            raise ValueError("Repeated spectator-repair gate")
        out.update(rows)
    return out


def controls():
    groups={"relative_norm_budget":repair.controls(),"finite_species":proca.controls()}
    return {"groups":groups,"rejected_inputs":sum(v["rejected_inputs"] for v in groups.values()),
        "no_old_action_background_state_or_scientific_library_changed":True}
