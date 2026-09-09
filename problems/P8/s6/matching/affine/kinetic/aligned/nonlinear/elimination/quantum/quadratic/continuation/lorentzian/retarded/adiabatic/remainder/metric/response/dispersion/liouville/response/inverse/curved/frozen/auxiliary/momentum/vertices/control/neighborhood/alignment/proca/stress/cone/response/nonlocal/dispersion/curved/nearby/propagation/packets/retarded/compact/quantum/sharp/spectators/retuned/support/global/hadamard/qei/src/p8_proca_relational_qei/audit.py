"""Native algebra and continuous-domain inputs for relational quadratic bounds."""
from functools import cache

from . import modes, observers, qei, stress


@cache
def residuals():
    out={}
    for module in (modes,observers,qei,stress):
        rows=module.checks()
        if set(out)&set(rows):
            raise ValueError("An exact identity name was reused")
        out.update(rows)
    return out


@cache
def gates():
    out={}
    for module in (modes,observers,qei,stress):
        rows=module.gates()
        if set(out)&set(rows):
            raise ValueError("A proof gate name was reused")
        out.update(rows)
    return {name:bool(value) for name,value in out.items()}


@cache
def controls():
    return {"rejected_inputs":observers.controls()["rejected_inputs"]+qei.controls()["rejected_inputs"]}
