"""Literal geometric, Euler, rational-evolution and continuous-cone audit."""
from functools import cache

from . import bounds, geometry, principal, rational


@cache
def residuals():
    out={}
    for group in (geometry.checks(),principal.checks(),rational.checks(),
                  rational.chain_checks(),bounds.checks()):
        if set(out).intersection(group):
            raise ValueError("Repeated actual nearby-cone identity")
        out.update(group)
    return out


def gates():
    return bounds.gates()


def controls():
    first,second=geometry.controls(),bounds.controls()
    return {"rejected_inputs":first["rejected_inputs"]+second["rejected_inputs"],
            "geometric_controls":first,"whole_box_controls":second,
            "no_parent_action_state_profile_or_scientific_library_changed":True}
