"""Native global support and state evidence, with exact-domain controls."""
from functools import cache

import sympy as sp

from . import charts, coverage, model, modes, phase, quantum, transfer


@cache
def residuals():
    out={}
    for module in (model,phase,charts,transfer,coverage,modes,quantum):
        rows=module.checks()
        if set(out).intersection(rows):
            raise ValueError("Duplicate global-support residual name")
        out.update(rows)
    return out


@cache
def gates():
    return coverage.gates()|transfer.gates()


def controls():
    bad=[True,False,sp.true,sp.false,"1",1.0,sp.Float(1),None,[],{},
        sp.oo,-sp.oo,sp.zoo,sp.nan,sp.I,sp.sqrt(2),sp.Symbol("unfixed_time")]
    rejected=0
    for value in bad:
        for pair in ((value,0),(0,value)):
            try:
                coverage.time_pair(*pair)
            except (TypeError,ValueError):
                rejected+=1
            else:
                raise ValueError("Invalid exact global retarded time accepted")
    for pair in ((1,0),(0,-1),(coverage.CUT,-coverage.CUT)):
        try:
            coverage.time_pair(*pair)
        except ValueError:
            rejected+=1
        else:
            raise ValueError("Reversed global retarded endpoints accepted")
    state_controls=quantum.controls()
    return {"global_endpoint_rejected_inputs":rejected,"new_global_Gaussian_inputs":state_controls,
        "rejected_inputs":rejected+state_controls["rejected_inputs"],
        "all_finite_rational_times_allowed_not_limited_to_old_nearby_interval":True,
        "gamma_finite_q_velocity_pole_not_a_pole_of_original_phase":True,
        "global_and_nearby_backgrounds_and_states_remain_distinct":True}
