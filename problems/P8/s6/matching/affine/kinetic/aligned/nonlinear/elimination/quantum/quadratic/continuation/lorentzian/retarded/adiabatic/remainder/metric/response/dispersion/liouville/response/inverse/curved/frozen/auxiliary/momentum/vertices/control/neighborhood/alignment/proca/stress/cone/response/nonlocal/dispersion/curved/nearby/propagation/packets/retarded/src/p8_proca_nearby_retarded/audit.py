"""Complete exact classical retarded-source audit without a UV-front transfer."""
from functools import cache

from . import fronts, remainder, source, support, transport


@cache
def residuals():
    out={}
    for module in (source,transport,fronts,remainder,support):
        rows=module.checks()
        if set(out).intersection(rows):
            raise ValueError("Repeated nearby retarded-source residual")
        out.update(rows)
    return out


@cache
def gates():
    out={}
    for module in (transport,fronts,remainder,support):
        rows=module.gates()
        if set(out).intersection(rows):
            raise ValueError("Repeated nearby retarded-source gate")
        out.update(rows)
    return out


def controls():
    groups={"high_frequency_domain":transport.controls(),"compact_probe_times":support.controls()}
    return {"rejected_inputs":sum(group["rejected_inputs"] for group in groups.values()),
            "groups":groups,"no_parent_action_state_profile_or_scientific_library_changed":True}
