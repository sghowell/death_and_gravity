"""Frozen oscillator-frame tree symbols; no adiabatic or chart-invariance claim."""

from functools import cache

import sympy as sp
from p8_control import oscillator
from p8_physical.vertices import nonexceptional, tensor_basis

from .frequency import PARTITIONS, Leg, reverse, vertex


def omega_squared(kind, wave, time_point, chart):
    out = oscillator.derive("tensor" if kind == "t" else chart)["omega_squared"]
    q = sum(sp.Rational(k)**2 for k in wave)
    out = out.subs(oscillator.q, q)
    return sp.factor(out if time_point is None else out.subs(oscillator.x, time_point))


@cache
def tree(legs, time_point, chart="gamma"):
    if len(legs) != 4:
        raise ValueError("Four external legs required")
    nonexceptional(legs)
    if sp.simplify(sum(sp.sympify(leg.frequency) for leg in legs)) != 0:
        raise ValueError("Frozen symbol requires total incoming frequency zero")
    contact_data = vertex(legs, time_point, chart)
    channels = []
    for left, right in PARTITIONS:
        wave = tuple(-sum(sp.Rational(legs[i].wave[j]) for i in left) for j in range(3))
        energy = -sum(sp.sympify(legs[i].frequency) for i in left)
        entries = []
        for kind, E in (("s", None), *(("t", E) for E in tensor_basis(wave))):
            internal = Leg(wave, kind, energy, E)
            inverse_propagator = sp.factor(energy**2-omega_squared(kind, wave, time_point, chart))
            if inverse_propagator == 0:
                raise ValueError("Internal propagator pole; specify a regulated observable")
            Vleft = vertex(tuple(legs[i] for i in left)+(internal,), time_point, chart)["kernel"]
            Vright = vertex(tuple(legs[i] for i in right)+(reverse(internal),), time_point, chart)["kernel"]
            contribution = sp.factor(-Vleft*Vright/inverse_propagator)
            entries.append({"kind": kind, "polarization": E, "inverse_propagator": inverse_propagator,
                            "left_vertex": Vleft, "right_vertex": Vright, "contribution": contribution})
        channels.append({"partition": (left, right), "entries": entries})
    total = sp.factor(contact_data["kernel"]+sum(e["contribution"] for c in channels for e in c["entries"]))
    shells = [sp.factor(sp.sympify(leg.frequency)**2-omega_squared(leg.kind, leg.wave, time_point, chart)) for leg in legs]
    return {"total": total, "contact": contact_data["kernel"], "contact_details": contact_data,
            "channels": channels, "external_shell_residuals": shells,
            "scope": "frozen oscillator-frame tree symbol, not a time-dependent physical amplitude"}


def free_redefinition_check():
    """Free massive chi=phi+g*phi^2/2: derivative contacts and exchange cancel."""
    g, mass2, s, t = sp.symbols("g mass2 s t")
    u = 4*mass2-s-t
    contact = g*g*mass2
    exchange = [-g*g*(channel-mass2) for channel in (s, t, u)]
    return {"contact": contact, "exchange": exchange, "residual": sp.expand(contact+sum(exchange))}
