"""The same mode-clock transformation before taking the dimensional limit."""
from functools import cache

import sympy as sp
from p8_vector_liouville import reduction
from p8_vector_metric_local import canonical, jets

D = jets.D


@cache
def pumps():
    a, am, bm, k, m = reduction.a, reduction.am, reduction.bm, reduction.k, reduction.m
    U = a**2*m**2*am
    return {"T": a**(D-3), "L": m**2*a**(D-1)*sp.sqrt(am*bm)*k**2/(k**2+U)}


@cache
def checks():
    a, N, am, k, m = reduction.a, reduction.N, reduction.am, reduction.k, reduction.m
    q = k**2/a**2
    kinetic = {"T": a**(D-2)/N, "L": a**D*m**2*am*q/(N*(q+m**2*am))}
    out = {}
    for sector in ("T", "L"):
        target = pumps()[sector]
        out[sector+"_full_dimensional_acoustic_pump"] = sp.simplify(
            kinetic[sector]*reduction.data(sector)["acoustic_rate"]-target)
        out[sector+"_physical_dimensional_limit"] = sp.simplify(
            target.subs(D, 3)-reduction.data(sector)["acoustic_pump_squared"])
        out[sector+"_pump_dimensional_first_jet"] = sp.simplify(
            sp.diff(target, D).subs(D, 3)-sp.log(a)*target.subs(D, 3))
    n, v, z = jets.n[0], jets.v[0], jets.z
    alpha, beta = jets.alpha[0], jets.beta[0]
    rate_sources = {"T": n-v, "L": (1+(beta-alpha)/2)*n-v}
    changed = {"T": (D-3)*v/2,
               "L": (D-1)*v/2+(alpha+beta)*n/4-(1-z)*(2*v+alpha*n)/2}
    for sector in ("T", "L"):
        out[sector+"_dimensional_source_normalization"] = sp.expand(
            canonical.data(sector)["delta_log_g"]+rate_sources[sector]/2-changed[sector])
        out[sector+"_source_dimensional_first_jet"] = sp.diff(changed[sector], D)-v/2
    return out
