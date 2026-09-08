"""Independent D-dimensional physical-Hamiltonian source variation."""
from functools import cache

import sympy as sp

from . import jets, variation


@cache
def checks():
    a, q, m2 = sp.symbols("scale physical_momentum_squared mass_squared", positive=True)
    e = sp.Symbol("source_parameter", real=True)
    am, bm = 1+e*jets.alpha[0]*jets.n[0], 1+e*jets.beta[0]*jets.n[0]
    out = {}
    for sector in ("T", "L"):
        if sector == "T":
            g2, bare = a**(jets.D-2), q+m2*bm
            A, B = sp.Integer(0), -m2*jets.beta[0]/bare
        else:
            g2, bare = a**jets.D*m2*am*q/(q+m2*am), bm*(q/am+m2)
            A, B = jets.alpha[0]*q/(am*(q+m2*am)), -jets.beta[0]/bm

        def project(value):
            return sp.factor(value.subs(q, m2*jets.z/(1-jets.z)))

        rate = (jets.time(g2)+a*jets.H[0]*sp.diff(g2, a)-2*jets.H[0]*q*sp.diff(g2, q))/(2*g2)
        targets = variation.data(sector)
        out[sector+"_independent_delta_log_g"] = sp.factor(project(sp.diff(g2, e).subs(e, 0)/(2*g2.subs(e, 0)))-targets["delta_log_g"])
        out[sector+"_independent_delta_log_frequency"] = sp.factor(project(sp.diff(bare, e).subs(e, 0)/(2*bare.subs(e, 0)))-targets["delta_log_frequency"])
        out[sector+"_independent_delta_rate"] = sp.factor(project(sp.diff(rate, e).subs(e, 0))-targets["delta_rate"])
        for name, value in (("A", A), ("B", B)):
            out[sector+"_independent_source_"+name] = sp.factor(project(value.subs(e, 0))-targets["source_"+name])
            out[sector+"_independent_delta_source_"+name] = sp.factor(project(sp.diff(value, e).subs(e, 0))-targets["delta_source_"+name])
    return out
