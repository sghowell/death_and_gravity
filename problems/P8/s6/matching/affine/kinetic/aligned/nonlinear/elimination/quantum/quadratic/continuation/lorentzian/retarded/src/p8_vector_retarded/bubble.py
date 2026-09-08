"""Unequal-time mass bubble for general external spatial momentum."""
from functools import cache

import sympy as sp

from . import modes

alpha = sp.symbols("alpha_0:2", real=True)
beta = sp.symbols("beta_0:2", real=True)
CHANNELS = ("TT", "TL", "LT", "LL")


def vertex(time):
    time = modes.label(time)
    return sp.diag(modes.mass**2*alpha[time], *([-modes.mass**2*beta[time]]*3))


def channel(value):
    if type(value) is not str or value not in CHANNELS:
        raise ValueError("Require one of TT, TL, LT, LL")
    return value


def direct(sector):
    return _direct(channel(sector))


@cache
def _direct(sector):
    first = (0, 1) if sector[0] == "T" else (2,)
    second = (0, 1) if sector[1] == "T" else (2,)
    out = 0
    for i in first:
        for j in second:
            pair = [(modes.polarization(time, 0)[i].T*vertex(time)*
                     modes.polarization(time, 1)[j])[0] for time in (0, 1)]
            out += pair[0]*sp.conjugate(pair[1])/2
    return modes.angular(out)


def closed(sector):
    return _closed(channel(sector))


@cache
def _closed(sector):
    x, y = [modes.data(time, 0) for time in (0, 1)]
    r, s = [modes.data(time, 1) for time in (0, 1)]
    c, m = modes.cosine, modes.mass
    if sector == "TT":
        return m**4*beta[0]*beta[1]*(1+c**2)*x["vT"]*r["vT"]*sp.conjugate(y["vT"]*s["vT"])/2
    if sector == "TL":
        return m**2*beta[0]*beta[1]*(1-c**2)*r["omega"]*s["omega"]*x["vT"]*r["vL"]*sp.conjugate(y["vT"]*s["vL"])/2
    if sector == "LT":
        return m**2*beta[0]*beta[1]*(1-c**2)*x["omega"]*y["omega"]*x["vL"]*r["vT"]*sp.conjugate(y["vL"]*s["vT"])/2
    pair = []
    for time, one, two in ((0, x, r), (1, y, s)):
        pair.append(alpha[time]*one["k"]*two["k"]*one["pL"]*two["pL"]/(one["omega"]*two["omega"])
                    +beta[time]*one["omega"]*two["omega"]*c*one["vL"]*two["vL"])
    return pair[0]*sp.conjugate(pair[1])/2


@cache
def checks():
    out = {sector+"_explicit_polarization_sum": modes.angular(direct(sector)-closed(sector))
           for sector in CHANNELS}
    W, X = modes.wightman(0), modes.wightman(1)
    V, Z = vertex(0), vertex(1)
    index_sum = sum(V[i, j]*Z[r, s]*W[i, r]*X[j, s]
                    for i in range(4) for j in range(4) for r in range(4) for s in range(4))/2
    out["two_Wick_pairings_and_covariant_index_sum"] = modes.angular(index_sum-sum(closed(sector) for sector in CHANNELS))
    for sector in ("TL", "LT"):
        out[sector+"_zero_external_momentum_vanishes"] = sp.expand(closed(sector).subs(modes.cosine, -1))
    x = modes.data(0, 0)
    pair = alpha[0]*x["k"]**2*x["pL"]**2/x["omega"]**2-beta[0]*x["omega"]**2*x["vL"]**2
    explicit = (modes.polarization(0, 0)[2].T*vertex(0)*modes.polarization(0, 1)[2])[0]
    replacement = {value: x[key] for key, value in modes.data(0, 1).items()}
    out["zero_external_longitudinal_canonical_pair_vertex"] = sp.factor(
        explicit.subs(replacement).subs({modes.cosine: -1, modes.sine: 0})-pair)
    return out
