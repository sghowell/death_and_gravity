"""Regular active-component physical-clock identity with independent NEC sources."""

from functools import cache

import sympy as sp
from p8_star.exact import nonnegative, positive, rational


@cache
def derive_three():
    gs = sp.symbols("G0:3", positive=True)
    ys = sp.symbols("y0:3", positive=True)
    cs = sp.symbols("c0:3", positive=True)
    nulls = sp.symbols("null0:3", nonnegative=True)
    h, hp = sp.symbols("H Hprime", real=True)
    yp = tuple(y*(c-1)*h for y, c in zip(ys, cs, strict=True))
    k = sum(g*y**2 for g, y in zip(gs, ys, strict=True))
    kp = sum(sp.diff(k, y)*rate for y, rate in zip(ys, yp, strict=True))
    weighted_null = sum(c*y**4*n for c, y, n in zip(cs, ys, nulls, strict=True))
    # Interaction null densities are cancelled by the literal potential sum.
    geometric = sum(g*y**2*(hp+(1-c)*h**2) for g, y, c in zip(gs, ys, cs, strict=True))
    return {"Gs": gs, "ys": ys, "cs": cs, "nulls": nulls, "H": h, "Hprime": hp,
            "yprimes": yp, "K": k, "Kprime": kp, "weighted_null": weighted_null,
            "weighted_geometric_null": geometric,
            "combined_residual": k*hp-h*kp/2+weighted_null/2}


def reconstruct(Gs, ys, cs, H, nulls, *, physical=0):
    """Conditional regular-sum rates for precisely the supplied active component.

    The physical metric must be included with y_r=c_r=1. This function does
    not verify Bianchi, sum invertibility or a background solution. Zero-beta
    disconnected metrics must be omitted; use disconnected_GR for them.
    """
    lists = (Gs, ys, cs, nulls)
    if not all(isinstance(values, (tuple, list)) for values in lists):
        raise TypeError("The finite component inputs must be lists or tuples")
    if len({len(values) for values in lists}) != 1 or not Gs:
        raise ValueError("A nonempty active component needs equal input lengths")
    if isinstance(physical, bool) or not isinstance(physical, int):
        raise TypeError("physical must be an integer component index")
    if not 0 <= physical < len(Gs):
        raise ValueError("physical must index the supplied component")
    gs = tuple(positive(value, "G_i") for value in Gs)
    ratios = tuple(positive(value, "y_i") for value in ys)
    cones = tuple(positive(value, "c_i") for value in cs)
    source_nulls = tuple(nonnegative(value, "rho_i+p_i") for value in nulls)
    h = rational(H)
    if ratios[physical] != 1 or cones[physical] != 1:
        raise ValueError("The physical proper-clock normalization is y_r=c_r=1")
    k = sum(g*y**2 for g, y in zip(gs, ratios, strict=True))
    yp = tuple(y*(c-1)*h for y, c in zip(ratios, cones, strict=True))
    kp = sum(2*g*y*rate for g, y, rate in zip(gs, ratios, yp, strict=True))
    null = sum(c*y**4*n for c, y, n in zip(cones, ratios, source_nulls, strict=True))
    hp = (h*kp-null)/(2*k)
    return {"K": k, "Kprime": kp, "Hprime": hp, "yprimes": yp,
            "weighted_null": null, "scaled_H_prime": -null/(2*k**sp.Rational(3, 2))}


def disconnected_GR(G, null):
    """The independently coupled beta_r=0 metric obeys flat GR, not sum locking."""
    g, nh = positive(G, "G_r"), nonnegative(null, "rho_r+p_r")
    return {"Hprime": -nh/(2*g)}


def checks():
    d = derive_three()
    y, c, g = sp.symbols("y c G", positive=True)
    h, hp = sp.symbols("H Hprime", real=True)
    yp = y*(c-1)*h
    leaf_derivative = (hp/y-h*yp/y**2)/(c*y)
    k = sp.Symbol("K", positive=True)
    kp = sp.Symbol("Kprime", real=True)
    return {
        "proper_clock_leaf_null": sp.factor(-2*g*leaf_derivative*c*y**4
                                              +2*g*y**2*(hp+(1-c)*h**2)),
        "K_chain_rule": sp.expand(d["Kprime"]-2*d["H"]*sum(
            gi*yi**2*(ci-1) for gi, yi, ci in zip(d["Gs"], d["ys"], d["cs"], strict=True))),
        "weighted_null_combination": sp.expand(d["weighted_geometric_null"]
                                               -d["K"]*d["Hprime"]+d["H"]*d["Kprime"]/2),
        "normalized_monotonicity": sp.simplify(hp/sp.sqrt(k)-h*kp/(2*k**sp.Rational(3, 2))
                                               -(k*hp-h*kp/2)/k**sp.Rational(3, 2)),
        "bounce_without_H_division": sp.factor(d["combined_residual"].subs(d["H"], 0)
                                               -d["K"]*d["Hprime"]-d["weighted_null"]/2),
    }
