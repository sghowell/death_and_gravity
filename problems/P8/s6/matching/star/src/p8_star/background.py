"""Dynamic-stratum identity; the global no-crossing argument is in the proof."""

from functools import cache

import sympy as sp

from .exact import nonnegative, positive, rational


@cache
def derive_two_leaf():
    g0 = sp.Symbol("G_u", nonnegative=True)
    g1, g2, r1, r2, c1, c2 = sp.symbols("G_1 G_2 R_1 R_2 c_1 c_2", positive=True)
    h, hp, j1, j2 = sp.symbols("H_u H_u_prime J_1 J_2", real=True)
    nh = sp.Symbol("n_h", nonnegative=True)
    ks = (g1/r1**2, g2/r2**2)
    q = ks[0]*(1-c1)+ks[1]*(1-c2)
    k = g0+sum(ks)
    kp = -2*h*q
    leaf1 = ks[0]*(hp+(1-c1)*h**2)-j1*(1-c1)/r1**3
    leaf2 = ks[1]*(hp+(1-c2)*h**2)-j2*(1-c2)/r2**3
    central = g0*hp+nh/2+j1*(1-c1)/r1**3+j2*(1-c2)/r2**3
    return {"G_u": g0, "Gs": (g1, g2), "Rs": (r1, r2), "cs": (c1, c2),
            "Js": (j1, j2), "H": h, "Hprime": hp, "nh": nh,
            "K": k, "Kprime": kp, "weighted_cone_sum": q,
            "leaf_residuals": (leaf1, leaf2), "central_residual": central,
            "combined_residual": k*hp-h*kp/2+nh/2}


def dynamic_reconstruction(Gs, Rs, cs, H, nh, *, central_G=0):
    """Conditional exact rates when precisely the supplied legs are dynamic.

    This does not test the field equations or decide that any supplied link has
    J!=0. Algebraic and endpoint-only legs must be omitted. Empty lists are
    allowed only for central_G>0, giving ordinary central Einstein evolution.
    """
    if not all(isinstance(values, (list, tuple)) for values in (Gs, Rs, cs)):
        raise TypeError("The finite leaf inputs must be lists or tuples")
    if not len(Gs) == len(Rs) == len(cs):
        raise ValueError("Einstein coefficients, ratios and cones must have equal length")
    gs = tuple(positive(value, "G_i") for value in Gs)
    rs = tuple(positive(value, "R_i") for value in Rs)
    cones = tuple(positive(value, "c_i") for value in cs)
    g0 = nonnegative(central_G, "G_u")
    h, null = rational(H), nonnegative(nh, "n_h")
    if not gs and g0 == 0:
        raise ValueError("No positive kinetic coefficient remains in this dynamic stratum")
    ks = tuple(g/r**2 for g, r in zip(gs, rs, strict=True))
    k = g0+sum(ks)
    q = sum(ki*(1-c) for ki, c in zip(ks, cones, strict=True))
    kp = -2*h*q
    hp = (-null/2-h**2*q)/k
    return {"K": k, "Kprime": kp, "Hprime": hp,
            "scaled_H_prime": -null/(2*k**sp.Rational(3, 2)),
            "Rprimes": tuple(r*(1-c)*h for r, c in zip(rs, cones, strict=True))}


def compact_comparison_bound(r_log_rate_bounds, c_log_rate_bounds):
    """Exact envelope for the written positive-part/Gronwall argument.

    Each entry is an externally justified supremum of |R_i'/R_i| or
    |c_i'/c_i| on the same compact regular interval. This interface does not
    estimate those suprema for an unspecified solution.
    """
    if not all(isinstance(values, (list, tuple)) for values in
               (r_log_rate_bounds, c_log_rate_bounds)):
        raise TypeError("Rate bounds must be finite lists or tuples")
    if len(r_log_rate_bounds) != len(c_log_rate_bounds):
        raise ValueError("The same genuine leaves need both rate bounds")
    bounds = [nonnegative(value, "logarithmic-rate bound")
              for values in (r_log_rate_bounds, c_log_rate_bounds) for value in values]
    return max(bounds, default=sp.S.Zero)


def checks():
    d = derive_two_leaf()
    r, c, g, h, hp, j = sp.symbols("R c G H Hprime J", nonzero=True)
    rp = r*(1-c)*h
    # D_tau_i=(R/c)D_T on the actual physical-clock chart.
    ray = g*(r/c)*(rp*h+r*hp)-r*j*(1/c-1)
    normalized = g/r**2*(hp+(1-c)*h**2)-j/r**3*(1-c)
    k = sp.Symbol("K", positive=True)
    kp = sp.Symbol("Kprime", real=True)
    chain_kprime = sum(sp.diff(d["K"], ri)*ri*(1-ci)*d["H"]
                       for ri, ci in zip(d["Rs"], d["cs"], strict=True))
    return {
        "single_leaf_clock_and_raychaudhuri": sp.cancel(ray*c/r**4-normalized),
        "two_leaf_central_null_cancellation": sp.factor(
            sum(d["leaf_residuals"])+d["central_residual"]-d["combined_residual"]),
        "scaled_H_derivative": sp.simplify(
            hp/sp.sqrt(k)-h*kp/(2*k**sp.Rational(3, 2))-(k*hp-h*kp/2)/k**sp.Rational(3, 2)),
        "K_derivative": sp.expand(chain_kprime-d["Kprime"]),
        "bounce_without_H_division": sp.factor(d["combined_residual"].subs(d["H"], 0)
                                                -d["K"]*d["Hprime"]-d["nh"]/2),
    }
