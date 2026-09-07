"""Literal aligned-coframe variations; no inverse of a summed coframe is used."""

from functools import cache

import sympy as sp
from p8_star.exact import positive, rational


@cache
def derive_three():
    lam = sp.Symbol("lambda", positive=True)
    aa = sp.symbols("a0:3", positive=True)
    nn = sp.symbols("n0:3", positive=True)
    beta = sp.symbols("beta0:3", real=True)
    velocities = sp.symbols("v0:3", real=True)
    sa = sum(b*a for b, a in zip(beta, aa, strict=True))
    sn = sum(b*n for b, n in zip(beta, nn, strict=True))
    sadot = sum(b*v for b, v in zip(beta, velocities, strict=True))
    lag = -lam*sn*sa**3
    rho = tuple(sp.factor(-sp.diff(lag, n)/a**3) for n, a in zip(nn, aa, strict=True))
    pressure = tuple(sp.factor(sp.diff(lag, a)/(3*n*a**2)) for n, a in zip(nn, aa, strict=True))
    nulls = tuple(sp.factor(r+p) for r, p in zip(rho, pressure, strict=True))
    bianchi = []
    for i in range(3):
        rate = sum(sp.diff(rho[i], a)*v for a, v in zip(aa, velocities, strict=True))/nn[i]
        bianchi.append(sp.factor(rate+3*velocities[i]*nulls[i]/(nn[i]*aa[i])))
    return {"lambda": lam, "a": aa, "n": nn, "beta": beta, "velocities": velocities,
            "S_a": sa, "S_N": sn, "S_a_dot": sadot, "L": lag,
            "rho": rho, "pressure": pressure, "nulls": nulls, "bianchi": tuple(bianchi)}


def evaluate(beta, a, n, *, lambda_=1):
    """Literal polynomial stresses, also defined at singular sums.

    A returned regular_sum=False is not a healthy-branch or bounce verdict.
    All actual metric coframes must remain positive; beta may have any sign.
    """
    if not all(isinstance(values, (list, tuple)) for values in (beta, a, n)):
        raise TypeError("Coefficients and coframes must be finite lists or tuples")
    if not len(beta) == len(a) == len(n) or not beta:
        raise ValueError("At least one metric and equal finite input lengths are required")
    betas = tuple(rational(value) for value in beta)
    scales = tuple(positive(value, "a_i") for value in a)
    lapses = tuple(positive(value, "n_i") for value in n)
    lam = positive(lambda_, "lambda")
    sa, sn = sum(b*x for b, x in zip(betas, scales, strict=True)), sum(
        b*x for b, x in zip(betas, lapses, strict=True))
    rho = tuple(lam*b*sa**3/x**3 for b, x in zip(betas, scales, strict=True))
    pressure = tuple(-lam*b*sn*sa**2/(l*x**2) for b, x, l in zip(betas, scales, lapses, strict=True))
    return {"S_a": sa, "S_N": sn, "L": -lam*sn*sa**3,
            "regular_sum": sa != 0 and sn != 0,
            "rho": rho, "pressure": pressure,
            "nulls": tuple(r+p for r, p in zip(rho, pressure, strict=True)),
            "active": tuple(i for i, b in enumerate(betas) if b != 0)}


def checks():
    d = derive_three()
    lam, sa, sn, sadot = (d[key] for key in ("lambda", "S_a", "S_N", "S_a_dot"))
    out = {}
    for i, (a, n, b, v) in enumerate(zip(d["a"], d["n"], d["beta"], d["velocities"], strict=True)):
        out[f"density_{i}"] = sp.factor(d["rho"][i]-lam*b*sa**3/a**3)
        out[f"pressure_{i}"] = sp.factor(d["pressure"][i]+lam*b*sn*sa**2/(n*a**2))
        out[f"unfactored_Bianchi_{i}"] = sp.factor(
            d["bianchi"][i]-3*lam*b*sa**2/(n*a**3)*(sadot-sn*v/n))
    out["weighted_interaction_null_sum"] = sp.factor(sum(
        n*a**3*null for n, a, null in zip(d["n"], d["a"], d["nulls"], strict=True)))
    return out
