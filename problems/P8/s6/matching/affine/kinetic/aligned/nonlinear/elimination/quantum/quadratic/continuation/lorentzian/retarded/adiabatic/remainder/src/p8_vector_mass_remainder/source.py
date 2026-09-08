"""Actual-clock C10 source jets and exact linear coefficient bounds."""
from functools import cache

import sympy as sp
from p8_vector_hadamard import series
from p8_vector_mass_adiabatic import jets as generic
from p8_vector_state import wkb

u, z = wkb.u, wkb.z
n = sp.symbols("prepared_source0:12", real=True)
H = wkb.background()["H"]
lam = wkb.background()["lambda"]
alpha, beta = 4/(9*(1+u**2)**3), 28/(81*(1+u**2)**3)


def kind(value):
    if type(value) is not str or value not in ("T", "L"):
        raise ValueError("Require physical T or L mode")
    return value


def clean(value):
    polynomial = sp.Poly(sp.expand(value), *n)
    return sp.Add(*(sp.factor(coefficient)*sp.prod(field**power for field, power in zip(n, powers))
                    for powers, coefficient in polynomial.terms()))


def time(value):
    return clean(sp.diff(value, u)+wkb.background()["z_prime"]*sp.diff(value, z)
                 +sum(sp.diff(value, n[j])*n[j+1] for j in range(len(n)-1)))


def project(value):
    mapping = {generic.z: z, generic.D: 3}
    for row, function in ((generic.H, H), (generic.alpha, alpha), (generic.beta, beta)):
        mapping.update({field: sp.diff(function, u, j) for j, field in enumerate(row)})
    mapping.update({field: n[j] for j, field in enumerate(generic.n)})
    return clean(value.subs(mapping))


def data(sector):
    return _data(kind(sector))


@cache
def _data(sector):
    d = H/2 if sector == "T" else H*(sp.Rational(1, 2)+z)
    r = beta*(1-z)*n[0]/2 if sector == "T" else (beta-alpha*z)*n[0]/2
    rg = sp.Integer(0) if sector == "T" else alpha*z*n[0]/2
    dd, dl = time(rg), time(r)
    if sector == "T":
        A, B, dA = sp.Integer(0), -beta*(1-z), sp.Integer(0)
        dB = -2*r*B
    else:
        A, B = alpha*z, -beta
        dA, dB = -alpha**2*z*(2-z)*n[0], beta**2*n[0]
    return {"rate": d, "r": r, "rg": rg, "delta_d": dd, "delta_lambda": dl,
            "delta_U": time(dd)+2*d*dd, "A": A, "B": B,
            "delta_A": dA, "delta_B": dB, "c1": d+lam/2, "delta_c1": dd+dl/2}


def baseline(sector, order):
    kind(sector)
    if type(order) is not int or not 0 <= order <= 4:
        raise ValueError("Require native reference coefficient order 0..4")
    return series.coefficient("transverse" if sector == "T" else "longitudinal", order)


def linear_bound(value):
    if not isinstance(value, sp.Expr) or value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("Require an exact finite linear source expression")
    polynomial = sp.Poly(sp.expand(value), *n)
    total, reconstructions, max_order = sp.Integer(0), [], 0
    for powers, coefficient in polynomial.terms():
        if coefficient == 0:
            continue
        if sum(powers) != 1:
            raise ValueError("Require a linear prepared-source coefficient")
        derivative = powers.index(1)
        if derivative > 10:
            raise ValueError("The C10 source norm does not control this derivative")
        if coefficient.free_symbols-{u, z}:
            raise ValueError("Only fixed-clock time and momentum fraction may enter source coefficients")
        bound = wkb.box_bound(sp.factor(coefficient))
        total += bound["absolute_upper"]
        reconstructions.append(bound["reconstruction"])
        max_order = max(max_order, derivative)
    return {"absolute_upper": total, "reconstructions": reconstructions, "highest_source_derivative": max_order}
