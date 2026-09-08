"""Prepared two-source physical metric variations on the fixed clock."""
from functools import cache

import sympy as sp
from p8_vector_dimensional import local
from p8_vector_hadamard import series
from p8_vector_metric_local import canonical, variation
from p8_vector_metric_local import jets as generic
from p8_vector_state import wkb

u, z = wkb.u, wkb.z
n = sp.symbols("prepared_lapse_source0:12", real=True)
v = sp.symbols("prepared_logscale_source0:12", real=True)
H, lam = wkb.background()["H"], wkb.background()["lambda"]
kind = generic.kind


def output(value):
    if type(value) is not str or value not in ("energy", "pressure"):
        raise ValueError("Require physical energy or pressure readout")
    return value


def clean(value):
    fields = n+v
    polynomial = sp.Poly(sp.expand(value), *fields)
    return sp.Add(*(sp.factor(coefficient)*sp.prod(field**power for field, power in zip(fields, powers))
                    for powers, coefficient in polynomial.terms()))


def time(value):
    return clean(sp.diff(value, u)+wkb.background()["z_prime"]*sp.diff(value, z)
                 +sum(sp.diff(value, row[j])*row[j+1] for row in (n, v) for j in range(len(row)-1)))


def project(value):
    actual = generic.actual(value.subs(generic.D, 3))
    mapping = {generic.z: z}
    mapping.update({field: n[j] for j, field in enumerate(generic.n)})
    mapping.update({field: v[j] for j, field in enumerate(generic.v)})
    return clean(actual.subs(mapping))


def data(sector):
    return _data(kind(sector))


@cache
def _data(sector):
    physical = canonical.data(sector)
    d = H/2 if sector == "T" else H*(sp.Rational(1, 2)+z)
    r, rg = project(physical["delta_log_frequency"]), project(physical["delta_log_g"])
    dd, dl = time(rg), time(r)
    return {"rate": d, "r": r, "rg": rg, "delta_d": dd, "delta_lambda": dl,
            "delta_U": time(dd)+2*d*dd, "c1": d+lam/2, "delta_c1": dd+dl/2}


def readout(sector, component):
    return _readout(kind(sector), output(component))


@cache
def _readout(sector, component):
    physical = canonical.data(sector)
    label, factor = ("N", -sp.Integer(1)) if component == "energy" else ("Z", sp.Rational(1, 3))
    A, B = [factor*project(value) for value in physical["weights"][label]]
    dA, dB = [factor*project(value) for value in physical["varied_weights"][label]]
    normalization = 3*v[0] if component == "energy" else n[0]+3*v[0]
    return {"A": A, "B": B, "delta_A": clean(dA-normalization*A),
            "delta_B": clean(dB-normalization*B)}


def adiabatic(sector, component, order):
    sector, component, order = kind(sector), output(component), generic.order(order)
    return _adiabatic(sector, component, order)


@cache
def _adiabatic(sector, component, order):
    label, factor = ("N", -sp.Integer(1)) if component == "energy" else ("Z", sp.Rational(1, 3))
    varied = factor*project(variation.data(sector)["coefficients"][label][order])
    old = local.energy_coefficients() if component == "energy" else local.pressure_coefficients()
    baseline = project(old["transverse" if sector == "T" else "longitudinal"][order])
    normalization = 3*v[0] if component == "energy" else n[0]+3*v[0]
    return clean(varied-normalization*baseline)


def baseline(sector, order):
    kind(sector)
    if type(order) is not int or not 0 <= order <= 4:
        raise ValueError("Require native reference coefficient order 0..4")
    return series.coefficient("transverse" if sector == "T" else "longitudinal", order)


def linear_bound(value):
    if not isinstance(value, sp.Expr) or value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("Require an exact finite linear two-source expression")
    fields = n+v
    polynomial = sp.Poly(sp.expand(value), *fields)
    total, reconstructions, max_order = sp.Integer(0), [], 0
    for powers, coefficient in polynomial.terms():
        if coefficient == 0:
            continue
        if sum(powers) != 1:
            raise ValueError("Require a linear prepared metric-source coefficient")
        derivative = powers.index(1) % len(n)
        if derivative > 10:
            raise ValueError("The C10 source norm does not control this derivative")
        if coefficient.free_symbols-{u, z}:
            raise ValueError("Only fixed-clock time and momentum fraction may enter source coefficients")
        bound = wkb.box_bound(sp.factor(coefficient))
        total += bound["absolute_upper"]
        reconstructions.append(bound["reconstruction"])
        max_order = max(max_order, derivative)
    return {"absolute_upper": total, "reconstructions": reconstructions, "highest_source_derivative": max_order}
