"""Both prepared acoustic clocks and their fixed-time covariance sources."""
from functools import cache

import sympy as sp
from p8_vector_liouville import clock, reduction, variation
from p8_vector_metric_response import source as physical

from . import transport

u, a, mass2 = clock.u, clock.a, clock.mass2
xi_T, xi_L = sp.symbols("prepared_transverse_time_shift prepared_longitudinal_time_shift", real=True)
fields = variation.fields+(xi_T, xi_L)


def clean(value):
    if not isinstance(value, (int, sp.Expr)) or isinstance(value, (bool, float)):
        raise TypeError("Require an exact source expression")
    value = sp.sympify(value)
    if value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("Require a finite exact acoustic source expression")
    poly = sp.Poly(sp.expand(value), *fields)
    if poly.is_zero:
        return sp.Integer(0)
    if any(sum(powers) != 1 for powers, _ in poly.terms()):
        raise ValueError("Require a homogeneous linear source or prepared time shift")
    if any(coefficient.free_symbols-{u, mass2, reduction.k} for _, coefficient in poly.terms()):
        raise ValueError("Require fixed-clock time and momentum coefficients")
    return sp.Add(*(sp.factor(coefficient)*sp.prod(field**power for field, power in zip(fields, powers))
                    for powers, coefficient in poly.terms()))


def project(value):
    mapping = {physical.z: reduction.k**2/(reduction.k**2+clock.data()["U"])}
    mapping.update({physical.n[j]: variation.n[j] for j in range(4)})
    mapping.update({physical.v[j]: variation.zeta[j] for j in range(4)})
    polynomial = sp.Poly(sp.expand(value.subs(mapping, simultaneous=True)), *fields)
    return sp.Add(*(sp.factor(coefficient)*sp.prod(field**power for field, power in zip(fields, powers))
                    for powers, coefficient in polynomial.terms()))


def data(sector):
    return _data(reduction.sector(sector))


@cache
def _data(sector):
    bg, src, varied = clock.data(), variation.source(), variation.data()
    D = reduction.k**2+bg["U"]
    if sector == "T":
        rate = variation.n[0]-variation.zeta[0]
        dU = bg["U"]*(2*variation.zeta[0]+src["beta"]*variation.n[0])
        bk, dbk, xi = sp.Integer(0), sp.Integer(0), xi_T
    else:
        rate = src["acoustic_log_rate"]
        dU = varied["delta_U"]
        bk = bg["b"]-bg["U1"]/(2*D)
        dbk = varied["delta_pump_rate"]-varied["delta_U1"]/(2*D)+bg["U1"]*dU/(2*D**2)
        xi = xi_L
    bk = sp.factor(bk)
    dU, dbk = clean(dU), clean(dbk)
    dUs = clean(dU-xi*clock.derivative(bg["U"]))
    dbks = clean(dbk-xi*clock.derivative(bk))
    return {"rate": 1/a, "log_rate_variation": clean(rate), "prepared_shift": xi,
            "D": D, "bk": bk, "delta_U_fixed_u": dU, "delta_bk_fixed_u": dbk,
            "delta_U_fixed_sigma": dUs, "delta_bk_fixed_sigma": dbks,
            "matrix": transport.matrix(bk, D),
            "delta_matrix_fixed_sigma": transport.varied_generator().subs({
                transport.dbk: dbks, transport.dU: dUs}, simultaneous=True)}


def derivative(value, sector):
    item = data(sector)
    return clean(variation.acoustic(value)+item["log_rate_variation"]*sp.diff(value, item["prepared_shift"]))


@cache
def checks():
    out = {}
    H = 4*u/(1+u**2)
    for sector in ("T", "L"):
        item = data(sector)
        old = physical.data(sector)
        rate = item["log_rate_variation"]
        old_rate = project(old["rate"])
        actual_dbk = a*(project(old["delta_d"])+variation.time(rate)/2-rate*(old_rate-H/2))
        out[sector+"_physical_frequency_source"] = clean(
            project(old["r"])-rate-item["delta_U_fixed_u"]/(2*item["D"]))
        out[sector+"_physical_moving_canonical_source"] = clean(actual_dbk-item["delta_bk_fixed_u"])
        out[sector+"_background_covariance_half_rate"] = sp.factor(a*(old_rate-H/2)-item["bk"])
        delta_potential = clean(item["delta_U_fixed_sigma"]-derivative(item["delta_bk_fixed_sigma"], sector)
                                -2*item["bk"]*item["delta_bk_fixed_sigma"])
        if sector == "T":
            expected = item["delta_U_fixed_sigma"]
        else:
            V = variation.data()
            expected = (V["delta_scalar_principal_fixed_acoustic_time"]
                        +V["delta_remainder_fixed_acoustic_time"]).subs(variation.xi, xi_L)
        out[sector+"_potential_response_reproduces_acoustic_covariance_source"] = clean(delta_potential-expected)
    out["two_clock_shifts_are_not_identified"] = clean(
        data("L")["log_rate_variation"]-data("T")["log_rate_variation"]
        -(variation.source()["beta"]-variation.source()["alpha"])*variation.n[0]/2)
    return out
