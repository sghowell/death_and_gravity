"""Prepared physical variation including acoustic-time transport."""
from functools import cache

import sympy as sp
from p8_vector_mass_adiabatic.bounds import rational_bound
from p8_vector_metric_dispersion import principal

from . import clock, reduction

u, a, mass2 = clock.u, clock.a, clock.mass2
n = sp.symbols("physical_lapse_jet0:4", real=True)
zeta = sp.symbols("physical_scale_jet0:4", real=True)
fields = n+zeta
xi = sp.Symbol("initial_zero_acoustic_time_shift", real=True)


def time(value):
    return sp.expand(sp.diff(value, u)+sum(
        sp.diff(value, row[j])*row[j+1] for row in (n, zeta) for j in range(3)))


def acoustic(value):
    return sp.expand(a*time(value))


def clean(value):
    if not isinstance(value, (int, sp.Expr)) or isinstance(value, (bool, float)):
        raise TypeError("Require an exact source-linear expression")
    value = sp.sympify(value)
    if value.has(sp.Float) or value.is_finite is False:
        raise ValueError("Require exact finite coefficients")
    poly = sp.Poly(sp.expand(value), *fields)
    if poly.is_zero:
        return sp.Integer(0)
    if any(sum(powers) != 1 for powers, _ in poly.terms()):
        raise ValueError("Require a homogeneous source-linear expression")
    return sp.Add(*(sp.factor(coefficient)*sp.prod(field**power for field, power in zip(fields, powers))
                    for powers, coefficient in poly.terms()))


def coefficient_bound(value):
    value = clean(value)
    if value.has(n[3], zeta[3]):
        raise ValueError("Third source derivatives lie outside the C2 domain")
    coefficients = [sp.factor(sp.expand(value).coeff(field)) for field in fields[:3]+fields[4:7]]
    if any(coefficient.free_symbols-{u} for coefficient in coefficients):
        raise ValueError("Require fixed rational background coefficients")
    return sum(rational_bound(coefficient) for coefficient in coefficients)


@cache
def source():
    h = (1+u**2)**3
    alpha, beta = 4/(9*h), 28/(81*h)
    c, rate_n = (alpha+beta)/4, 1+(beta-alpha)/2
    return {"h": h, "alpha": alpha, "beta": beta,
            "pump": zeta[0]+c*n[0],
            "acoustic_log_rate": rate_n*n[0]-zeta[0],
            "mass_log_variation": 2*zeta[0]+alpha*n[0],
            "physical_to_acoustic_sources": sp.ImmutableMatrix([[c, 1], [rate_n, -1]])}


@cache
def data():
    bg, src = clock.data(), source()
    chi, rate, j = src["pump"], src["acoustic_log_rate"], src["mass_log_variation"]
    U, U1, U2, b, curv = (bg[key] for key in ("U", "U1", "U2", "b", "pump_curvature"))
    dU = U*j
    dU1 = clean(acoustic(dU)-rate*U1)
    dU2 = clean(acoustic(dU1)-rate*U2)
    db = clean(acoustic(chi)-rate*b)
    dcurv = clean(acoustic(acoustic(chi))+2*b*acoustic(chi)-b*acoustic(rate)-2*rate*curv)
    numerator = b*U1+U2/2
    dnumerator = clean(db*U1+b*dU1+dU2/2)
    varied_rows = (dnumerator, clean(-numerator*dU-sp.Rational(3, 2)*U1*dU1),
                   clean(sp.Rational(3, 2)*U1**2*dU))
    history_rows = (clock.derivative(numerator),
                    -numerator*U1-sp.Rational(3, 2)*U1*U2, sp.Rational(3, 2)*U1**3)
    D = reduction.k**2+U
    fixed_time = sum(row/D**power for power, row in enumerate(varied_rows, 1))
    history_rate = sum(row/D**power for power, row in enumerate(history_rows, 1))
    principal_fixed_time = clean(dU-dcurv)
    return {"delta_U": clean(dU), "delta_U1": dU1, "delta_U2": dU2,
            "delta_pump_rate": db, "delta_pump_curvature": dcurv,
            "varied_remainder_numerators": varied_rows,
            "history_remainder_numerators": tuple(sp.factor(row) for row in history_rows),
            "delta_remainder_fixed_physical_time": fixed_time,
            "delta_remainder_fixed_acoustic_time": fixed_time-xi*history_rate,
            "delta_scalar_principal_fixed_physical_time": principal_fixed_time,
            "delta_scalar_principal_fixed_acoustic_time": principal_fixed_time-xi*clock.derivative(bg["scalar_principal_potential"])}


@cache
def estimates():
    item = data()
    varied = tuple(coefficient_bound(row/mass2**power) for power, row in enumerate(
        item["varied_remainder_numerators"], 1))
    history = tuple(rational_bound(sp.factor(row/mass2**power)) for power, row in enumerate(
        item["history_remainder_numerators"], 1))
    total = sum(varied)+2*sum(history)
    return {"varied_numerator_C2_bounds": varied, "history_numerator_bounds": history,
            "initial_zero_time_shift_C0_upper_per_source_C0": sp.Integer(2),
            "correction_C2_coefficient": total, "coarse_correction_C2_coefficient": sp.Integer(18000),
            "scope": "The first variation of the mode-potential remainder at fixed acoustic time, not the renormalized stress or a coupled inverse."}


def bound(momentum, mass):
    momentum, mass = clock.exact_rational(momentum, True), clock.exact_rational(mass)
    fraction = mass**2/(momentum**2+mass**2)
    return {"momentum": momentum, "mass": mass,
            "prepared_remainder_C2_to_C0_upper": estimates()["correction_C2_coefficient"]*fraction,
            "coarse_prepared_remainder_C2_to_C0_upper": 18000*fraction}


@cache
def checks():
    src, item, bg = source(), data(), clock.data()
    e = sp.Symbol("variation_parameter", real=True)
    rate = src["acoustic_log_rate"]
    chi = src["pump"]
    # Directly vary the physical-to-acoustic derivative, not just the pump.
    def perturbed_derivative(value):
        return a*sp.exp(-e*rate)*time(value)
    Ae = a*sp.exp(e*chi)
    Ue = bg["U"]*sp.exp(e*src["mass_log_variation"])
    be = perturbed_derivative(Ae)/Ae
    ce = perturbed_derivative(perturbed_derivative(Ae))/Ae
    out = {
        "varied_acoustic_pump_rate": clean(sp.diff(be, e).subs(e, 0)-item["delta_pump_rate"]),
        "varied_acoustic_pump_curvature": clean(sp.diff(ce, e).subs(e, 0)-item["delta_pump_curvature"]),
        "varied_acoustic_mass_first": clean(sp.diff(perturbed_derivative(Ue), e).subs(e, 0)-item["delta_U1"]),
        "varied_acoustic_mass_second": clean(sp.diff(perturbed_derivative(perturbed_derivative(Ue)), e).subs(e, 0)-item["delta_U2"]),
    }
    d0, d1, d2, db = sp.symbols("delta_U delta_U1 delta_U2 delta_b", real=True)
    original = reduction.longitudinal_potential()["momentum_remainder"]
    varied = sp.diff(original.subs({reduction.U: reduction.U+e*d0, reduction.U1: reduction.U1+e*d1,
                                   reduction.U2: reduction.U2+e*d2, reduction.b: reduction.b+e*db}, simultaneous=True), e).subs(e, 0)
    replacement = {reduction.U: bg["U"], reduction.U1: bg["U1"], reduction.U2: bg["U2"], reduction.b: bg["b"],
                   d0: item["delta_U"], d1: item["delta_U1"], d2: item["delta_U2"], db: item["delta_pump_rate"]}
    out["independent_rational_remainder_variation"] = clean(
        varied.subs(replacement, simultaneous=True)-item["delta_remainder_fixed_physical_time"])
    replacement = {reduction.U: bg["U"], reduction.U1: bg["U1"], reduction.U2: bg["U2"], reduction.b: bg["b"]}
    remainder = original.subs(replacement)
    target_history = sum(row/(reduction.k**2+bg["U"])**power
                         for power, row in enumerate(item["history_remainder_numerators"], 1))
    out["independent_fixed_acoustic_time_history"] = sp.factor(clock.derivative(remainder)-target_history)
    T = src["physical_to_acoustic_sources"]
    out["regular_acoustic_source_chart_determinant"] = sp.factor(
        T.det()+(1+4/(27*src["h"])))
    out["highest_derivative_pump_reproduces_rank_one_pole"] = sp.ImmutableMatrix(
        2*sp.ImmutableMatrix([sp.diff(chi, n[0]), sp.diff(chi, zeta[0])])
        *sp.ImmutableMatrix([[sp.diff(chi, n[0]), sp.diff(chi, zeta[0])]])
        -principal.pole().subs(principal.h, src["h"]))
    for field, coefficient in ((n[2], -(src["alpha"]+src["beta"])*a**2/4), (zeta[2], -a**2)):
        out[str(field)+"_principal_second_derivative"] = sp.factor(
            sp.expand(item["delta_scalar_principal_fixed_physical_time"]).coeff(field)-coefficient)
    out["omitted_acoustic_rate_bounce_control"] = sp.factor(
        (bg["b"]*acoustic(rate)+2*rate*bg["pump_curvature"]).subs(u, 0)
        -sp.Rational(616, 81)*n[0]+8*zeta[0])
    return out
