"""Ordinary mode conservation and the actual homogeneous clock source."""
from functools import cache

import sympy as sp
from p8_vector_clock_matching import bounds as matched
from p8_vector_state import wkb
from p8_vector_subtraction import tail

from . import majorants


@cache
def ordinary_adiabatic_conservation():
    generic, weights, H = tail.algebra(), tail.physical_weights(), wkb.background()["H"]
    out = {}
    for kind in ("transverse", "longitudinal"):
        data, pressure = majorants.coefficients(kind), weights[kind+"_pressure"]
        mapping = {tail.p: data["P2"], tail.r: data["P4"], tail.b2: data["B2"],
                   tail.b4: data["B4"], tail.c: pressure["c1"]}
        for n, name in enumerate(("zero", "second", "fourth")):
            energy = sp.factor(generic[name].subs({**mapping, tail.wa: 1, tail.wb: 1}))
            p = sp.factor(generic[name].subs({**mapping, tail.wa: pressure["A"], tail.wb: pressure["B"]}))
            out[kind+"_ordinary_full_adiabatic_conservation_"+str(2*n)] = sp.factor(
                wkb.D(energy)+(1-2*n)*wkb.background()["lambda"]*energy+3*H*p)
    return out


@cache
def clock_variation():
    u, N, velocity = sp.symbols("u N clock_velocity", real=True)
    alpha, beta, La, Lb, H = [sp.Function(name)(u) for name in ("alpha", "beta", "La", "Lb", "H")]
    qa, qb = sp.Function("qa")(u), sp.Function("qb")(u)
    x = -velocity**2/N**2
    am, bm = 1+alpha*(x+1)/2+qa*(x+1)**2, 1+beta*(x+1)/2+qb*(x+1)**3
    clock = {N: 1, velocity: 1}
    # La,Lb stand for arbitrary on-clock mass-insertion readouts.
    density = N*(La*(am-1)+Lb*(bm-1))
    mass_energy = -sp.diff(density, N).subs(clock)
    momentum = sp.diff(density, velocity).subs(clock)
    # The explicit clock-coordinate derivative of each mass function
    # vanishes on the clock. Derivatives of La,Lb belong to integration
    # by parts, not the partial mass-function derivative.
    explicit_clock = (La*sp.diff(am, u)+Lb*sp.diff(bm, u)).subs(clock)
    euler = explicit_clock-sp.diff(momentum, u)-3*H*momentum
    return {"homogeneous_clock_mass_jet_a": sp.diff(am, velocity).subs(clock)+alpha,
            "homogeneous_clock_mass_jet_b": sp.diff(bm, velocity).subs(clock)+beta,
            "homogeneous_clock_mass_explicit_partial_vanishes": explicit_clock,
            "clock_canonical_momentum_equals_mass_lapse_energy": sp.simplify(momentum-mass_energy),
            "clock_source_equals_negative_mass_energy_balance": sp.simplify(euler+sp.diff(mass_energy, u)+3*H*mass_energy)}


@cache
def local_coefficients():
    H, u = wkb.background()["H"], wkb.u
    out = {}
    for n, d in matched.actual_coefficients().items():
        out[n] = {"energy_derivative": sp.factor(sp.diff(d["energy"], u)),
                  "pressure_derivative": sp.factor(sp.diff(d["pressure"], u)),
                  "negative_clock_source": sp.factor(sp.diff(d["energy"], u)+3*H*(d["energy"]+d["pressure"]))}
    return out


@cache
def local_envelopes():
    return {n: {name: wkb.box_bound(value)["absolute_upper"] for name, value in d.items()}
            for n, d in local_coefficients().items()}


def checks():
    out = {**ordinary_adiabatic_conservation(), **clock_variation(),
           **ordinary_exact_mode_balance(), **local_mass_balance()}
    for n, d in local_coefficients().items():
        out.update({name+"_box_reconstruction_"+str(n): wkb.box_bound(value)["reconstruction"] for name, value in d.items()})
    return out


@cache
def ordinary_exact_mode_balance():
    D, z, H = sp.symbols("D z H")
    lam = -H*z
    out = {}
    for kind, rate, pa, pb in (
            ("transverse", H*(D-2)/2, (D-2)/D, (2*z-D+2)/D),
            ("longitudinal", H*((D-2)/2+z), (D-2+2*z)/D, -(D-2)/D)):
        out[kind+"_exact_ordinary_dimensional_kinetic_balance"] = sp.factor(-2*rate+D*H*pa)
        out[kind+"_exact_ordinary_dimensional_potential_balance"] = sp.factor(2*lam+2*rate+D*H*pb)
    return out


@cache
def local_mass_balance():
    H, u = wkb.background()["H"], wkb.u
    ordinary_energy = {0: -sp.Rational(5, 2), 1: -10*H**2,
                       2: 2*(6*H**2*sp.diff(H, u)+2*H*sp.diff(H, u, 2)-sp.diff(H, u)**2)}
    out = {}
    for n, actual in matched.actual_coefficients().items():
        mass_energy = actual["energy"]-ordinary_energy[n]
        out["ordinary_local_conservation_"+str(n)] = sp.factor(
            sp.diff(ordinary_energy[n], u)+3*H*(ordinary_energy[n]+actual["pressure"]))
        out["local_clock_source_is_mass_energy_balance_"+str(n)] = sp.factor(
            local_coefficients()[n]["negative_clock_source"]-sp.diff(mass_energy, u)-3*H*mass_energy)
    return out
