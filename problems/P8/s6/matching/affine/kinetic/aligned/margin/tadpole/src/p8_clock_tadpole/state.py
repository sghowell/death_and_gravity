"""Explicit, non-circular definition of the selected c-number profiles."""
from functools import cache

import sympy as sp
from p8_vector_clock_matching import bounds as local_matching
from p8_vector_hadamard import preparation as fixed_preparation
from p8_vector_state import comparison, wkb
from p8_vector_subtraction import tail


@cache
def integrands():
    u, z = wkb.u, wkb.z
    a = (1+u**2)**2
    omega = sp.Symbol("positive_mode_frequency", positive=True)
    vr, pr = {}, {}
    for kind in ("transverse", "longitudinal"):
        vr[kind] = sp.Symbol(kind+"_exact_mode_squared", nonnegative=True)
        pr[kind] = sp.Symbol(kind+"_exact_physical_momentum_squared", nonnegative=True)
    result = {}
    for observable in ("energy", "pressure"):
        modes = {}
        for kind in ("transverse", "longitudinal"):
            name = kind+"_"+observable
            weights, ad = tail.physical_weights()[name], tail.reference_terms()[name]
            physical = (weights["A"]*pr[kind]+weights["B"]*omega**2*vr[kind])/(2*a**3)
            subtraction = omega*(ad["zero"]+ad["second"]/omega**2+ad["fourth"]/omega**4)/(4*a**3)
            modes[kind] = {"physical_quadratic_form": physical, "subtraction": subtraction,
                           "physical_momentum_rate": sp.factor(weights["c1"]-wkb.background()["lambda"]/2)}
        finite = sum((2 if kind == "transverse" else 1)*(data["physical_quadratic_form"]-data["subtraction"])
                     for kind, data in modes.items())
        result[observable] = {"modes": modes, "combined_finite_integrand": finite}
    return {"u": u, "z": z, "a": a, "omega": omega, "observables": result}


def definition(planck_time_product, mass_time_product):
    prior = comparison.physical_bounds(planck_time_product, mass_time_product)
    L, m = prior["M_tau"], prior["m0_tau"]
    coefficients = local_matching.actual_coefficients()
    local = {name: sum(data[name]*m**(4-2*n) for n, data in coefficients.items())/(64*sp.pi**2*L**2)
             for name in ("energy", "pressure")}
    return {"M_tau": L, "m0_tau": m,
            "initial_u": -sp.Rational(1, 2),
            "finite_integral_prefactor": 1/L**2,
            "radial_measure": "k_com^2 dk_com/(2pi^2)",
            "local_finite_profile_terms": local,
            "mode_definition": "For each polarization and k_com>=0, nu=sqrt(m²+k_com²/(25/16)²); use the unchanged S6.55 initial_data(kind,m,nu). Set v0=(2W_B)^-1/2 and v0'=(-W_B'/(2W_B)-iW_B)v0. Propagate v''+[omega²-U_kind]v=0 on the original complete clock history, with omega²=m²+k_com²/a² and z=k_com²/(a² omega²).",
            "profile_definition": "Each normalized selected profile is L^-2 times the radial integral of the combined finite integrand in integrands(), plus the displayed matched local term. The two transverse polarizations and the actual lapse-dependent mass readouts are retained.",
            "non_circularity": "Compute these c-number functions once from the predecessor action, state and history. The added scalar action has no vector field, changes neither the clock vector operator nor its Cauchy data, and is not reevaluated under metric or state variations.",
            "global_definition_boundary": "The exact oscillator solution and every finite-order high-frequency expansion extend over every compact real-time interval. Their arbitrary-order remainder gives a C-infinity finite subtracted stress on the full original history. Quantitative bounds are supplied only on |u|<=1/2; no uniform all-time smallness claim is made."}


def initial_data_examples():
    return {kind: fixed_preparation.initial_data(kind, 1000, 3000)
            for kind in ("transverse", "longitudinal")}


def checks():
    data = integrands()
    out = {}
    H = wkb.background()["H"]
    for name, obs in data["observables"].items():
        out[name+"_transverse_physical_canonical_rate"] = sp.factor(obs["modes"]["transverse"]["physical_momentum_rate"]-H/2)
        out[name+"_longitudinal_physical_canonical_rate"] = sp.factor(obs["modes"]["longitudinal"]["physical_momentum_rate"]-H*(sp.Rational(1, 2)+wkb.z))
    return out
