"""Frozen actual-clock replay and independent Gamma continuation controls."""
from functools import cache

import sympy as sp
from p8_vector_state import wkb
from p8_vector_subtraction import tail

from . import local


@cache
def physical_checks():
    H, u = wkb.background()["H"], wkb.u
    mapping = {local.dimension: 3, local.z: wkb.z,
               local.alpha: 4/(9*(1+u**2)**3), local.beta: 28/(81*(1+u**2)**3)}
    mapping.update({symbol: sp.diff(H, u, order) for order, symbol in
                    enumerate((local.H, local.Hd, local.Hdd, local.Hddd, local.Hdddd))})
    out = {}
    for kind in ("transverse", "longitudinal"):
        for key in ("U", "P2", "P4", "B2"):
            out[kind+"_actual_clock_"+key] = sp.factor(local.reference(kind)[key].subs(mapping)-wkb.frequency(kind)[key])
        for observable, generator in (("energy", local.energy_coefficients), ("pressure", local.pressure_coefficients)):
            actual = generator()[kind]
            prior = tail.reference_terms()[kind+"_"+observable]
            for order, key in enumerate(("zero", "second", "fourth")):
                out[kind+"_actual_"+observable+"_"+key] = sp.factor(actual[order].subs(mapping)-prior[key])
    return out


@cache
def gamma_checks():
    e = sp.Symbol("epsilon_DR", positive=True)
    out = {}
    for k in range(3):
        residue = sp.Rational((-1)**k, sp.factorial(k))
        actual = sp.gamma(e-k)
        out["Gamma_pole_"+str(k)] = sp.limit(e*actual, e, 0)-residue
        out["Gamma_finite_"+str(k)] = sp.simplify(sp.limit(actual-residue/e, e, 0)
                                               -residue*(sp.harmonic(k)-sp.EulerGamma))
    for j in range(6):
        out["Gamma_angular_ratio_"+str(j)] = sp.combsimp(sp.gamma(j+local.dimension/2)/sp.gamma(local.dimension/2)
                                                       -sp.rf(local.dimension/2, j))
    return out
