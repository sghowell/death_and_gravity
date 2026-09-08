"""Actual rolling local-coefficient envelopes in the candidate prescription."""
from functools import cache

import sympy as sp
from p8_vector_dimensional import local
from p8_vector_state import wkb
from p8_vector_subtraction import tail

from . import counterterms


@cache
def actual_coefficients():
    H, u = wkb.background()["H"], wkb.u
    mapping = {local.alpha: 4/(9*(1+u**2)**3), local.beta: 28/(81*(1+u**2)**3), local.ell: 0}
    mapping.update({symbol: sp.diff(H, u, order) for order, symbol in
                    enumerate((local.H, local.Hd, local.Hdd, local.Hddd, local.Hdddd))})
    return {n: {key: sp.factor(d[key].subs(mapping)) for key in ("energy", "pressure")}
            for n, d in counterterms.finite_coefficients().items()}


@cache
def continuous_envelopes():
    coefficients = actual_coefficients()
    out = {"energy": {0: sp.Integer(4)}, "pressure": {0: sp.Rational(5, 2)}}
    for n in (1, 2):
        for key, envelope in out.items():
            envelope[n] = wkb.box_bound(coefficients[n][key])["absolute_upper"]
    return out


def physical_bounds(planck_time_product, reference_mass_time_product):
    prior = tail.physical_bounds(planck_time_product, reference_mass_time_product)
    L, R = prior["M_tau"], prior["m0_tau"]
    envelopes = continuous_envelopes()
    out = {"M_tau": L, "m0_tau": R}
    for key in ("energy", "pressure"):
        b = envelopes[key]
        local_bound = (b[0]*R**4+b[1]*R**2+b[2])/(576*L**2)
        out[key+"_finite_local_over_reference_density"] = local_bound
        out[key+"_candidate_matched_total_over_reference_density"] = prior["full_subtracted_"+key+"_integral_over_reference_density"]+local_bound
    out["unique_scheme_or_full_quantum_bounce_claim"] = False
    return out


@cache
def checks():
    coefficients = actual_coefficients()
    h = (1+wkb.u**2)**3
    out = {"actual_frozen_finite_flat_energy": sp.factor(coefficients[0]["energy"]+sp.Rational(5, 2)+22/(27*h)),
           "actual_frozen_finite_flat_pressure": coefficients[0]["pressure"]-sp.Rational(5, 2)}
    for n in (1, 2):
        for key in ("energy", "pressure"):
            out[key+"_local_coefficient_box_reconstruction_"+str(n)] = wkb.box_bound(coefficients[n][key])["reconstruction"]
    return out
