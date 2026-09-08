"""Retained vector-block principal speed, not a full quantum cone verdict."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_aligned_quantum import kernel
from p8_vector_variation import modes


@cache
def principal():
    N, am, bm, q, m2 = modes.N, modes.am, modes.bm, modes.q, modes.mass2
    H, Hd, Nd, Ndd, amd, amdd = sp.symbols("H H_dot N_dot N_ddot am_dot am_ddot", real=True)
    data = modes.canonical()
    dL = (H+amd/am-Nd/N+(2*H*q-m2*amd)/(q+m2*am))/2
    dT = (H-Nd/N)/2
    flows = {H: Hd, N: Nd, Nd: Ndd, am: amd, amd: amdd, q: -2*H*q}
    derivative = lambda value: sum(sp.diff(value, var)*flow for var, flow in flows.items())
    out, checks = {}, {}
    for kind, rate in (("transverse", dT), ("longitudinal", dL)):
        correction = sp.factor(derivative(rate)+rate**2)
        speed = sp.limit((data[kind]["bare_frequency_squared"]-correction)/(N**2*q), q, sp.oo)
        target = 1 if kind == "transverse" else bm/am
        out[kind] = {"canonical_logarithmic_rate": rate, "squared_principal_speed": speed}
        checks[kind+"_bounded_rate_correction_does_not_change_principal_speed"] = sp.limit(correction/q, q, sp.oo)
        checks[kind+"_physical_metric_principal_normalization"] = sp.factor(speed-target)
    return {"modes": out, "checks": checks}


@cache
def actual():
    p = kernel.geometry.P
    mass = kernel.masses()
    speed = sp.factor(mass["b"]/mass["a"])
    polynomial = 24*p**3+12*p**2+6*p-5
    factor = -3*p*polynomial/((8*p+5)*(22*p**3+3*p-11))
    return {"p": p, "longitudinal_squared_principal_speed": speed,
            "positive_sign_factor": factor, "sign_polynomial": polynomial,
            "p_derivative_at_clock": sp.factor(sp.diff(speed, p).subs(p, sp.Rational(1, 2))),
            "checks": {"actual_longitudinal_speed_sign_factorization": sp.factor(speed-1-(2*p-1)*factor),
                       "actual_clock_front_is_luminal": sp.factor(speed.subs(p, sp.Rational(1, 2))-1),
                       "actual_clock_slope": sp.factor(sp.diff(speed, p).subs(p, sp.Rational(1, 2))-sp.Rational(16, 81))}}


def speed(p_affine):
    value = exact(p_affine, "p_affine")
    if not bool(sp.Rational(9, 20) <= value <= sp.Rational(21, 40)):
        raise ValueError("Require exact 9/20<=p_affine<=21/40")
    data = actual()
    c2 = sp.factor(data["longitudinal_squared_principal_speed"].subs(data["p"], value))
    return {"p_affine": value, "longitudinal_squared_principal_speed": c2,
            "retained_longitudinal_block_subluminal_or_luminal": bool(c2 <= 1),
            "full_coupled_quantum_or_UV_front_claim": False}


def proof_checks():
    lo, hi = sp.Rational(9, 20), sp.Rational(21, 40)
    return {"positive_sign_polynomial_on_full_probe_interval": bool(24*lo**3+12*lo**2+6*lo-5 > 0),
            "temporal_mass_numerator_negative": bool(22*hi**3+3*hi-11 < 0),
            "temporal_mass_denominator_negative": bool(2*hi**3-1 < 0),
            "spatial_mass_numerator_negative": bool(72*hi**2-88*lo-55 < 0),
            "original_clock_tube_contained_in_probe_interval": bool(lo**2 < sp.Rational(9, 40) and hi**2 > sp.Rational(11, 40)),
            "strictly_positive_actual_off_clock_front_slope": bool(sp.Rational(16, 81) > 0)}
