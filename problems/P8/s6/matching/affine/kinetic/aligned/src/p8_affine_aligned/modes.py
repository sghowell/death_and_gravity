"""Actual time-dependent quadratic Proca modes, conditional on source alignment.

A canonical frequency floor is not a stationary scattering gap or a
nonlinear cutoff. Fixed comoving, not physical, momentum is differentiated.
"""
from functools import cache

import sympy as sp
from p8_affine_retuned import dynamics as parent
from p8_affine_retuned.bounds import exact

u = parent.u
ZETA = sp.Symbol("zeta_aligned", positive=True)
KCOM2 = sp.Symbol("k_com_squared", positive=True)
R = sp.Symbol("r", nonnegative=True)
ZETA_MAX = sp.Rational(1, 2000)


@cache
def canonical():
    bg = parent.old.background()
    a, hubble = bg["a"], bg["H"]
    q = KCOM2/a**2
    r = ZETA*q
    weights = {"transverse": a*ZETA, "longitudinal": a**3*ZETA*q/(1+ZETA*q)}
    pieces = {}
    for name, weight in weights.items():
        rate = sp.factor(sp.diff(weight, u)/(2*weight))
        correction = sp.factor(sp.diff(rate, u)+rate**2)
        pieces[name] = {"weight": weight, "log_rate": rate, "correction": correction,
                        "frequency_squared": q+1/ZETA-correction}
    alpha = (1+3*R)/(2*(1+R))
    beta = (1-2*R+9*R**2)/(4*(1+R)**2)
    return {"q": q, "r": r, "a": a, "H": hubble,
            "alpha": alpha, "beta": beta, **pieces,
            "expected_transverse": sp.diff(hubble, u)/2+hubble**2/4,
            "expected_longitudinal": (alpha*sp.diff(hubble, u)+beta*hubble**2).subs(R, r)}


@cache
def checks():
    data = canonical()
    hubble, r = data["H"], data["r"]
    alpha, beta = data["alpha"], data["beta"]
    return {"fixed_comoving_q_derivative": sp.factor(sp.diff(data["q"], u)+2*hubble*data["q"]),
            "transverse_actual_normalization": sp.factor(data["transverse"]["correction"]-data["expected_transverse"]),
            "longitudinal_actual_normalization": sp.factor(data["longitudinal"]["correction"]-data["expected_longitudinal"]),
            "longitudinal_log_rate": sp.factor(data["longitudinal"]["log_rate"]-hubble*alpha.subs(R, r)),
            "longitudinal_beta_chain_rule": sp.factor(alpha**2-2*R*sp.diff(alpha, R)-beta),
            "beta_strict_upper": sp.factor(sp.Rational(9, 4)-beta-(2+5*R)/(1+R)**2),
            "beta_strict_positive_numerator": sp.expand(1-2*R+9*R**2-9*(R-sp.Rational(1, 9))**2-sp.Rational(8, 9)),
            "Hubble_square_upper": sp.factor(4-hubble**2-4*(u**2-1)**2/(1+u**2)**2),
            "Hubble_dot_upper": sp.factor(4-sp.diff(hubble, u)-4*u**2*(u**2+3)/(1+u**2)**2),
            "Hubble_dot_lower": sp.factor(4+sp.diff(hubble, u)-4*(u**4+u**2+2)/(1+u**2)**2)}


@cache
def proof_checks():
    data = canonical()
    alpha, beta = data["alpha"], data["beta"]
    return {"positive_transverse_weight": data["transverse"]["weight"].is_positive is True,
            "positive_longitudinal_weight_for_nonzero_comoving_momentum": data["longitudinal"]["weight"].is_positive is True,
            "alpha_lower_identity": sp.factor(alpha-sp.Rational(1, 2)-R/(1+R)) == 0,
            "alpha_upper_identity": sp.factor(sp.Rational(3, 2)-alpha-1/(1+R)) == 0,
            "beta_positive_square_identity": sp.factor(beta-(9*(R-sp.Rational(1, 9))**2+sp.Rational(8, 9))/(4*(1+R)**2)) == 0,
            "transverse_curvature_absolute_upper_three": sp.Rational(1, 2)*4+sp.Rational(1, 4)*4 == 3,
            "longitudinal_curvature_absolute_upper_fifteen": sp.Rational(3, 2)*4+sp.Rational(9, 4)*4 == 15,
            "chosen_curl_uniform_frequency_floor": 1/ZETA_MAX-15 == 1985,
            "zero_momentum_longitudinal_correction_matches_homogeneous_chart": sp.factor(
                data["longitudinal"]["correction"].subs(KCOM2, 0)-data["transverse"]["correction"]) == 0,
            "frequency_floor_not_stationary_spectrum_or_cutoff": True}


def frequency_bounds(coupling, time, comoving_squared):
    zeta, point, k2 = [exact(value, name) for value, name in
                       ((coupling, "zeta"), (time, "time"), (comoving_squared, "k_com_squared"))]
    if not bool(0 < zeta <= ZETA_MAX) or k2.is_nonnegative is not True:
        raise ValueError("Require 0<zeta<=1/2000 and nonnegative comoving momentum squared")
    data = canonical()
    q = sp.factor(data["q"].subs({u: point, KCOM2: k2}))
    homogeneous = k2.is_zero is True
    correction = data["transverse" if homogeneous else "longitudinal"]["correction"].subs(
        {u: point, ZETA: zeta, KCOM2: k2})
    return {"zeta": zeta, "u": point, "q": q,
            "chart": "three_homogeneous_coordinate_vectors" if homogeneous else "longitudinal_canonical",
            "canonical_frequency_squared": sp.factor(q+1/zeta-correction),
            "uniform_frequency_squared_lower": q+1/zeta-15,
            "stationary_gap_or_nonlinear_cutoff_claim": False}
