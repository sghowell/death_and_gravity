"""Exact proper-clock H2 constants for the new actual-solution interval."""

from functools import cache

import sympy as sp
from p8a_preparation import bounds as prior

from . import scattering


@cache
def calibration():
    data = prior.calibration()
    mode = scattering.calibration()
    ell = data["length"]/2
    # ds=a dx, so the proper envelope is <=3 ell; pi>3 twice.
    f0 = ell**2/4
    f1 = (ell+sp.Rational(3, 8)*ell**2)/2
    f2 = 1+ell/2+sp.Rational(117, 448)*ell**2
    weighted = ell/4+sp.Rational(39, 224)*ell**2
    flat_root = f2+sp.Rational(3, 2)*weighted
    c0, c1, c2 = mode["forward_two_derivatives"]
    p0, p1 = mode["local_one_derivative"]
    back = mode["backward_coefficient"]
    t, m = mode["history"], mode["derivative_cap"]
    # sqrt(ell), sqrt(ell/3)<=10^-5, 2/sqrt(pi)<2,
    # 2sqrt(pi)<4 and sqrt(T)<2. Every replacement is upward.
    sqrt_span = sp.Rational(1, 10**5)
    contributions = {
        "forward": 2*sqrt_span*(c0*f2+2*c1*f1+c2*f0),
        "local": 2*sqrt_span*(p0*f1+p1*f0),
        "real_history_Parseval": 8*back*m*f0,
        "history_remainder": 2*back*t*mode["qprime_minus_uprime_inverse_k"]*f0,
        "infrared": 2*mode["infrared_error"]*f0,
    }
    root = flat_root+sum(contributions.values())
    margins = {"sqrt_span": sqrt_span**2-ell,
               "proper_Hdot": sp.Rational(1, 7)
                   -(data["geometry"]["u_cap"]+sp.Rational(1, 2))/4,
               "root_coarsening": 1+sp.Rational(1, 10**8)-root,
               "rounded_coefficient": 2-root**2}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("the new all-sampler coefficient failed")
    return {"conformal_span": ell, "proper_span_upper": 3*ell,
            "proper_H_cap": sp.Rational(1, 4), "proper_Hdot_cap": sp.Rational(1, 7),
            "sampler_norms": [f0, f1, f2], "weighted_first_norm": weighted,
            "auxiliary_root": flat_root, "spectral_error_roots": contributions,
            "total_root": root, "derived_coefficient": root**2,
            "rounded_coefficient": sp.Integer(2), "strict_margins": margins}


def absolute_bound(proper_second_derivative_norm_squared, *, hbar=1):
    """Bound magnitude for A.11 free-half samplers and all Hadamard targets.

The fixed named prescription and actual-solution domain are hypotheses;
this function does not certify an arbitrary metric or sampler support.
    """
    norm = scattering.nonnegative(proper_second_derivative_norm_squared)
    hbar = scattering.nonnegative(hbar)
    return hbar*calibration()["rounded_coefficient"]*norm/(16*sp.pi**2)


def spectral_identities():
    k, alpha = sp.symbols("k alpha", positive=True)
    a, b = sp.symbols("a b", real=True)
    return {
        "forward_alpha_moment": sp.integrate((alpha+k)**-4, (alpha, 0, sp.oo))-1/(3*k**3),
        "forward_double_moment": sp.integrate(1/(3*k**2), (k, 1, sp.oo))-sp.Rational(1, 3),
        "local_alpha_moment": sp.integrate((alpha+k)**-2, (alpha, 0, sp.oo))-1/k,
        "local_double_moment": sp.integrate(k**-2, (k, 1, sp.oo))-1,
        "backward_remainder_moment": sp.integrate(k**-3, (k, 1, sp.oo))-sp.Rational(1, 2),
        "infrared_radial_moment": sp.integrate(k, (k, 0, 1))-sp.Rational(1, 2),
        "auxiliary_positive_form": sp.expand(a*a+2*b*b+sp.Rational(8, 3)*a*b
            -(a+sp.Rational(4, 3)*b)**2-sp.Rational(2, 9)*b*b),
    }


def clock_identities():
    s = sp.Symbol("s", real=True)
    a, psi = sp.Function("a", positive=True)(s), sp.Function("psi")(s)
    dx = lambda value: a*sp.diff(value, s)
    f = a**-sp.Rational(3, 2)*psi
    h, hubble = sp.diff(a, s), sp.diff(a, s)/a
    derivative = sp.diff(hubble, s)
    return {
        "sampler_zero": sp.simplify(f/sp.sqrt(a)-psi/a**2),
        "sampler_first": sp.simplify(dx(f)/sp.sqrt(a)
            -(sp.diff(psi, s)-sp.Rational(3, 2)*hubble*psi)/a),
        "sampler_second": sp.simplify(dx(dx(f))/sp.sqrt(a)-sp.diff(psi, s, 2)
            +2*hubble*sp.diff(psi, s)-(sp.Rational(3, 4)*hubble**2-sp.Rational(3, 2)*derivative)*psi),
        "weighted_first": sp.simplify(dx(h*f)/sp.sqrt(a)-hubble*sp.diff(psi, s)
            -(derivative-hubble**2/2)*psi),
    }
