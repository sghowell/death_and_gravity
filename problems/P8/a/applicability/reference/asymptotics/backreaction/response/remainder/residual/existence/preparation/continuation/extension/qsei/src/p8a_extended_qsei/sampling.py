"""Recalibrated C1 scattering and proper H2 sampling on the longer solution."""

from functools import cache

import sympy as sp
from p8a_extension import bounds as actual
from p8a_see_qsei import scattering

from . import reference


@cache
def scattering_calibration():
    geo = actual.calibration()["geometry"]
    t, m, h0 = geo["history"], geo["history_uprime_cap"], geo["h_max"]
    h1 = h2 = sp.Rational(1, 3)
    # Invoke only the general C1-potential lemma, with fresh A14 inputs.
    # No A12 short-domain numerical sampling or reference bound is used.
    mode = scattering.majorants(t, m*t, m, h0, h1, h2)
    margins = {"local_hprime": h1-h0**2-geo["u_cap"],
               "local_hsecond": h2-m-2*h0*(geo["u_cap"]+h0**2),
               "strict_modulus_exponent": sp.Rational(1, 2)-mode["modulus_exponent"]}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("the longer-domain C1 scattering inputs failed")
    return {**mode, "strict_margins": margins}


@cache
def calibration():
    data, mode = actual.calibration(), scattering_calibration()
    # The source turns off at L0/2, NOT L/2. L is a safe sampler envelope
    # for every compact interval in (L0/2,L); its exact width is recorded.
    ell = data["length"]
    f0 = ell**2/4
    f1 = (ell+sp.Rational(3, 8)*ell**2)/2
    f2 = 1+ell/2+sp.Rational(117, 448)*ell**2
    weighted = ell/4+sp.Rational(39, 224)*ell**2
    auxiliary = f2+sp.Rational(3, 2)*weighted
    c0, c1, c2 = mode["forward_two_derivatives"]
    p0, p1 = mode["local_one_derivative"]
    back, t, m = mode["backward_coefficient"], mode["history"], mode["derivative_cap"]
    sqrt_span = sp.Rational(1, 1000)
    roots = {
        "forward": 2*sqrt_span*(c0*f2+2*c1*f1+c2*f0),
        "local": 2*sqrt_span*(p0*f1+p1*f0),
        "real_history_Parseval": 8*back*m*f0,
        "history_remainder": 2*back*t*mode["qprime_minus_uprime_inverse_k"]*f0,
        "infrared": 2*mode["infrared_error"]*f0,
    }
    root = auxiliary+sum(roots.values())
    signed_credit = 16*reference.calibration()["rounded_negative_EED_magnitude"]*ell**4
    coefficient = root**2+signed_credit
    margins = {
        "source_free_width": ell-data["source_off_from"],
        "proper_Hdot": sp.Rational(1, 7)-(data["geometry"]["u_cap"]+sp.Rational(1, 2))/4,
        "root_coarsening": 1+sp.Rational(2, 10**6)-root,
        "absolute_coefficient": 2-coefficient,
    }
    if sqrt_span**2 != ell or any(value <= 0 for value in margins.values()):
        raise ValueError("the extended absolute all-sampler coefficient failed")
    return {
        "source_off_from": data["source_off_from"], "future_end": ell,
        "source_free_conformal_width": data["source_free_span"],
        "conformal_envelope": ell, "proper_span_upper": 3*ell,
        "sqrt_conformal_envelope": sqrt_span,
        "proper_H_cap": sp.Rational(1, 4), "proper_Hdot_cap": sp.Rational(1, 7),
        "sampler_norms": [f0, f1, f2], "weighted_first_norm": weighted,
        "auxiliary_root": auxiliary, "spectral_error_roots": roots,
        "total_spectral_root": root, "difference_coefficient": root**2,
        "signed_reference_penalty": signed_credit,
        "derived_absolute_coefficient": coefficient, "rounded_absolute_coefficient": sp.Integer(2),
        "strict_margins": margins,
    }


def validate_support_interval(left, right):
    """Validate an enclosing compact conformal interval, relative to x0.

    The actual sampler must have its support in this interval. This
    routine does not inspect a sampler or certify a different geometry.
    """
    left, right = map(scattering.nonnegative, (left, right))
    data = calibration()
    if not data["source_off_from"] < left < right < data["future_end"]:
        raise ValueError("support must be compact inside the actual source-free extension")
    return True


def absolute_bound(proper_second_derivative_norm_squared, *, hbar=1):
    """Magnitude for A14-source-free proper H2 samplers/all Hadamard targets."""
    norm, hbar = map(scattering.nonnegative, (proper_second_derivative_norm_squared, hbar))
    return hbar*calibration()["rounded_absolute_coefficient"]*norm/(16*sp.pi**2)


def identities():
    scale, ell, d2, eed = sp.symbols("T0 ell D2 C", positive=True)
    # ||f''||²=T0^-3 D2² and ||f||²<=T0 ell^4 D2².
    signed = eed/scale**4*(scale*ell**4*d2**2)
    normalized = (16*eed*ell**4)/16*(d2**2/scale**3)
    root, penalty = sp.symbols("root penalty", nonnegative=True)
    return {
        "proper_reference_penalty_rescaling": sp.simplify(signed-normalized),
        "absolute_difference_plus_signed_reference": sp.expand((root**2+penalty)*d2**2-root**2*d2**2-penalty*d2**2),
    }
