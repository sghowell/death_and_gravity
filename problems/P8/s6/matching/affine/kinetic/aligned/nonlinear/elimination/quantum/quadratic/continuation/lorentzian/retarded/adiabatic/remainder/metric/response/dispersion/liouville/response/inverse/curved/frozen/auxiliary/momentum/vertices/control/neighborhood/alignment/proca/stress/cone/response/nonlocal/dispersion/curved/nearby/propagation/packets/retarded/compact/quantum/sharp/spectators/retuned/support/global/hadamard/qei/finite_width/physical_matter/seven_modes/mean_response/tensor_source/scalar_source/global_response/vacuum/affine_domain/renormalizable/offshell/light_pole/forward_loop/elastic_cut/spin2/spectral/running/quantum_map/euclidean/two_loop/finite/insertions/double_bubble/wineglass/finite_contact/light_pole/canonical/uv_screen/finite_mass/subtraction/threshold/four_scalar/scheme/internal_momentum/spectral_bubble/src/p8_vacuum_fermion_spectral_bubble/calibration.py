"""Actual bound for the precisely named first spectral outer-bubble family."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_local_matching import calibration as fermion

from . import outer


@cache
def data():
    f = fermion.data()
    p = model.data()["actual_parameters"]
    m = f["fermion_mass"]
    Y = f["rational_Yukawa_squared_upper"]
    L = p["bare_polynomial_quartic"]
    lam = p["lambda"]
    result = outer.enclosure(m, Y, L, 144)
    E = result["forward_second_coefficient_upper"]
    relative = E / (4 * lam)
    linear = result["channel_subtracted_linear_envelope_coefficient"]
    threshold = result["spectral_threshold"]
    checks = {
        "same_fixed_fermion_mass": m - 10**200,
        "same_spectral_threshold": threshold - 4 * m * m,
        "same_forward_tree_reference": 4 * lam - 4 * sp.Rational(1, 10**600),
        "family_upper_in_exact_reference_couplings": sp.factor(
            E - 6 * Y * L * L / (3 * 144**2 * m**4)
        ),
        "relative_error_reference": relative * (4 * lam) - E,
        "linear_channel_bound_to_second_bound_ratio": sp.factor(
            linear / E - 3 * threshold / 2
        ),
    }
    bounds = {
        "strict_actual_positive_quartic_and_Y": bool(L > 0 and Y > 0),
        "spectral_threshold_strictly_above_sixteen": bool(threshold > 16),
        "positive_literal_family_second_coefficient": result[
            "strictly_positive_family_coefficient"
        ],
        "actual_family_second_coefficient_upper_below_one_e_minus_1418": bool(
            0 < E < sp.Rational(1, 10**1418)
        ),
        "actual_family_relative_upper_below_one_e_minus_819": bool(
            0 < relative < sp.Rational(1, 10**819)
        ),
        "finite_channel_linear_envelope_below_one_e_minus_1017": bool(
            0 < linear < sp.Rational(1, 10**1017)
        ),
    }
    return {
        "actual_reference_parameters": {"mF": m, "Y_upper": Y, "L": L, "lambda": lam},
        "literal_family_enclosure": result,
        "family_second_coefficient_upper_relative_to_positive_tree": relative,
        "bounds": bounds,
        "scope": "The positive order-L^2 Y family is bounded in the same leading reference parameters. This number is not a complete enlarged-model two-loop error budget or a strict V verdict.",
        "checks": checks,
    }
