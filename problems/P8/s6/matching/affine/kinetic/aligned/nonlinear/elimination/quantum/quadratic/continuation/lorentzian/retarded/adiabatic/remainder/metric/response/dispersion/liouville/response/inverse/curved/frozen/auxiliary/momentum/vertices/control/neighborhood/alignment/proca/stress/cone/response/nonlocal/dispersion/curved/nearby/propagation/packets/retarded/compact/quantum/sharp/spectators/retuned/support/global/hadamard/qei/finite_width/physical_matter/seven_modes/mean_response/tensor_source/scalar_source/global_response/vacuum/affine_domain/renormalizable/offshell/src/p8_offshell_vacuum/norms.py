"""Exact coefficient majorants on a stated finite canonical jet cube."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import calibration

from . import jets


def coefficient_groups(expr):
    lam, gamma, c = sp.symbols("quartic_lambda quartic_gamma redundant_c", real=True)
    fields = tuple(sorted(expr.free_symbols.intersection(jets.BY_SYMBOL), key=str))
    result = {lam: sp.Integer(0), gamma: sp.Integer(0), c: sp.Integer(0)}
    for _, coefficient in sp.Poly(expr, *fields).terms():
        poly = sp.Poly(coefficient, lam, gamma, c)
        for powers, number in poly.terms():
            if sum(powers) != 1:
                raise ValueError(
                    "Require exactly linear actual quartic coefficient groups"
                )
            result[(lam, gamma, c)[powers.index(1)]] += abs(number)
    return result


def actual_norm(expr):
    d = calibration.point()
    lam, gamma, c = sp.symbols("quartic_lambda quartic_gamma redundant_c", real=True)
    values = {lam: d["lambda"], gamma: d["gamma"], c: d["cubic_squared"] / d["D"] ** 2}
    return sum(value * values[key] for key, value in coefficient_groups(expr).items())


@cache
def data():
    d = jets.data()
    R = d["cubic_field_redefinition"]
    L = d["literal_centered_resolvent_truncation"]
    C = actual_norm(R)
    CL = actual_norm(L)
    # Every monomial in R has exactly three field-jet factors. Each
    # coordinate differentiation increases its l1 coefficient norm by at
    # most three, including repeated derivative indices.
    C6 = 3**6 * C
    free = sp.Rational(37, 2) * C * C
    quartic = CL * ((1 + C6) ** 4 - 1)
    total = free + quartic
    B = sp.Symbol("canonical_jet_radius", nonnegative=True)
    general = sp.Rational(37, 2) * C * C * B**6 + CL * B**4 * ((1 + C6 * B**2) ** 4 - 1)
    derivative_checks = {}
    for axis in range(4):
        value = jets.derivative(R, axis)
        derivative_checks["first_derivative_" + str(axis) + "_majorant"] = (
            actual_norm(value) <= 3 * C
        )
    return {
        "field_map_coefficient_groups": {
            str(k): v for k, v in coefficient_groups(R).items()
        },
        "quartic_action_coefficient_groups": {
            str(k): v for k, v in coefficient_groups(L).items()
        },
        "actual_C_R": C,
        "actual_C_truncated_quartic": CL,
        "all_first_six_derivatives_majorant": C6,
        "unit_jet_cube_free_sixth_order_error": free,
        "unit_jet_cube_quartic_substitution_error": quartic,
        "unit_jet_cube_total_field_redefinition_error": total,
        "general_radius_pointwise_error": general,
        "jet_domain": "Every independent canonical field derivative of order <=10 has absolute value <=B; actual map jets needed by the sixth-derivative quartic action then change by <=3^6 C_R B^3. This is a pointwise bound, not global Sobolev invertibility.",
        "checks": {
            "field_redefinition_is_homogeneous_cubic": sp.expand(
                sum(
                    v * sp.diff(R, v)
                    for v in R.free_symbols.intersection(jets.BY_SYMBOL)
                )
                - 3 * R
            ),
            "original_truncation_is_homogeneous_quartic": sp.expand(
                sum(
                    v * sp.diff(L, v)
                    for v in L.free_symbols.intersection(jets.BY_SYMBOL)
                )
                - 4 * L
            ),
            "free_action_remainder_majorant_constant": sp.Rational(37, 2)
            - (4 * 3**2 + 1) / sp.Integer(2),
        },
        "bounds": {
            "field_map_majorant_below_one_e_minus_403": 0 < C < sp.Rational(1, 10**403),
            "first_six_map_derivatives_below_one_e_minus_400": 0
            < C6
            < sp.Rational(1, 10**400),
            "unit_jet_free_error_below_one_e_minus_805": 0
            < free
            < sp.Rational(1, 10**805),
            "unit_jet_total_error_below_one_e_minus_800": 0
            < total
            < sp.Rational(1, 10**800),
            **derivative_checks,
        },
    }
