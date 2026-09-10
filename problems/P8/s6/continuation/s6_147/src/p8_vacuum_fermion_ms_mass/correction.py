"""Entire finite scalar-mass-ratio correction after the mass integration."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_mixed_quartic.bounds import dyadic
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Y, Q):
    m, Y, Q = map(rational, (m, Y, Q))
    if m < 720 or Y < 0 or not 0 < Q <= 144:
        raise ValueError("Need m>=720, Y>=0 and 0<Q<=144")
    T = 4 * m * m
    nm, nT = dyadic(m * m), dyadic(T)
    pref = 6 * Y**2 / Q**2
    leading = pref * (78 + 32 * nm)
    correction = pref * (175 + 50 * nT) / T
    return {
        "threshold": T,
        "dyadic_m_squared": nm,
        "dyadic_threshold": nT,
        "leading_scalar_mass_difference_absolute_upper": leading,
        "finite_mass_ratio_remainder_absolute_upper": correction,
        "complete_scalar_mass_difference_absolute_upper": leading + correction,
        "scope": "The full MS b=1 minus b=0 mass reference, with finite pole products already retained in the leading part.",
    }


@cache
def data():
    b, T = s.symbols("b T", positive=True)
    B0 = s.integrate(b, (b, 0, 1))
    Blog = s.integrate(-b * s.log(b), (b, 0, 1))
    integral = (300 * B0 + 100 * s.log(T) * B0 + 100 * Blog) / T
    return {
        "unintegrated_reference_correction_bound": "(b/T)[300+100 log(T/b)] in CY/Q units, 0<b<=1.",
        "integrated_dimensionless_correction_bound": integral,
        "analytic_endpoint_scope": "After multiplying the full regulated reference by epsilon^2, the b integral is uniformly analytic near zero; logarithmic coefficients and b^(-epsilon) are integrable. The exact finite remainder has an additional b/T factor.",
        "checks": {
            "integrable_b_moment": B0 - s.Rational(1, 2),
            "integrable_b_log_moment": Blog - s.Rational(1, 4),
            "complete_integrated_correction_bound": s.simplify(
                integral - (175 + 50 * s.log(T)) / T
            ),
            "constant_after_integration": 300 * s.Rational(1, 2)
            + 100 * s.Rational(1, 4)
            - 175,
            "log_after_integration": 100 * s.Rational(1, 2) - 50,
            "zero_Yukawa_zero_mass_difference": enclosure(720, 0, 144)[
                "complete_scalar_mass_difference_absolute_upper"
            ],
        },
    }
