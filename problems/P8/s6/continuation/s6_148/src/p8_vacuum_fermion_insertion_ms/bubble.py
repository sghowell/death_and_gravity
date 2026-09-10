"""General heavy-mass bubble MS anchor, retaining the physical light denominator."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_mixed_quartic.bounds import dyadic
from p8_vacuum_fermion_outer_ms import laurent
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Y, M, Q):
    m, Y, M, Q = map(rational, (m, Y, M, Q))
    T = 4 * m * m
    if m < 2 or Y < 0 or not 1 <= M <= T / 16 or not 0 < Q <= 144:
        raise ValueError("Need m>=2,Y>=0,1<=M<=T/16 and 0<Q_lower<=144")
    n = dyadic(T / M)
    delta = 18 / T + 2 * M * (n + 1) / T
    pref = 12 * Y / Q**2
    return {
        "threshold": T,
        "ratio_M_over_T": M / T,
        "dyadic_T_over_M": n,
        "finite_bubble_correction_absolute_upper": pref * delta,
        "complete_finite_bubble_absolute_upper": pref * (4 + delta),
        "uniform_derivative_per_unit_g_upper": 4 * pref / T,
    }


@cache
def data():
    e = s.Symbol("epsilon", real=True)
    z, k, r = s.symbols("z k r", positive=True)
    A = (1 - k) ** (-2)
    h = A - 1
    J = (1 - r ** (1 - e)) / ((1 - r) * (1 - e))
    L = -r * s.log(r) / (1 - r)
    derivative = s.diff(z ** (2 * e) * (1 - z) ** (-e) * (A * J - 1 / (1 - e)), e).subs(
        e, 0
    )
    expected = (2 * s.log(z) - s.log(1 - z) + 1) * h - A * L
    return {
        "general_mass_regulated_parameter_factor": J,
        "exact_leading_finite_MS_reference": laurent.data()[
            "leading_finite_MS_coefficient"
        ],
        "finite_difference": "Integral_0^1 (1-z)^(3/2)/z {[3-4log2+2logz-log(1-z)]h(z/T)-A(z/T)L(Mz/T)} dz, A(k)=(1-k)^(-2),h=A-1,L(k)=-klogk/(1-k). Multiply by C/Q.",
        "difference_simple_pole": "Integral_0^1 (1-z)^(3/2)/z h(z/T) dz /epsilon; unchanged by M and not omitted.",
        "uniform_finite_difference_bound": "18/T+2(M/T)[log(T/M)+1] for 1<=M<=T/16.",
        "uniform_derivative_bound": "For 0<=s<=1, |dF_MS(M;s)/ds|<=4C/(QT), same bound as the inherited finite positive on-shell slope. The local MS poles do not depend on s.",
        "checks": {
            "parameter_factor_at_e_zero": s.factor(J.subs(e, 0) - 1),
            "parameter_factor_derivative": s.simplify(
                s.diff(J, e).subs(e, 0) - (1 - L)
            ),
            "full_first_regulator_derivative": s.simplify(derivative - expected),
            "nonzero_pole_independent_of_heavy_mass": s.factor(
                A * J.subs(e, 0) - 1 - h
            ),
            "light_and_heavy_log_denominator_majorant": 2
            - s.Rational(16, 15) ** 3
            - s.Rational(2654, 3375),
            "light_ratio_group_constant": 3 * (3 + 2 + 1) - 18,
            "heavy_log_group_constant": 2 * (1 + 1) - 4,
            "prior_leading_finite_below_four_gap": 4
            - (s.Rational(47, 18) + s.Rational(16, 12))
            - s.Rational(1, 18),
            "same_M_one_remainder_formula": 18
            + 2 * (s.Symbol("log_T") + 1)
            - (20 + 2 * s.Symbol("log_T")),
            "zero_Y_bubble": enclosure(2, 0, 1, 144)[
                "complete_finite_bubble_absolute_upper"
            ],
        },
    }
