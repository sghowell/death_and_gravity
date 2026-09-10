"""Exact finite mass-ratio correction and a uniform compact-integral bound."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Yhi, Qlo):
    m, Yhi, Qlo = map(rational, (m, Yhi, Qlo))
    if m < 2 or Yhi < 0 or not 0 < Qlo <= 144:
        raise ValueError("Need m>=2, Y>=0 and conservative 0<Qlo<=144")
    T = 4 * m * m
    n = max(0, int(T.p).bit_length() - int(T.q).bit_length())
    while sp.Integer(2) ** n < T:
        n += 1
    while n > 0 and sp.Integer(2) ** (n - 1) >= T:
        n -= 1
    delta = (20 + 2 * n) / T
    pref = 12 * Yhi / Qlo**2
    return {
        "spectral_threshold": T,
        "mass_ratio": 1 / T,
        "least_dyadic_exponent": n,
        "dimensionless_finite_correction_absolute_upper": delta,
        "finite_F_absolute_upper": pref * (4 + delta),
        "finite_F_correction_absolute_upper": pref * delta,
        "C_over_Q_upper": pref,
        "scope": "Finite outer MS reference only, with its full small mass-ratio correction retained.",
    }


@cache
def data():
    e = sp.Symbol("epsilon", real=True)
    k, z, r = sp.symbols("k z r", positive=True)
    A = 1 / (1 - k) ** 2
    h = A - 1
    L = -k * sp.log(k) / (1 - k)
    J = (1 - k ** (1 - e)) / ((1 - k) * (1 - e))
    derivative = sp.diff(
        z ** (2 * e) * (1 - z) ** (-e) * (A * J - 1 / (1 - e)), e
    ).subs(e, 0)
    expected = (2 * sp.log(z) - sp.log(1 - z) + 1) * h - A * L
    finite = (
        (1 - z) ** sp.Rational(3, 2)
        / z
        * ((3 - 4 * sp.log(2) + 2 * sp.log(z) - sp.log(1 - z)) * h - A * L)
    )
    real_log_primitive = (1 - z) * sp.log(1 - z) + z
    checks = {
        "outer_Feynman_parameter_at_zero_regulator": sp.factor(J.subs(e, 0) - 1),
        "outer_parameter_regulator_derivative": sp.simplify(
            sp.diff(J, e).subs(e, 0) - (1 - L)
        ),
        "nonzero_difference_pole_weight": sp.factor(A * J.subs(e, 0) - 1 - h),
        "full_difference_first_derivative": sp.simplify(derivative - expected),
        "exact_rational_pole_weight": sp.factor(h - (2 * k - k * k) / (1 - k) ** 2),
        "pole_weight_below_three_k": sp.factor(
            3 * k - h - k * (1 - 5 * k + 3 * k * k) / (1 - k) ** 2
        ),
        "pole_weight_positive_gap_on_domain": sp.expand(
            1
            - 5 * k
            + 3 * k * k
            - (sp.Rational(11, 16) + 5 * (sp.Rational(1, 16) - k) + 3 * k * k)
        ),
        "log_weight_denominator_power": sp.factor(
            A * L - (-k * sp.log(k)) / (1 - k) ** 3
        ),
        "log_weight_prefactor_below_two": 2
        - sp.Rational(16, 15) ** 3
        - sp.Rational(2654, 3375),
        "integrable_log_z": sp.integrate(-sp.log(z), (z, 0, 1)) - 1,
        "real_log_one_minus_z_primitive": sp.diff(real_log_primitive, z)
        + sp.log(1 - z),
        "integrable_log_one_minus_z": sp.limit(real_log_primitive, z, 1, dir="-")
        - real_log_primitive.subs(z, 0)
        - 1,
        "finite_first_group_bound": 3 * (3 + 2 + 1) - 18,
        "complete_correction_constant": 18 * r
        + 2 * r * (-sp.log(r) + 1)
        - r * (20 - 2 * sp.log(r)),
        "zero_Yukawa_zero_reference": enclosure(2, 0, 144)["finite_F_absolute_upper"],
    }
    return {
        "outer_parameter_J": J,
        "difference_pole_integrand": "(1-z)^(3/2)/z * h(rz), h(k)=(1-k)^(-2)-1",
        "finite_correction_integrand_before_k_substitution": finite,
        "finite_correction": "Substitute k=r z and integrate z from zero to one. Multiply by C/Q.",
        "analytic_regulator_domain": "The subtracted compact integral is holomorphic for |epsilon|<1/8; its endpoint powers and logarithms have integrable majorants.",
        "uniform_correction_upper": "|(F_MS-F_MS_r0)/(C/Q)| <= r[20-2log r], 0<r<=1/16",
        "scope": "The mass-ratio correction has a nonzero simple-pole coefficient as well as a finite part. Its pole is subtracted, not discarded.",
        "checks": checks,
    }
