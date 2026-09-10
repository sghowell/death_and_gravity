"""Two-region paired kernel estimate with separate legitimate soft radii."""

from functools import cache

import sympy as sp


@cache
def data():
    x, y, m, b, z = sp.symbols("x y m b z", nonnegative=True)
    S = m * m + x
    B1, B0, S1, S0 = sp.symbols("B1 B0 S1 S0", commutative=False)
    c = sp.Symbol("real_q_dot_l", real=True)
    checks = {
        "same_boson_chord_at_both_ends": (x + z) - (y + z) - (x - y),
        "paired_kernel_difference_identity": sp.expand(
            B1 * S1 - B0 * S0 - (B1 * (S1 - S0) + (B1 - B0) * S0)
        ),
        "chord_denominator_difference": sp.factor(
            1 / (x + y - 2 * c + b)
            - 1 / (y + b)
            - (2 * c - x) / ((x + y - 2 * c + b) * (y + b))
        ),
        "high_region_chord_gap": sp.Rational(1, 1)
        - sp.Rational(1, 2)
        - sp.Rational(1, 2),
        "high_region_numerator_rounding": 3 - sp.Rational(5, 2) - sp.Rational(1, 2),
        "high_region_boson_difference_constant": 3 * 4 - 12,
        "two_shifted_fermion_product_difference_constant": 4 * sp.Rational(1, 20)
        + 2 * sp.Rational(1, 20)
        - sp.Rational(3, 10),
        "fermion_product_difference_rounding_margin": 4
        - sp.Rational(3, 10)
        - sp.Rational(37, 10),
        "first_high_difference_bound": 4 * 4 - 16,
        "sum_high_kernel_constants": 16 + 12 - 28,
        "rounded_high_kernel_margin": 32 - 28 - 4,
        "raw_low_six_propagator_trace": 4 * 2**6 - 256,
        "subtraction_low_four_shifted_trace": 4 * 2**4 - 64,
        "local_anchor_four_shifted_trace": 4 * 2**4 - 64,
        "joint_radius_prefix_bound": sp.Rational(18, 360) - sp.Rational(1, 20),
        "m_720_minimum_radius": sp.Rational(720, 360) - 2,
        "log_five_strict_upper_finite_exponential_sum": sum(
            sp.Rational(2**k, sp.factorial(k)) for k in range(4)
        )
        - sp.Rational(19, 3),
        "low_ordered_log_ratio": sp.factor((m * m + 4 * S) / S - (4 + m * m / S)),
        "subtraction_log_upper_argument": sp.expand(5 * S - (m * m + 4 * S) - x),
    }
    return {
        "regions": "LOW: y=l^2 <=4(m^2+q^2); HIGH: y>=4(m^2+q^2), independent of complex external momenta",
        "raw_LOW_Cauchy_radius": "min(sqrt(m^2+q^2),sqrt(m^2+l^2))/360",
        "subtraction_LOW_Cauchy_radius": "sqrt(m^2+q^2)/360, because the subtracted kernel has no external momentum",
        "paired_HIGH_Cauchy_radius": "sqrt(m^2+q^2)/360; y>=4(m^2+q^2) also controls the short fermion arc",
        "HIGH_paired_kernel_norm": "<=32 sqrt(m^2+q^2)/y^(5/2), before scalar/gauge coupling and Dirac caps",
        "HIGH_proof": "Boson norms <=4/y and 1/y, their difference <=12 sqrt(S_q)/y^(3/2); the two-fermion product difference <=4 sqrt(S_q)/S_l^(3/2).",
        "scope": "Never integrate the raw HIGH vertex independently. The paired ultraviolet improvement is essential; in LOW the raw and subtraction terms use different proved Cauchy radii.",
        "checks": checks,
    }
