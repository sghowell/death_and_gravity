"""Joint finite scalar-boson-mass correction to the MS zero-soft slope."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_mixed_quartic import joint


@cache
def data():
    x = s.Symbol("chord_squared", positive=True)
    ab = tuple(
        (s.Rational(a + 2, 2), s.Rational(b, 2)) for a, b in ((3, 1), (1, 3), (2, 2))
    ) + tuple(
        (s.Rational(a, 2), s.Rational(b + 2, 2)) for a, b in ((3, 1), (1, 3), (2, 2))
    )
    g = s.Rational(7, 4)
    checks = {
        "exact_massive_minus_massless_chord": s.factor(
            1 / (x + 1) - 1 / x + 1 / (x * (x + 1))
        ),
        "fractional_chord_power_gap": g - 1 - s.Rational(3, 4),
        "four_fermion_trace_bound": 4 * 2**4 - 64,
        "two_minimum_radius_terms": 2 * 80000 - 160000,
        "all_three_cyclic_words_prefactor": 3 * 64 * 360**2 * 160000 - 3981312000000,
        "rounded_prefactor_positive_slack": 4 * 10**12 - 3981312000000 - 18688000000,
        "mass_power_after_two_loop_integration": 4 - 3 - g + s.Rational(3, 4),
        "Gamma_upper_constant": 4 + s.factorial(3) - 10,
    }
    for i, (a, b) in enumerate(ab):
        checks[f"sunset_total_exponent_{i}"] = a + b - 3
        checks[f"sunset_beta_denominator_{i}"] = a + b + 2 * g - 4 - s.Rational(5, 2)
    return {
        "exponents": ab,
        "gamma": g,
        "exact_joint_constants": tuple(
            s.simplify(joint.sunset_constant(a, b, g)) for a, b in ab
        ),
        "individual_constant_upper": 80000,
        "exact_prefactor": 3981312000000,
        "rounded_prefactor": 4 * 10**12,
        "bound": "|delta scalar MS slope| <=4e12 N Y_hi^2/(Q_lo^2 m^(3/2)).",
        "forest_difference": "Proper MS kinetic/mass/Yukawa UV poles are boson-mass independent and cancel in this difference. The whole-fermion-cycle counterterm has zero p^2 coefficient. The differentiated raw difference is absolutely convergent before any regulator limit.",
        "method": "Bound the individual soft Taylor coefficient of degree two with R=min(sqrt(S_k),sqrt(S_l))/360. Keep the unshifted chord difference exact, then majorize 1/[x(x+1)] by x^(-7/4).",
        "checks": checks,
    }
