"""Exact-optical matching error tradeoff with explicit unproved amplitude premises."""

from functools import cache

import sympy as s
from p8_vacuum_affine_full_flat_tensor_response import analytic as tensor

from . import amplitude as a

AGREEMENT_ENERGY = s.Integer(10) ** 125
ETA = s.Rational(1, 2)
COEFFICIENT_TOLERANCE = s.Integer(1)
COEF = s.Rational(45, 131072)


@cache
def data():
    S = a.S
    v = s.Symbol("nonnegative_s_minus_sixteen", nonnegative=True)
    eta, energy = s.symbols("relative_angular_error com_energy", nonnegative=True)
    K = s.Symbol("annulus_upper_s", positive=True)
    lower = s.Rational(9, 16) * S * (S - 4) ** 2 - 8 - S**3 / 4
    polynomial = s.Poly(s.expand(lower.subs(S, v + 16)), v)
    low = COEF * (1 - eta) ** 2 * a.GAMMA**2 * energy**8 / s.pi**2
    small = (
        4 * a.LAMBDA * AGREEMENT_ENERGY**4
        + 3 * a.GAMMA * AGREEMENT_ENERGY**6 / 4
        + 8 * a.GAMMA
    )
    actual_cut_lower = COEF * (1 - ETA) ** 2 * a.GAMMA**2 * AGREEMENT_ENERGY**8 / 10
    return {
        "explicit_premises": "An actual nongravitational scalar S matrix with physical pole mass1 and canonical positive LSZ normalization, the S6 V crossing/analyticity/positive-unitarity, specified pole subtraction and sufficiently bounded high-energy arc exists. Its pole-subtracted forward amplitude is analytic at the crossing point s2; no unresolved cut crosses that neighborhood. Its physical b2 is half the second derivative of the full pole-subtracted amplitude, not a freely adjustable tree coefficient. No interacting massless cut is left. These premises, including physical mass and residue matching, are not proved for the original parent.",
        "full_amplitude_error_premise": "For every s in[E^2/2,E^2], ||A_exact(s,.)-A_tree(s,.)||_L2(-1,1)<=eta||A_tree(s,.)||_L2, 0<=eta<1. This includes all loops, heavy states, higher operators and matching corrections. It is not assumed as an established error bound.",
        "central_angular_bound": "For s>=16 and |cos(theta)|<=1/2, the complete original massive A_tree>=gamma s^3/4; hence its full angular squared norm is at least gamma^2 s^6/16.",
        "exact_optical_dispersion_lower": low,
        "necessary_error_tradeoff": "For E^2/2>=16, b2/(4lambda)>9(1-eta)^2(E/10^125)^8, using pi^2<10 and the actual fixed lambda,gamma. Consequently |b2-4lambda|<=4lambda*delta entails 1+delta>9(1-eta)^2(E/10^125)^8.",
        "named_matching_disjunction": "At E10^125, eta<=1/2 and b2<=8lambda cannot coexist with the stated exact S-matrix/dispersion premises: the original elastic annulus alone forces b2>9lambda. At least one premise/tolerance fails; this does not determine which.",
        "named_energy": AGREEMENT_ENERGY,
        "named_angular_relative_tolerance": ETA,
        "named_b2_relative_tolerance": COEFFICIENT_TOLERANCE,
        "named_strict_cut_lower": actual_cut_lower,
        "small_original_tree_amplitude_upper": small,
        "tensor_pole_comparison": "S230 has |lambda_pole(0)|>10^398 in the same anchored units. This is far above both the nominal scalar real-part ceiling and the named conditional agreement energy. No uniform quantum decoupling at those energies, modification of the tensor poles or physical cutoff is inferred from this numerical scale comparison.",
        "checks": {
            "central_massive_positive_polynomial": s.expand(
                lower - polynomial.as_expr().subs(v, S - 16)
            ),
            "annulus_fourth_moment": s.integrate(S**3, (S, K / 2, K)) - 15 * K**4 / 64,
            "identical_optical_and_crossing_coefficient": s.Rational(2)
            * s.Rational(3, 4)
            / 64
            / s.Integer(16)
            * s.Rational(15, 64)
            - COEF,
            "actual_named_cut_is_nine_lambda": actual_cut_lower - 9 * a.LAMBDA,
            "general_dimensionless_tradeoff": s.cancel(
                COEF * a.GAMMA**2 * AGREEMENT_ENERGY**8 / (10 * 4 * a.LAMBDA) - 9
            ),
            "tensor_pole_lower_scale_source": tensor.data()[
                "actual_q0_lambda_modulus_bounds"
            ][0]
            - s.Integer(10) ** 398,
        },
        "gates": {
            "central_massive_polynomial_has_strict_positive_coefficients": all(
                c > 0 for c in polynomial.all_coeffs()
            ),
            "beta_greater_than_three_quarters_at_s16": 1 - s.Rational(4, 16)
            > s.Rational(3, 4) ** 2,
            "annulus_above_exact_massive_threshold": AGREEMENT_ENERGY**2 / 2 > 16,
            "named_cut_exceeds_b2_tolerance": actual_cut_lower > 8 * a.LAMBDA,
            "original_tree_still_tiny_at_named_dispersion_energy": small
            < s.Rational(1, 10**46),
            "conditional_matching_disjunction_not_physical_cutoff": True,
        },
    }
