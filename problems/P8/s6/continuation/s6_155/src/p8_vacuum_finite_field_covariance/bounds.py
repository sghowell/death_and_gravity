"""Rational control of first finite-field terms, not the new second slope."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_finite_contact.calibration import rational


def enclosure(L, g, M, Y, canonical_one_loop_relative, Q=144):
    L, g, M, Y, E, Q = map(rational, (L, g, M, Y, canonical_one_loop_relative, Q))
    if L <= 0 or g <= 0 or Y <= 0 or E < 0 or not 32 < M <= 10**400 or not 0 < Q <= 144:
        raise ValueError(
            "Need positive L,g,Y, nonnegative E, 32<M<=10^400 and 0<Q<=144"
        )
    N = s.Integer(6)
    r0 = g / (2 * Q * M)
    k0 = r0 + 4 * N * Y / (3 * Q)
    k1 = 1200 * r0 + 4 * N * Y / Q
    amputated = E + 2 * k0
    return {
        "first_normalization_zero_absolute_upper": k0,
        "first_normalization_epsilon_absolute_upper": k1,
        "one_loop_amputated_relative_upper": amputated,
        "first_normalization_second_amplitude_relative_upper": 3 * k0**2
        + 2 * k0 * amputated,
        "coordinate_MS_commutator_b2_relative_upper": L * k1 / Q,
        "first_normalization_tree_square_relative_upper": 3 * k0**2,
    }


@cache
def data():
    m = s.Symbol("m", positive=True)
    b = 1 / (4 * m**2 - 1)
    fp_allowance = b / 3 + b**2 / 4 + s.Rational(2, 3)
    return {
        "fermion_log_upper": b,
        "fermion_first_epsilon_over_4NY_div_Q_upper": fp_allowance,
        "checks": {
            "first_log_moment_prefactor": 2 * s.Rational(1, 6) - s.Rational(1, 3),
            "second_log_moment_prefactor": s.Rational(3, 2) * s.Rational(1, 6)
            - s.Rational(1, 4),
            "Gamma_finite_epsilon_allowance": s.Rational(16, 24) - s.Rational(2, 3),
            "strict_mass_36_denominator": 4 * 36**2 - 1 - 5183,
            "scalar_scale_log_upper": 400 * 3 - 1200,
        },
        "scope": "At mu=m=10^200 and m>=36, ell_x<=1/(4m^2-1). The two log moments and pi<4 give |fp1|<4NY/Q. The scalar coefficient obeys |r1|<=1200g/(2QM). The one-loop amputated coefficient is bounded by the known canonical one plus its single removed field factor. No new second slope or physical truncation is bounded here.",
    }
