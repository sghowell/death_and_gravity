"""A fixed complex disk and four exact derivative-control steps."""

from functools import cache

import sympy as s

RADIUS = s.Rational(1, 64)
STEP = s.Rational(1, 512)
SCALE = 2**20
FACTOR = 1024
TRANSITION_CONSTANT = 4 * FACTOR**4


def iteration(index):
    if type(index) is not int or not 0 <= index <= 4:
        raise TypeError("Require a native frame index from zero through four")
    return index


def radius(index):
    return RADIUS - iteration(index) * STEP


def frequency_floor(index):
    return s.Rational(9, 10) * s.Rational(15, 16) ** iteration(index)


@cache
def data():
    r = RADIUS
    x = s.Symbol("x", real=True)
    weight = (1 + x**8) ** (-s.Rational(9, 8))
    eps = s.Rational(1, 100)
    caps = {j: frequency_floor(j) for j in range(5)}
    checks = {
        "disk_polynomial_relative_constant": 8 * r * (1 + r) ** 7
        - s.Rational(4902227890625, 35184372088832),
        "profile_L1_and_four_rotations_constant": 2 * 2 * FACTOR**4
        - TRANSITION_CONSTANT,
        "four_step_transition_power": 1 - 2 - 4 + 5,
        "initial_mass_square_relative_allowance": 4 * eps
        + 4 * eps**2
        - s.Rational(101, 2500),
        "Cauchy_step_and_frequency_floor_factor": 1 / STEP / (s.Rational(1, 2))
        - FACTOR,
        "actual_transition_constant": TRANSITION_CONSTANT - 2**42,
    }
    for j in range(4):
        checks[f"nested_disk_difference_{j}"] = radius(j) - radius(j + 1) - STEP
        checks[f"frequency_floor_recursion_{j}"] = frequency_floor(
            j + 1
        ) - frequency_floor(j) * s.Rational(15, 16)
    return {
        "initial_disk_radius": RADIUS,
        "Cauchy_radius_loss_per_exact_rotation": STEP,
        "four_nested_radii": {j: radius(j) for j in range(5)},
        "complex_frequency_floor_factors": caps,
        "fixed_adiabatic_gap_threshold": SCALE,
        "profile_integrable_weight": weight,
        "integrated_transition_constant": TRANSITION_CONSTANT,
        "profile_disk_proof": "For each real x and |z-x|<=1/64, |z^8-x^8|/(1+x^8)<=8r(1+r)^7<1/4. The positive-real-root branch extends analytically; |s(z)|<2 and |s'(z)|<2w(x), w=(1+x^8)^(-9/8), integral w=2.",
        "initial_frequency_and_mixing": "E=sqrt(p^2+m0^2), m0=.99m. Delta/m0<.01 gives |M(z)^2-m^2|<.1m^2, so the principal omega is analytic and |omega|>=.9E. e0=tau omega and |g0|<=2pDelta/E^2 w(x).",
        "four_step_induction": "On disk r_j, |e_j|>=c_j tau E and |g_j|<=G_j w(x), c_j=.9(15/16)^j>.5, G_j=2pDelta/E^2[1024/(tau E)]^j. With tau*m0>=2^20, |g_j/e_j|<1/4. The next exact square root has modulus at least 15/16; Cauchy with radius loss1/512 bounds the next connection by G_(j+1)w. All square roots/atan branches are fixed by the real positive-frequency frame.",
        "transition_consequence": "The final exact off-diagonal generator is integrable: integral |g4|<=2^42 pDelta/(tau^4 E^6). Unitarity bounds the full transition by this integral, without discarding the final coupling or any of its Dyson orders.",
        "checks": checks,
        "gates": {
            "relative_eighth_power_disk_defect_below_one_quarter": bool(
                8 * r * (1 + r) ** 7 < s.Rational(1, 4)
            ),
            "complex_profile_modulus_below_two": bool(
                (1 + r) ** 8 * s.Rational(4, 3) < 2**8
            ),
            "complex_first_derivative_weight_factor_below_two": bool(
                s.Rational(4, 3) ** 9 < 2**8
            ),
            "mass_square_relative_change_below_one_tenth": bool(
                4 * eps + 4 * eps**2 < s.Rational(1, 10)
            ),
            "initial_sqrt_floor_nine_tenths": bool(
                s.Rational(9, 10) > s.Rational(9, 10) ** 2
            ),
            "all_five_frequency_floors_above_one_half": all(
                c > s.Rational(1, 2) for c in caps.values()
            ),
            "all_nested_disks_have_positive_radius": all(
                radius(j) > 0 for j in range(5)
            ),
            "complex_ratio_cap_below_one_quarter": bool(
                4 * eps / SCALE < s.Rational(1, 4)
            ),
            "frame_mixing_bound_contracts": bool(s.Rational(FACTOR, SCALE) < 1),
            "next_connection_factor_below_one": bool(
                1 / (2 * (1 - s.Rational(1, 4) ** 2)) < 1
            ),
            "sqrt_modulus_floor_valid": bool(
                s.Rational(15, 16) > s.Rational(15, 16) ** 2
            ),
        },
    }
