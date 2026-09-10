"""Uniform finite-mass two-gauge cut bound with both crossing orientations."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_vacuum_gauge_yukawa_screen import cuts as leading
from p8_vacuum_gauge_yukawa_screen import threshold

DEFAULT_INVARIANT = sp.Rational(3, 2)
DEFAULT_MASS = sp.Integer(10) ** 200


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational input")
    return sp.Rational(value)


def point(invariant=DEFAULT_INVARIANT, mass=DEFAULT_MASS):
    s, m = rational(invariant), rational(mass)
    if not 1 <= s <= 3:
        raise ValueError("Require the symmetric controlled interval 1<=s<=3")
    if m < 24:
        raise ValueError("Require the strict box routing domain mF>=24")
    eps = 1000000 / m
    D = sp.expand(s * s - (4 - s) ** 2)
    error = 18 * (eps + eps * eps / 8)
    lower = abs(D) - error
    return {
        "s": s,
        "u": 4 - s,
        "fermion_mass": m,
        "Neumann_ratio_upper": 12 / m,
        "per_box_error_over_z_upper": eps,
        "leading_boundary_difference_over_dA_z2_div_9pi": D,
        "absolute_boundary_difference_error_over_dA_z2_div_9pi_upper": error,
        "absolute_boundary_difference_over_dA_z2_div_9pi_lower": lower,
        "relative_boundary_difference_error_upper": error / abs(D) if D else None,
        "strict_nonzero_boundary_difference_at_this_order": bool(lower > 0),
        "scope": "The complete first two-gauge discontinuity, with finite-mass one-loop transitions; not an all-loop amplitude or a physical scalar cross section below s=4.",
    }


@cache
def data():
    eps, z, s = sp.symbols(
        "positive_box_epsilon positive_box_scale real_s", positive=True
    )
    a, Y, m = sp.symbols(
        "positive_gauge_squared positive_Yukawa_squared positive_fermion_mass",
        positive=True,
    )
    q = 16 * sp.pi**2
    color = 8
    product_error = (32 * eps + 4 * eps * eps) * z * z
    phase = color / (32 * sp.pi)
    rho_error = sp.factor(phase * product_error)
    total_error = 2 * rho_error
    rho0 = color * z * z * s * s / (9 * sp.pi)
    raw_relative = sp.factor(total_error / (4 * color * z * z / (9 * sp.pi)))
    inherited = leading.data()
    parent_C = threshold.data()["leading_Phi_squared_F_squared_coefficient"]
    expected_C = -a * Y / (3 * q * m * m)
    y = sp.Symbol("positive_Yukawa", positive=True)
    if parent_C.free_symbols != {a, y, m}:
        raise ValueError("The inherited threshold parameter dictionary changed")
    actual_C = sp.simplify(parent_C.subs(y, sp.sqrt(Y)))
    checks = {
        "inherited_threshold_dictionary": sp.factor(actual_C - expected_C),
        "per_pair_leading_amplitude_normalization": -4 * (-z / 3) * s
        - sp.Rational(4, 3) * z * s,
        "four_polarization_product_error": product_error
        - 4 * (2 * 4 * z * eps * z + (eps * z) ** 2),
        "optical_phase_and_color_factor": phase
        - color * inherited["integrated_massless_two_body_phase_space"] / 4,
        "single_channel_error": sp.factor(
            rho_error - color * z * z / sp.pi * (eps + eps * eps / 8)
        ),
        "single_channel_leading_cut": sp.factor(
            inherited["s_channel_Im_amplitude"].subs(
                {inherited["coefficient"]: -z / 3, inherited["s"]: s}
            )
            - rho0
        ),
        "two_crossed_channel_errors": sp.factor(
            total_error - 2 * color * z * z / sp.pi * (eps + eps * eps / 8)
        ),
        "interior_leading_difference": sp.factor(
            (rho0 - rho0.subs(s, 4 - s)).subs(s, sp.Rational(3, 2))
            + 4 * color * z * z / (9 * sp.pi)
        ),
        "interior_relative_error": sp.factor(
            raw_relative - sp.Rational(9, 2) * (eps + eps * eps / 8)
        ),
        "strict_five_epsilon_margin": sp.factor(
            5 * eps - raw_relative - eps * (8 - 9 * eps) / 16
        ),
        "crossing_antisymmetric_leading_difference": sp.expand(
            s * s - (4 - s) ** 2 - 8 * (s - 2)
        ),
        "leading_crossing_center_cancellation": (s * s - (4 - s) ** 2).subs(s, 2),
        "two_one_loop_transitions_and_two_particle_sewing": 1 + 1 + (2 - 1) - 3,
        "three_particle_sewing_earliest_loop": 1 + 1 + (3 - 1) - 4,
    }
    return {
        "box_scale_z": a * Y / (q * m * m),
        "normalized_box_error_epsilon": 1000000 / m,
        "leading_operator_coefficient": expected_C,
        "leading_parallel_transition_magnitude": sp.Rational(4, 3) * z * s,
        "four_pair_product_absolute_error_upper": product_error,
        "optical_identical_particle_color_phase_factor": phase,
        "leading_single_channel_rho": rho0,
        "single_channel_rho_absolute_error_upper": rho_error,
        "crossed_channel_difference_absolute_error_upper": total_error,
        "interior_relative_error_upper": raw_relative,
        "interior_simplified_relative_error_upper": 5 * eps,
        "first_parent_loop_order": 3,
        "first_transition_inventory": [
            "Six cyclic one-loop fermion boxes, with two active equal-mass opposite-Yukawa flavors",
            "Single-gauge transition vanishes by color; the possible fermion triangle attached to a three-gauge vertex also has a single color trace",
            "No direct H Yukawa or Phi^2 F^2 tree vertex; Phi parity eliminates a Phi single-particle exchange",
            "An H-to-two-gauge transition needs at least two loops; a heavy-H exchange is not a missing one-loop transition",
            "Mass, wavefunction and parameter counterterms enter the transition at the next perturbative order",
            "Two one-loop transitions sewn across two gauge lines give three parent loops; three-gauge final states first need at least four",
        ],
        "scope": "Sew the holomorphic continuations of the two transition factors. Below scalar threshold their product is not replaced by an absolute square of complex external momenta. The s and u boundary prescriptions are opposite. Local counterterms have no discontinuity; this certificate is not the entire three-loop amplitude or its all-order continuation.",
        "checks": checks,
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("s"),
    )
    return (
        [
            (f"inexact_invariant_{j}", point, (v, DEFAULT_MASS))
            for j, v in enumerate(invalid)
        ]
        + [
            (f"inexact_mass_{j}", point, (DEFAULT_INVARIANT, v))
            for j, v in enumerate(invalid)
        ]
        + [
            (f"outside_channel_interval_{j}", point, (v, DEFAULT_MASS))
            for j, v in enumerate((0, sp.Rational(1, 2), sp.Rational(7, 2), 4))
        ]
        + [
            (f"outside_mass_domain_{j}", point, (DEFAULT_INVARIANT, v))
            for j, v in enumerate((-1, 0, 23))
        ]
    )
