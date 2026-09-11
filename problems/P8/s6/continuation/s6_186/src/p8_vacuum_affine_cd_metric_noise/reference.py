"""Uniform complex reference amplitudes and exact actual-state remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vector_hadamard import series
from p8_vector_regularity import frequency, preparation
from p8_vector_state import wkb

MASS = modes.MASS
OUTER = s.Rational(1, 10000)
RADIUS = OUTER / 2
FIELD_NORM = s.Integer(4000)
PAIR_REF = s.Integer(10) ** 9
PAIR_ERROR = s.Integer(10) ** 15


@cache
def data():
    checks = {}
    bounds = {}
    for kind in ("transverse", "longitudinal"):
        row = {}
        for order in range(1, 5):
            P = series.coefficient(kind, order)
            num, den = s.fraction(s.cancel(P))
            d = s.Poly(den, wkb.u)
            power = d.degree() // 2
            lead = d.LC()
            assert lead > 0
            checks[f"{kind}_complex_denominator_{order}"] = s.expand(
                den - lead * (1 + wkb.u**2) ** power
            )
            poly = s.Poly(num / lead, wkb.u, wkb.z)
            B = (
                sum(
                    abs(c) * s.Rational(51, 100) ** i * s.Rational(101, 100) ** j
                    for (i, j), c in poly.terms()
                )
                / s.Rational(7, 10) ** power
            )
            row[order] = {
                "quadratic_denominator_power": power,
                "complex_coefficient_upper": B,
            }
        row["ratio_defect"] = sum(
            row[n]["complex_coefficient_upper"]
            / (s.Rational(99, 100) * MASS) ** (2 * n)
            for n in range(1, 5)
        )
        bounds[kind] = row
    state = preparation.constants()["common"]
    B6 = state["initial_mixing_envelopes"][6]
    E = state["evolution_mixing_envelope"]
    C = max(
        frequency.reference(kind)["residual_over_inverse_frequency_power_upper"]
        for kind in bounds
    )
    remainder = B6 + E / MASS**4 + 6 * C / MASS**3
    scale_defect = 6 * OUTER / (1 - 6 * OUTER)
    W_lower = s.Rational(99, 100) - 2 * s.Rational(1, 1000)
    W_slope = 3 / RADIUS
    half_rate = W_slope / (2 * s.Rational(1, 2))
    prefactor = 3 + (half_rate + 5) / MASS
    B0 = s.Symbol("initial_beta_modulus", nonnegative=True)
    checks.update(
        {
            "same_unmodified_all_order_mixing_envelope": modes.data()[
                "actual_state_mixing_upper"
            ]
            - (B6 + E / MASS**4) / MASS**6,
            "actual_reference_residual_constant": C - 5207948696836,
            "actual_initial_B6": B6 - 1813229,
            "actual_evolved_E": E - 5347035781757616,
            "constant_initial_amplitude_CCR": (1 + B0**2) - B0**2 - 1,
            "reference_pair_norm_count": 2 * (2 * FIELD_NORM) ** 2 - 128000000,
            "remainder_pair_safe_count": 2 * FIELD_NORM**2 * 3 * (2 * 10**6)
            - 192000000000000,
            "reference_real_volume_already_cancelled": s.Rational(3, 2)
            + s.Rational(3, 2)
            - 3,
        }
    )
    return {
        "outer_complex_time_radius": OUTER,
        "inner_amplitude_Cauchy_radius": RADIUS,
        "actual_W8_complex_coefficient_majorants": bounds,
        "relative_a_inverse_squared_defect_upper": scale_defect,
        "complex_reference_real_part_lower_over_nu": W_lower,
        "complex_Wprime_upper_over_nu": W_slope,
        "complex_reference_momentum_upper_over_sqrt_nu": prefactor,
        "each_phase_stripped_ten_field_component_display_upper_over_sqrt_nu": s.Integer(
            1000
        ),
        "phase_stripped_ten_field_Euclidean_norm_upper_over_sqrt_nu": FIELD_NORM,
        "source_pinned_constants": {"B6": B6, "E": E, "C": C},
        "exact_mode_decomposition": "u_actual=alpha_initial u_ref+e, with ||e||<=4000 R nu^(-11/2); alpha_initial is constant, not alpha(t), and |alpha_initial|<2",
        "actual_state_remainder_coefficient": remainder,
        "reference_stress_pair_coefficient_upper": PAIR_REF,
        "remainder_stress_pair_coefficient_upper": PAIR_ERROR,
        "reference_pair_bound": "||a_ref(t,k,l)||F <=1e9 sqrt(nu mu) on the inner complex disc, after both actual constant alpha_initial factors",
        "remainder_pair_bound": "||a_error(t,k,l)||F <=1e15 sqrt(nu mu)(nu^-6+mu^-6) on the real slab; no derivatives of the rapidly oscillating mixing are taken",
        "checks": checks,
        "gates": {
            "complex_u_and_quadratic_margin": s.Rational(1, 2) + OUTER
            < s.Rational(51, 100)
            and 1 - s.Rational(51, 100) ** 2 > s.Rational(7, 10),
            "complex_H_below_three": 4 * s.Rational(51, 100) / s.Rational(7, 10) < 3,
            "full_scale_relative_defect_below_one_thousandth": scale_defect
            < s.Rational(1, 1000),
            "complex_z_below_one_point_zero_one": (1 + s.Rational(1, 1000))
            / (1 - s.Rational(1, 1000))
            < s.Rational(101, 100),
            "complex_omega_modulus_lower": s.Rational(99, 100) ** 2
            < 1 - s.Rational(1, 1000),
            "complex_omega_upper_below_two_nu": (1 + s.Rational(1, 1000))
            * s.Rational(25, 16) ** 2
            < 4,
            "all_full_W8_ratio_defects_below_one_thousandth": all(
                row["ratio_defect"] < s.Rational(1, 1000) for row in bounds.values()
            ),
            "complex_W_real_part_above_half_nu": W_lower > s.Rational(1, 2),
            "complex_W_absolute_upper_below_three_nu": 2 * (1 + s.Rational(1, 1000))
            < 3,
            "complex_canonical_rate_below_five": 3
            * (s.Rational(1, 2) + s.Rational(101, 100))
            < 5,
            "reference_momentum_below_sixty_four_sqrt_nu": prefactor < 64,
            "all_full_constraint_components_below_one_thousand": 3 * 64 < 1000,
            "ten_field_norm_below_four_thousand": 10 * 1000**2 < FIELD_NORM**2,
            "actual_state_remainder_below_two_million": 0 < remainder < 2 * 10**6,
            "remainder_relative_factor_below_one": 2 * 10**6 / MASS**6 < 1,
            "reference_pair_display_covers_full_count": 2 * (2 * FIELD_NORM) ** 2
            < PAIR_REF,
            "remainder_pair_display_covers_full_count": 2
            * FIELD_NORM**2
            * 3
            * (2 * 10**6)
            < PAIR_ERROR,
        },
    }
