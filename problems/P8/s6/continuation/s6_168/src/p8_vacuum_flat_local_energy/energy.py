"""Uniform exact rational allowances in the explicitly fixed one-loop scheme."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_flat_dirac_production import calibration
from p8_vacuum_superadiabatic_state_energy import energy as previous

from . import frames


def enclosures(mean_mass, amplitude, timescale, multiplicity=6):
    m, d, tau, _floor = previous.parameters(mean_mass, amplitude, timescale)
    if m < 36:
        raise ValueError("Require mF>=36 in the fixed light-mass-one reference units")
    pieces = frames.enclosures(m, d, tau, multiplicity)
    N = multiplicity
    Q = s.Integer(144)
    pieces.update(
        {
            "fixed_mass_anchor_potential_upper": 2 * N * d * d / (3 * Q),
            "full_quartic_potential_upper": 8 * N * d**4 / (3 * Q),
            "complete_higher_even_potential_upper": N
            * d**6
            / (15 * Q * m * m * (1 - (d / m) ** 2)),
            "MSbar_derivative_energy_upper": N * d * d / (Q * tau * tau),
        }
    )
    pieces["complete_absolute_one_loop_energy_upper"] = sum(pieces.values())
    return pieces


@cache
def data():
    old = calibration.data()
    m, d, tau = (
        old["quadratic_subsystem_mean_mass"],
        old["profile_amplitude_rational_upper"],
        old["same_profile_time_scale"],
    )
    pieces = enclosures(m, d, tau)
    remainder = sum(
        v
        for k, v in pieces.items()
        if k.startswith(("state_to", "fourth_to", "first_frame"))
    )
    ratio = d / m
    return {
        "same_quadratic_mean_mass": m,
        "same_mass_profile_amplitude_upper": d,
        "same_clock_profile_time_scale": tau,
        "active_color_flavor_multiplicity": 6,
        "strict_loop_denominator_lower": 144,
        "actual_uniform_exact_rational_enclosures": pieces,
        "complete_exact_state_subtracted_remainder_upper": remainder,
        "absolute_one_loop_energy_over_named_kappa": pieces[
            "complete_absolute_one_loop_energy_upper"
        ]
        / analytic.KAPPA,
        "decimal_diagnostics_only": {k: str(s.N(v, 28)) for k, v in pieces.items()},
        "scope": "For either the same free in- or out-Hadamard state, at every real time, the absolute value of the specified quadratic one-loop local T00 plus its fixed MS/vacuum/saturated-mass references is bounded by the displayed sum. This includes vacuum polarization and the exact state remainder, not just particle occupations. It does not include scalar/heavy/gauge quantum increments, higher fermion loops, pressure, curved response or the physical EFT cutoff.",
        "reference_ratio_warning": "The kappa comparison is a named dimensionful reference ratio, not a relative error against the bounce density, which is zero at H=0.",
        "checks": {
            "same_mass_reference": m - 10**200,
            "same_amplitude_upper": d - 3 * 10**197,
            "same_time_scale": tau - s.Rational(1, 10**100),
            "exact_sum_not_a_selected_subset": sum(
                v
                for k, v in pieces.items()
                if k != "complete_absolute_one_loop_energy_upper"
            )
            - pieces["complete_absolute_one_loop_energy_upper"],
            "zero_profile_constant_vacuum_reference": enclosures(m, 0, tau)[
                "complete_absolute_one_loop_energy_upper"
            ],
            "vacuum_subtracted_remainder_included": remainder
            - sum(
                pieces[k]
                for k in (
                    "state_to_fourth_frame_energy_upper",
                    "fourth_to_first_frame_energy_upper",
                    "first_frame_second_order_Taylor_remainder_upper",
                )
            ),
        },
        "gates": {
            "entire_accepted_parameter_domain_derivative_factor_below_one": bool(
                2 * s.Rational(99, 10000) / (1 - s.Rational(99, 10000))
                + s.Rational(2, 3)
                < 1
            ),
            "actual_absolute_logarithm_bound_below_one_hundredth": bool(
                2 * ratio / (1 - ratio) < s.Rational(1, 100)
            ),
            "derivative_log_plus_two_thirds_below_one": bool(
                2 * ratio / (1 - ratio) + s.Rational(2, 3) < 1
            ),
            "actual_subtracted_state_remainder_below_one_e411": bool(
                remainder < 10**411
            ),
            "actual_complete_absolute_one_loop_energy_below_one_e789": bool(
                pieces["complete_absolute_one_loop_energy_upper"] < 10**789
            ),
            "actual_named_reference_ratio_below_one_e_minus11": bool(
                pieces["complete_absolute_one_loop_energy_upper"] / analytic.KAPPA
                < s.Rational(1, 10**11)
            ),
            "all_nonzero_actual_allowances_positive": all(
                v > 0 for v in pieces.values()
            ),
        },
    }
