"""Actual free flat-clock SAT8 transition and energy allowances."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_clock_transparent_map import counting as clock_map
from p8_vacuum_protected_yukawa_profile import calibration as parent

from . import energy, transition


@cache
def data():
    p = parent.data()
    m = p["same_fermion_mass"]
    R = p["fixed_profile_radius"]
    y = p["first_power_Yukawa_rational_upper"]
    delta = y * R
    tau = R / s.sqrt(analytic.KAPPA)
    result = energy.enclosure(m, delta, tau, 6)
    rho = result["complete_free_out_particle_energy_density_upper"]
    ratio = rho / analytic.KAPPA
    return {
        "quadratic_subsystem_mean_mass": m,
        "profile_amplitude_rational_upper": delta,
        "same_profile_time_scale": tau,
        "same_active_color_flavor_multiplicity": 6,
        "complete_quadratic_transition_and_out_energy_enclosure": result,
        "out_energy_upper_divided_by_named_kappa_reference_scale": ratio,
        "state_and_parameter_scope": "The specified free/quadratic Dirac operator uses the SAT8 numerical interaction-boundary parameters and the clock argument aligned by S6.164. This does not identify a higher-loop fermion pole or an interacting gauge-invariant asymptotic particle state.",
        "reference_scale_scope": "rho_out/kappa compares with the named kappa*m_Phi^4 scale in m_Phi=1 units. It is NOT a relative error against the bounce's total energy density, which vanishes at H=0, and is not a semiclassical background-response bound.",
        "checks": {
            "same_mass_reference": m - 10**200,
            "same_actual_profile_amplitude_upper": delta - 3 * 10**197,
            "same_actual_clock_transition_scale": tau - s.Rational(1, 10**100),
            "same_protected_pointwise_mass_floor": result["strict_mass_lower"]
            - p["strict_pointwise_paired_mass_enclosure"]["both_mass_strict_lower"],
            "all_energy_allowances_included": rho
            - result["first_transition_part_of_energy_upper"]
            - result["higher_transition_part_of_energy_upper"],
            "reference_scale_not_bounce_density": ratio * analytic.KAPPA - rho,
            "clock_source_variations_preserved": clock_map.data()[
                "clock_variations_identical_through_order"
            ]
            - 7,
            "zero_momentum_has_zero_transition": transition.amplitude_terms(
                0, result["strict_mass_lower"], delta, tau
            )["exact_transition_amplitude_upper"],
        },
        "bounds": {
            "strict_small_relative_mass_modulation": bool(
                delta / result["strict_mass_lower"] < s.Rational(1, 100)
            ),
            "uniform_mixing_cap_below_one_hundredth": bool(
                result["max_integrated_mixing_upper"] < s.Rational(1, 100)
            ),
            "out_frequency_ratio_below_two": bool(
                m + delta < 2 * result["strict_mass_lower"]
            ),
            "uniform_exact_transition_amplitude_below_one_e_minus_9": bool(
                result["uniform_exact_transition_amplitude_upper"]
                < s.Rational(1, 10**9)
            ),
            "free_out_particle_energy_below_one_e_784": bool(rho < 10**784),
            "named_reference_scale_ratio_below_one_e_minus_16": bool(
                ratio < s.Rational(1, 10**16)
            ),
        },
    }
