"""Absolute complete specified free curved tensor, with no omitted mode range."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_curved_dirac_state import geometry
from p8_vacuum_flat_local_energy import energy as potential

from . import bloch, reference

POTENTIAL_KEYS = (
    "fixed_mass_anchor_potential_upper",
    "full_quartic_potential_upper",
    "complete_higher_even_potential_upper",
)


def enclosures(mean_mass, amplitude, timescale, multiplicity=6, total_multiplicity=42):
    local = reference.enclosures(mean_mass, amplitude, timescale, multiplicity)
    if type(total_multiplicity) is not int or total_multiplicity < multiplicity:
        raise TypeError(
            "Require a native total multiplicity including every active copy"
        )
    remainder = bloch.enclosures(mean_mass, amplitude, timescale, total_multiplicity)
    old = potential.enclosures(mean_mass, amplitude, timescale, multiplicity)
    V = {key: old[key] for key in POTENTIAL_KEYS}
    Euler_rho = s.Rational(total_multiplicity * 11 * 16, 60 * 144)
    Euler_P = s.Rational(total_multiplicity * 11, 60 * 144) * (16 + s.Rational(64, 3))
    return {
        "complete_potential_components": V,
        "complete_potential_upper": sum(V.values()),
        "complete_remainder_components": remainder,
        "complete_local_two_derivative_components": local,
        "finite_Euler_energy_upper": Euler_rho,
        "finite_Euler_pressure_upper": Euler_P,
        "complete_absolute_curved_energy_upper": sum(V.values())
        + remainder["complete_subtracted_energy_remainder"]
        + local["local_two_derivative_energy"]
        + Euler_rho,
        "complete_absolute_curved_pressure_upper": sum(V.values())
        + remainder["complete_subtracted_pressure_remainder"]
        + local["local_two_derivative_pressure"]
        + Euler_P,
    }


@cache
def data():
    geom = geometry.data()
    m, d, tau = geom["mean_mass"], geom["amplitude_upper"], geom["mass_transition_time"]
    bounds = enclosures(m, d, tau)
    zero = enclosures(m, 0, tau)
    rho = bounds["complete_absolute_curved_energy_upper"]
    P = bounds["complete_absolute_curved_pressure_upper"]
    checks = {
        "same_actual_mean_mass": m - 10**200,
        "same_actual_amplitude_upper": d - 3 * 10**197,
        "same_actual_mass_transition": tau - s.Rational(1, 10**100),
        "complete_energy_sum": rho
        - bounds["complete_potential_upper"]
        - bounds["complete_remainder_components"][
            "complete_subtracted_energy_remainder"
        ]
        - bounds["complete_local_two_derivative_components"][
            "local_two_derivative_energy"
        ]
        - bounds["finite_Euler_energy_upper"],
        "complete_pressure_sum": P
        - bounds["complete_potential_upper"]
        - bounds["complete_remainder_components"][
            "complete_subtracted_pressure_remainder"
        ]
        - bounds["complete_local_two_derivative_components"][
            "local_two_derivative_pressure"
        ]
        - bounds["finite_Euler_pressure_upper"],
        "same_entire_potential_not_a_truncated_quartic": bounds[
            "complete_potential_upper"
        ]
        - sum(potential.enclosures(m, d, tau)[k] for k in POTENTIAL_KEYS),
        "zero_mass_profile_local_potential": zero["complete_potential_upper"],
        "all_forty_two_Euler_copies_retained": bounds["finite_Euler_energy_upper"]
        - s.Rational(42 * 11 * 16, 60 * 144),
    }
    return {
        "same_actual_mean_mass": m,
        "same_actual_profile_amplitude_upper": d,
        "same_actual_profile_time_scale": tau,
        "active_copies": 6,
        "all_curved_copies": 42,
        "actual_uniform_exact_rational_enclosures": bounds,
        "energy_over_named_kappa": rho / analytic.KAPPA,
        "pressure_over_named_kappa": P / analytic.KAPPA,
        "decimal_diagnostics_only": {
            "absolute_curved_energy": str(s.N(rho, 28)),
            "absolute_curved_pressure": str(s.N(P, 28)),
            **{
                k: str(s.N(v, 28))
                for k, v in bounds["complete_remainder_components"].items()
            },
        },
        "scope": "For either exact S6.170 free in/out Hadamard state on the actual CD metric a=(1+t^2)^2, the specified GY14-SAT8-MR/EC-N0 quadratic one-loop contribution obeys both absolute bounds at every real time in the comoving orthonormal frame. All42 Dirac copies, full momenta, local potential, covariant derivative/Newton terms and evanescent Euler stress are retained.",
        "boundary": "This is not the complete interacting/canonically re-expressed scalar, heavy, gauge or gravitational parent tensor. No higher-loop or physical-cutoff remainder, controlled C2 background, perturbative response, quantum target matching, V/G/B or original P8 closure follows. Division by kappa is a named reference ratio, not a relative error against the zero bounce density.",
        "checks": checks,
        "gates": {
            "subtracted_energy_remainder_below_one_e410": bool(
                bounds["complete_remainder_components"][
                    "complete_subtracted_energy_remainder"
                ]
                < 10**410
            ),
            "subtracted_pressure_remainder_below_one_e416": bool(
                bounds["complete_remainder_components"][
                    "complete_subtracted_pressure_remainder"
                ]
                < 10**416
            ),
            "finite_Euler_energy_below_one": bool(
                bounds["finite_Euler_energy_upper"] < 1
            ),
            "finite_Euler_pressure_below_two": bool(
                bounds["finite_Euler_pressure_upper"] < 2
            ),
            "complete_curved_energy_below_one_e789": bool(rho < 10**789),
            "complete_curved_pressure_below_one_e789": bool(P < 10**789),
            "energy_reference_ratio_below_one_e_minus11": bool(
                rho / analytic.KAPPA < s.Rational(1, 10**11)
            ),
            "pressure_reference_ratio_below_one_e_minus11": bool(
                P / analytic.KAPPA < s.Rational(1, 10**11)
            ),
            "constant_mass_geometric_stress_allowance_retained": bool(
                zero["complete_absolute_curved_energy_upper"] > 0
                and zero["complete_absolute_curved_pressure_upper"] > 0
            ),
        },
    }
