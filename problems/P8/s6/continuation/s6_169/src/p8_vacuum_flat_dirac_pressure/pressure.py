"""Uniform pressure and every homogeneous flat stress component."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_flat_local_energy import energy as prior_energy
from p8_vacuum_flat_local_energy import frames as prior_frames
from p8_vacuum_superadiabatic_state_energy import energy as parameters

POTENTIAL_KEYS = (
    "fixed_mass_anchor_potential_upper",
    "full_quartic_potential_upper",
    "complete_higher_even_potential_upper",
)


def enclosures(mean_mass, amplitude, timescale, multiplicity=6):
    m, d, tau, m0 = parameters.parameters(mean_mass, amplitude, timescale)
    old = prior_energy.enclosures(m, d, tau, multiplicity)
    old_frames = prior_frames.enclosures(m, d, tau, multiplicity)
    N = multiplicity
    Q = s.Integer(144)
    ell = 2 * d / (m - d)
    pieces = {
        "state_and_diagonal_pressure_remainder_upper": sum(old_frames.values()) / 3,
        "linear_off_diagonal_pressure_remainder_upper": N
        * m
        * 1024**3
        * d
        / (27 * tau**4 * m0**2),
        "cubic_off_diagonal_pressure_remainder_upper": 2
        * N
        * m
        * 1024
        * d**3
        / (27 * tau**4 * m0**4),
        "full_reference_subtracted_potential_upper": sum(
            old[key] for key in POTENTIAL_KEYS
        ),
        "restored_MS_derivative_pressure_upper": N
        / (Q * tau * tau)
        * (6 * (m + d) * d * ell + d * d * (ell / 3 + s.Rational(2, 3))),
    }
    pieces["complete_absolute_pressure_upper"] = sum(pieces.values())
    pieces["inherited_absolute_energy_upper"] = old[
        "complete_absolute_one_loop_energy_upper"
    ]
    pieces["every_homogeneous_flat_stress_component_upper"] = max(
        pieces["complete_absolute_pressure_upper"],
        pieces["inherited_absolute_energy_upper"],
    )
    return pieces


@cache
def data():
    old = prior_energy.data()
    m, d, tau = (
        old["same_quadratic_mean_mass"],
        old["same_mass_profile_amplitude_upper"],
        old["same_clock_profile_time_scale"],
    )
    pieces = enclosures(m, d, tau)
    remainder = sum(v for k, v in pieces.items() if "remainder_upper" in k)
    return {
        "same_mass_amplitude_and_time_scale": (m, d, tau),
        "all_actual_uniform_rational_allowances": pieces,
        "complete_subtracted_pressure_remainder_upper": remainder,
        "all_stress_components_over_named_kappa": pieces[
            "every_homogeneous_flat_stress_component_upper"
        ]
        / analytic.KAPPA,
        "decimal_diagnostics_only": {k: str(s.N(v, 30)) for k, v in pieces.items()},
        "homogeneous_stress_tensor": "The same in/out states are rotationally and translationally invariant. Momentum density and off-diagonal stresses vanish by the angular integral. In the fixed Minkowski orthonormal frame <T_mu_nu>=diag(rho,P,P,P); these estimates bound each component, not every Lorentz-boosted component by the same number.",
        "restored_local_pressure_bound": "Use |log(M^2/m^2)|<=ell=2Delta/(m-Delta), |M|<=m+Delta, |Mdot|<=Delta/tau and |Mddot|<=9Delta/tau^2. The paired derivative pressure is bounded by N/[Q tau^2]{6(m+Delta)Delta ell+Delta^2(ell/3+2/3)}, Q>144.",
        "same_finite_mass_reference": "All potential contributions to pressure are minus their energy contributions. Thus the complete S6.168 paired potential, inert vacuum cancellation and fixed -fF(1)chi^2/2 reference use precisely the same entire-tail allowance, with no new energy fit.",
        "curvature_choice": "Fix additional finite curvature counterterms to zero in the common MSbar prescription. This extends the S6.168 T00 calculation, which was insensitive to homogeneous flat curvature improvements. It does not claim that T00 alone had fixed that choice.",
        "scope": "This is the full homogeneous flat expectation tensor of the specified quadratic one-loop Dirac/reference contribution in the named MS field coordinate. It is not the full physical interacting parent stress, a curved state, scalar/heavy/gauge energy, a cutoff or higher-loop remainder, or a bound after the full canonical-field/heavy-background dictionary.",
        "reference_ratio": "The named kappa*m_Phi^4 ratio is not a relative error against zero bounce density. All absolute bounds use light-mass-one units.",
        "checks": {
            "same_mean_mass": m - 10**200,
            "same_amplitude_upper": d - 3 * 10**197,
            "same_clock_time_scale": tau - s.Rational(1, 10**100),
            "complete_pressure_sum": sum(
                pieces[k]
                for k in (
                    "state_and_diagonal_pressure_remainder_upper",
                    "linear_off_diagonal_pressure_remainder_upper",
                    "cubic_off_diagonal_pressure_remainder_upper",
                    "full_reference_subtracted_potential_upper",
                    "restored_MS_derivative_pressure_upper",
                )
            )
            - pieces["complete_absolute_pressure_upper"],
            "inherited_energy_unchanged": pieces["inherited_absolute_energy_upper"]
            - old["actual_uniform_exact_rational_enclosures"][
                "complete_absolute_one_loop_energy_upper"
            ],
            "same_complete_potential": pieces[
                "full_reference_subtracted_potential_upper"
            ]
            - sum(
                old["actual_uniform_exact_rational_enclosures"][key]
                for key in POTENTIAL_KEYS
            ),
            "constant_mass_vacuum_reference_zero": enclosures(m, 0, tau)[
                "every_homogeneous_flat_stress_component_upper"
            ],
            "linear_off_diagonal_bound_mass_dimension": 1 + 1 + 4 - 2 - 4,
            "cubic_off_diagonal_bound_mass_dimension": 1 + 3 + 4 - 4 - 4,
        },
        "gates": {
            "actual_subtracted_pressure_remainder_below_one_e410": bool(
                remainder < 10**410
            ),
            "actual_restored_derivative_pressure_below_one_e595": bool(
                pieces["restored_MS_derivative_pressure_upper"] < 10**595
            ),
            "actual_complete_absolute_pressure_below_one_e789": bool(
                pieces["complete_absolute_pressure_upper"] < 10**789
            ),
            "actual_all_homogeneous_stress_components_below_one_e789": bool(
                pieces["every_homogeneous_flat_stress_component_upper"] < 10**789
            ),
            "actual_reference_ratio_below_one_e_minus11": bool(
                pieces["every_homogeneous_flat_stress_component_upper"] / analytic.KAPPA
                < s.Rational(1, 10**11)
            ),
            "actual_pressure_bound_dominates_the_energy_allowance": bool(
                pieces["complete_absolute_pressure_upper"]
                > pieces["inherited_absolute_energy_upper"]
            ),
        },
    }
