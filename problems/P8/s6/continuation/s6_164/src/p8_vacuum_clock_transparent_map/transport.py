"""Named vacuum transport and an exact kinematic clock argument."""

from functools import cache

import sympy as s
from p8_vacuum_protected_yukawa_profile import calibration as profile_calibration
from p8_vacuum_protected_yukawa_profile import transport as profile_transport

from . import calibration, counting, gate


@cache
def data():
    old = profile_transport.data()
    p = profile_calibration.data()
    named = dict(old["identical_named_coefficient_bounds"])
    old_classical = named.pop("classical_full_analytic_target_match_coefficient")
    q = gate.data()["gate"]
    a = calibration.data()["clock_transparent_full_analytic_target_match_coefficient"]
    derivatives = p["diagnostic_linear_argument_mass_derivative_ratios"]
    return {
        "identical_named_vacuum_coefficient_bounds": named,
        "old_map_classical_full_target_error": old_classical,
        "new_map_classical_full_target_error": a,
        "pointwise_paired_mass_enclosure": p["strict_pointwise_paired_mass_enclosure"],
        "on_unit_clock_mass_derivative_ratio_bounds": derivatives,
        "clock_argument_dictionary": "For the covariant canonical clock Psi=sqrt(kappa)t in cosmic time, physical X=(partial Psi)^2/kappa=1. Therefore Fhat=Psi along the whole trajectory, and its time derivatives are exactly those of Psi. The protected Yukawa argument is now the same direct argument used in the S6.163 mass-derivative estimate, independently of the old cubic map's extrapolation.",
        "vacuum_prescription": "Use the full SAT8 action, all fixed counterterms and Jacobian under Fhat, with physical source J Fhat and unchanged J_H H. The gate is polynomial and the map difference starts at degree19. The lower source/action vertices and all named through-two-loop vacuum coefficients coincide with the previous map. No independent ordinary-Psi normalization is inferred.",
        "background_limit": "This is only a kinematic clock/source dictionary. The ordinary canonical polynomial parent does not thereby acquire the target bounce as a solution. No scalar/gravity parent matching, nonlinear invertibility, full Dirac state, particle-production error or interaction cutoff is supplied.",
        "checks": {
            "all_named_vacuum_bounds_unchanged": sum(
                v != old["identical_named_coefficient_bounds"][k]
                for k, v in named.items()
            ),
            "clock_gate_zero": q.subs(gate.X, 1),
            "clock_first_gate_derivative_zero": s.diff(q, gate.X).subs(gate.X, 1),
            "vacuum_cubic_map_unchanged": q.subs(gate.X, 0) - 1,
            "new_classical_error_same_declared_composition": a
            - calibration.data()[
                "clock_transparent_full_analytic_target_match_coefficient"
            ],
            "first_clock_mass_ratio_same_now_matched_argument": derivatives[
                "first_mass_derivative_over_mass_squared_upper"
            ]
            - p["diagnostic_linear_argument_mass_derivative_ratios"][
                "first_mass_derivative_over_mass_squared_upper"
            ],
            "second_clock_mass_ratio_same_now_matched_argument": derivatives[
                "second_mass_derivative_over_mass_cubed_upper"
            ]
            - p["diagnostic_linear_argument_mass_derivative_ratios"][
                "second_mass_derivative_over_mass_cubed_upper"
            ],
            "clock_variations_through_seven": counting.data()[
                "clock_variations_identical_through_order"
            ]
            - 7,
        },
        "bounds": {
            "named_total_b2_relative_below_one_e_minus_6": bool(
                named["total_b2_relative_upper"] < s.Rational(1, 10**6)
            ),
            "named_second_cut_relative_below_one_e_minus_407": bool(
                named["second_elastic_cut_relative_upper"] < s.Rational(1, 10**407)
            ),
            "new_common_class_full_target_match_below_one_e_minus_800": bool(
                a < s.Rational(1, 10**800)
            ),
            "unit_clock_first_mass_ratio_below_one_e_minus_102": bool(
                derivatives["first_mass_derivative_over_mass_squared_upper"]
                < s.Rational(1, 10**102)
            ),
            "unit_clock_second_mass_ratio_below_one_e_minus_201": bool(
                derivatives["second_mass_derivative_over_mass_cubed_upper"]
                < s.Rational(1, 10**201)
            ),
        },
    }
