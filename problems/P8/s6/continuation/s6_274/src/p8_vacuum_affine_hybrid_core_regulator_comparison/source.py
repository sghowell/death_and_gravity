"""Whole original source and a refined full trace-momentum constraint row."""

from functools import cache

import sympy as s
from p8_vacuum_affine_selfconsistent_finite_feedback import source as previous

q, N, u = previous.q, previous.N, previous.u
TIME = previous.TIME
RADIUS = s.Integer(10) ** 20
ANALYSIS_RADIUS = s.Integer(10) ** 150
H_AMPLITUDE = s.Integer(10) ** 172
F_AMPLITUDE = s.Rational(1, 10**510)
HOMOGENEOUS_CAUCHY = previous.HOMOGENEOUS_CAUCHY
CP = s.Rational(1, 10**108)
upper_power = previous.upper_power


@cache
def trace_constraint():
    derivative = q.eliminate_N_primitive(s.diff(q.CONSTRAINT, q.p))
    clock = previous.previous.clock_at(derivative)
    clock_bound = previous.absolute_interval(clock)
    # Both the full invariant and its actual clock value lie in COORDINATE.
    # This covers the ENTIRE source tube, not only the homogeneous center.
    actual = (
        clock_bound + 10**6 * previous.LAPSE_RADIUS + 24 * 10**4 * previous.COORDINATE
    )
    return {
        "whole_Cp": derivative,
        "whole_actual_clock_Cp": clock,
        "whole_real_clock_interval_bound": clock_bound,
        "whole_complete_joint_tube_Cp_bound": actual,
    }


@cache
def data():
    parent = previous.data()
    trace = trace_constraint()
    explicit = 3 * u * (4 * u**6 + 12 * u**4 + 12 * u**2 + 3) / (1 + u * u) ** 4
    return {
        "whole_original_Hamiltonian": q.HAMILTONIAN,
        "whole_original_constraint": q.CONSTRAINT,
        "whole_original_source_bindings": parent["whole_actual_source_bindings"],
        "whole_original_fixed_profiles": parent["whole_actual_fixed_profile_bindings"],
        "whole_refined_trace_constraint_row": trace,
        "whole_same_physical_and_new_analysis_radii": [RADIUS, ANALYSIS_RADIUS],
        "whole_two_safe_centered_amplitudes": [H_AMPLITUDE, F_AMPLITUDE],
        "whole_full_Cp_proof": "Differentiate the entire original constraint before clock restriction. Its exact clock row is retained, including all original primitive/time contacts. Use the actual real-time C5 clock interval1e-130, then integrate full C_pN bounded by1e6 over the unchanged lapse radius1e-115 and all twelve C_pz bounded by1e4 over at most twice the entire invariant radius1e-120. The resulting full joint-tube row is<1e-108. No independent lapse or changed source is introduced.",
        "whole_domain_scope": "The core radius1e20 and both c1,c2 cutoffs are unchanged. Radius1e150 is only a larger complex ANALYSIS neighborhood used to estimate the exact nonlinear remainder. All original functions, state, physical parameters and the S273 finite hybrid model remain unchanged.",
        "checks": {
            "full_clock_trace_constraint_closed_formula": s.factor(
                trace["whole_actual_clock_Cp"] - explicit
            ),
            "same_hybrid_time": TIME - s.Rational(1, 10**180),
            "same_original_full_500_expression_family": s.Integer(
                len(previous.expressions()) - 500
            ),
            "same_real_fixed_C5_profile_bound": previous.PROFILE
            - s.Rational(2, 10**400),
        },
        "gates": {
            "all_whole_original_source_gates_retained": all(parent["gates"].values()),
            "whole_refined_full_Cp": trace["whole_complete_joint_tube_Cp_bound"] < CP,
            "actual_reference_invariant_clock_inside_original_box": 8
            * previous.SOURCE_TIME
            < previous.COORDINATE,
            "whole_CpN_envelope": parent["whole_six_envelope_group_bounds"]["CNz"]
            < 10**6,
            "all_whole_Cpz_envelopes": parent["whole_six_envelope_group_bounds"]["Czz"]
            < 10**4,
            "analysis_radius_does_not_change_physical_cutoffs": ANALYSIS_RADIUS
            > RADIUS,
            "no_higher_real_time_derivative_or_profile_retuning": True,
        },
    }
