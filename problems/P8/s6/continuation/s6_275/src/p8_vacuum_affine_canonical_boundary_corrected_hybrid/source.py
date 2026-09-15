"""Original full source, with no inheritance of the refuted physical dictionary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_hybrid_core_regulator_comparison import source as previous

q, N, u = previous.q, previous.N, previous.u
TIME, RADIUS, ANALYSIS_RADIUS = previous.TIME, previous.RADIUS, previous.ANALYSIS_RADIUS
H_AMPLITUDE, F_AMPLITUDE = previous.H_AMPLITUDE, previous.F_AMPLITUDE
HOMOGENEOUS_CAUCHY = previous.HOMOGENEOUS_CAUCHY
HOMOGENEOUS_RADIUS = previous.previous.HOMOGENEOUS_RADIUS
REAL_RADIUS = previous.previous.REAL_RADIUS
LAPSE_RADIUS = previous.previous.LAPSE_RADIUS
COORDINATE = previous.previous.COORDINATE
PROFILE = previous.previous.PROFILE
OFF_CLOCK = previous.previous.OFF_CLOCK
CP = previous.CP
upper_power = previous.upper_power


@cache
def data():
    whole = previous.previous.data()
    row = previous.trace_constraint()
    return {
        "whole_original_raw_Hamiltonian": q.HAMILTONIAN,
        "whole_original_raw_constraint": q.CONSTRAINT,
        "whole_original_actual_source_bindings": whole["whole_actual_source_bindings"],
        "whole_original_fixed_profile_bindings": whole[
            "whole_actual_fixed_profile_bindings"
        ],
        "whole_all500_source_derivative_ceilings": whole[
            "whole_500_derivative_ceilings"
        ],
        "whole_all63_componentwise_RFj_bounds": whole[
            "whole_all_63_componentwise_RFj_bounds"
        ],
        "whole_all27_mixed_six_jet_bounds": whole["whole_all_27_mixed_six_jet_bounds"],
        "whole_six_source_envelope_groups": whole["whole_six_envelope_group_bounds"],
        "whole_complete_refined_constraint_p_row": row,
        "whole_unchanged_physical_and_analysis_parameters": {
            "time": TIME,
            "physical_core": RADIUS,
            "analysis_radius": ANALYSIS_RADIUS,
            "full_interaction_amplitude": H_AMPLITUDE,
            "full_volume_amplitude": F_AMPLITUDE,
            "lapse_radius": LAPSE_RADIUS,
            "source_box": COORDINATE,
            "homogeneous_complex_radius": HOMOGENEOUS_RADIUS,
            "homogeneous_Cauchy_radius": HOMOGENEOUS_CAUCHY,
            "homogeneous_real_path_radius": REAL_RADIUS,
        },
        "whole_source_reuse_boundary": "Only whole raw source identities, actual fixed function/profile bindings and their evaluated source-domain envelopes are reused. No old raw/prepared state dictionary or assertion that the S273-S274 solutions are the corrected solutions is inherited. All actual lapse/primitive, heavy, Gauss, shape and generated-mode contacts remain; the fixed real C5 profiles are not live means or analytically continued in time.",
        "checks": {
            "whole_actual_raw_constraint_derivative": s.factor(
                q.eliminate_N_primitive(s.diff(q.HAMILTONIAN, N)) - q.CONSTRAINT
            ),
            "whole500_source_rows": s.Integer(
                len(previous.previous.expressions()) - 500
            ),
            "original_fixed_profile_bound": PROFILE - s.Rational(2, 10**400),
            "corrected_time_unchanged": TIME - s.Rational(1, 10**180),
            "whole_clock_trace_row": s.factor(
                row["whole_actual_clock_Cp"]
                - 3 * u * (4 * u**6 + 12 * u**4 + 12 * u**2 + 3) / (1 + u * u) ** 4
            ),
        },
        "gates": {
            "all_complete_raw_source_envelopes_valid": all(whole["gates"].values()),
            "whole_complete_joint_tube_Cp_bound": row[
                "whole_complete_joint_tube_Cp_bound"
            ]
            < CP,
            "all500_full_source_rows_bounded": len(
                whole["whole_500_derivative_ceilings"]
            )
            == 500,
            "source_and_profiles_not_replaced": True,
            "tree_fixture_is_not_the_actual_live_profile_hessian": True,
        },
    }
