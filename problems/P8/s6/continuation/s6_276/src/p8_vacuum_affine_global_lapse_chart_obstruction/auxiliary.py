"""Joint original trace/temporal elimination and a removable sequential pole."""

from functools import cache

import sympy as s
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical as original
from p8_vacuum_affine_nonlinear_reference_volume import source as global_source


@cache
def data():
    N, R, K, T = original.N, original.R, original.K, original.T
    p, G = original.p, original.G
    aa, U, r = original.aa, original.U, original.r
    B, F, lower = original.B, original.F, original.lower
    lag = aa * K * K + B * K + F + U * (T - r * K - lower) ** 2 / 2
    mixed = lag + G * T
    Tfirst = r * K + lower - G / U
    routh = s.factor(mixed.subs(T, Tfirst))
    Kjoint = (p - B - r * G) / (2 * aa)
    Tjoint = s.factor(Tfirst.subs(K, Kjoint))
    rules = {K: Kjoint, T: Tjoint}
    full = N * (p * K - lag - T * G + original.other)
    expected = original.full_trace()["whole_reduced_Hamiltonian"]
    matrix = s.hessian(lag, (K, T))
    A = s.hessian(full, (K, T))
    cross = s.Matrix([s.diff(full, N, K), s.diff(full, N, T)]).subs(
        rules, simultaneous=True
    )
    velocities = s.Matrix([s.diff(Kjoint, N), s.diff(Tjoint, N)])
    schur = (
        s.diff(full, N, 2).subs(rules, simultaneous=True)
        - (cross.T * A.inv() * cross)[0]
    )
    upper = (4 + s.sqrt(7)) / 3
    at_cross = matrix.subs(R, upper)
    packet = global_source.global_R()
    checks = {
        "entire_T_first_Routhian": s.factor(
            routh - (aa * K * K + (B + r * G) * K + F + lower * G - G * G / (2 * U))
        ),
        "entire_T_first_then_K_Hamiltonian": s.factor(
            full.subs(rules, simultaneous=True) - expected
        ),
        "entire_joint_determinant": s.factor(matrix.det() - 2 * aa * U),
        "exact_determinant_value": s.factor(matrix.det() + 2 / (3 * s.sqrt(R))),
        "old_sequential_first_pivot": s.factor(matrix[0, 0] - 2 * aa * original.gamma),
        "new_Routh_velocity_pivot": s.factor(s.diff(routh, K, 2) - 2 * aa),
        "full_stationary_K": s.factor(s.diff(full, K).subs(rules, simultaneous=True)),
        "full_stationary_T": s.factor(s.diff(full, T).subs(rules, simultaneous=True)),
        "entire_first_lapse_envelope": s.factor(
            s.diff(full, N).subs(rules, simultaneous=True) - s.diff(expected, N)
        ),
        "entire_stationary_velocity_contacts": (A * velocities + cross).applyfunc(
            s.factor
        ),
        "entire_second_lapse_Schur": s.factor(schur - s.diff(expected, N, 2)),
        "sequential_Gamma_at_crossing": s.simplify(original.gamma.subs(R, upper)),
        "sequential_first_pivot_at_crossing": s.simplify(at_cross[0, 0]),
        "joint_determinant_at_crossing": s.simplify(
            at_cross.det() + 2 / (3 * s.sqrt(upper))
        ),
        "complete_current_R_on_clock": s.simplify(
            packet["whole_original_R"].subs(global_source.X, 1) - 1
        ),
        "old_pivot_product_cancellation": s.factor(
            (-N * U / original.gamma) * (2 * aa * original.gamma) + 2 * N * aa * U
        ),
    }
    return {
        "whole_original_mixed_lagrangian": mixed,
        "whole_original_Routhian": routh,
        "whole_joint_stationary_solution": s.Matrix([Kjoint, Tjoint]),
        "whole_original_reduced_Hamiltonian": expected,
        "whole_joint_stationary_Hessian": matrix,
        "whole_all_lapse_envelopes": (s.diff(expected, N), s.diff(expected, N, 2)),
        "whole_sequential_crossing_R": upper,
        "whole_actual_global_R_floor": packet["whole_global_X_nonnegative_R_floor"],
        "whole_global_positive_R_argument": packet["whole_global_argument"],
        "whole_scope": "For every positive R the full mixed K,T stationary system is invertible with one positive and one negative eigenvalue. This does NOT prove the remaining lapse Schur pivot is nonzero, a quantum determinant prescription at a rank-changing lapse fold, global spatial reduction or original V/G/B closure.",
        "checks": checks,
        "gates": {
            "whole_U_positive": bool(U.is_positive),
            "whole_joint_negative_determinant": bool((-2 * aa * U).is_positive),
            "whole_actual_R_strictly_above_half": packet[
                "whole_global_X_nonnegative_R_floor"
            ]
            == s.Rational(1, 2),
            "actual_global_R_arguments": all(packet["gates"].values()),
            "actual_R_N_one_third_lower_crosses": s.Rational(9, 2) > upper,
            "joint_crossing_determinant_nonzero": s.simplify(at_cross.det()) != 0,
            "T_first_pivot_at_crossing_nonzero": s.simplify((2 * aa).subs(R, upper))
            != 0,
            "full_other_channels_retained": all(
                expected.has(q)
                for q in (
                    original.pm,
                    original.ph,
                    original.h,
                    original.j,
                    original.curvature,
                    original.shear,
                    original.electric,
                    original.magnetic,
                    original.gm,
                    original.gh,
                    G,
                )
            ),
            "no_Gamma_pole_physical_no_go_inference": True,
        },
    }
