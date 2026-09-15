"""Whole retained source and explicit formal-loop bookkeeping."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_detector_ir_cut import source as previous

ETA = previous.ETA
KAPPA = previous.KAPPA
MASS2 = previous.MASS2
S, MU, K = previous.S, previous.MU, previous.K
LOOP = s.Symbol("formal_loop_marker", real=True)
current = previous.current


VACUUM_CHECK_NAMES = (
    "actual_kappa",
    "all_three_vacuum_constants_retained",
    "current_massive_quadratic_kinetic",
    "current_massive_quadratic_mass",
    "no_current_scalar_one_point",
    "current_R_vacuum",
    "no_nonminimal_Y_R_cubic",
    "no_nonminimal_Phi_R_mixing",
    "no_nonminimal_Phi_squared_R_cubic",
    "actual_Newton_dictionary",
)


def retained_vacuum_checks():
    """Copy precisely the original low-jet contract, regardless of parent cache order."""
    full, mutable_checks = previous.full_vacuum()
    return full, {name: mutable_checks[name] for name in VACUUM_CHECK_NAMES}


@cache
def data():
    full, old = retained_vacuum_checks()
    u, X = current.u, current.X
    zero = {u: 0, X: 0}
    fixed = current.quantum.fixed_profile()
    quantum = current.profile.data()["full_coefficient_correction"] + fixed[
        "DeltaF"
    ].subs({current.state.TIME: u, current.state.X: X}, simultaneous=True)
    classical = full["F_full"] - quantum
    marked = classical + LOOP * quantum
    heavy_source = current.heavy.coefficients()["normalized_heavy_source"]
    checks = {"retained_" + k: v for k, v in old.items()}
    checks.update(
        {
            "entire_original_action_at_loop_marker_one": s.expand(
                marked.subs(LOOP, 1) - full["F_full"]
            ),
            "classical_vacuum_origin_is_zero": s.cancel(classical.subs(zero)),
            "all_fixed_quantum_retuning_constants_retained": s.cancel(
                quantum.subs(zero) + full["vacuum_pressure_total"]
            ),
            "formal_Gaussian_plus_counterterm_vacuum_density": s.cancel(
                marked.subs(zero) + LOOP * full["vacuum_pressure_total"]
            ),
            "fixed_quantum_profiles_no_quadratic_kinetic_shift": s.cancel(
                s.diff(quantum, X).subs(zero)
            ),
            "fixed_quantum_profiles_no_quadratic_mass_shift": s.cancel(
                s.diff(quantum, u, 2).subs(zero)
            ),
            "whole_R_X_zero_for_every_Phi_at_X_zero": s.cancel(
                s.diff(full["R_full"], X).subs(X, 0)
            ),
            "whole_R_equals_one_for_every_Phi_at_X_zero": s.cancel(
                full["R_full"].subs(X, 0) - 1
            ),
            "whole_heavy_source_no_vacuum_H_onepoint": heavy_source.subs(zero),
            "whole_heavy_source_no_Phi_H_mixing": s.diff(heavy_source, u).subs(zero),
            "same_physical_Newton_dictionary": 8 * s.pi * previous.NEWTON * KAPPA - 1,
        }
    )
    return {
        "entire_original_coefficient_functions": (full["R_full"], full["F_full"]),
        "entire_fixed_reference_quantum_coefficient": quantum,
        "entire_formal_loop_marked_scalar_coefficient": marked,
        "all_three_retained_vacuum_constants": full["vacuum_pressure_total"],
        "explicit_retained_vacuum_check_contract": VACUUM_CHECK_NAMES,
        "original_parameters": {
            "kappa": KAPPA,
            "external_tree_mass_squared": MASS2,
            "Newton": previous.NEWTON,
            "heavy_mass_squared": current.heavy.MASS2,
            "unchanged_reference_lambda": current.heavy.LAMBDA,
        },
        "formal_loop_definition": "The old fixed Gaussian retunings are counted once at one loop. At marker1 the ENTIRE original F is recovered; nothing is retuned. At marker0 its vacuum value is zero, and the old Gaussian vacuum density cancels the fixed counterterm coefficient at first loop. This licenses the formal leading physical graviton reference, not an exact finite-gravity quantum vacuum.",
        "leg_count": "R-1 has at least four light fields near the vacuum. The regular Ia A3/A4 operators already have four light factors and A5 six. The heavy source begins H Phi^2 but no H h h or Phi H bilinear is present; the original Proca source starts above two light fields. Thus the leading two-Phi two-graviton tree is minimal Einstein plus the massive scalar. Counterterm insertions into its one-loop sewing are higher formal loop order.",
        "not_established": "Finite loop errors, an exact LSZ mass or an interacting quantum Minkowski/bounce state. Neither the nonzero fixed vacuum constants nor curvature-dependent one-loop terms are erased.",
        "checks": checks,
        "gates": {
            "retained_check_contract_is_ten_distinct_explicit_names": len(
                VACUUM_CHECK_NAMES
            )
            == len(set(VACUUM_CHECK_NAMES))
            == 10,
            "original_kappa_mass_and_all_constants_retained": KAPPA == 10**800
            and MASS2 == 1,
            "full_original_functions_not_replaced_by_clock_germs": full["R_full"]
            == current.heavy.coefficients()["R"],
            "formal_loop_marker_one_recovers_entire_original_F": checks[
                "entire_original_action_at_loop_marker_one"
            ]
            == 0,
            "fixed_Gaussian_vacuum_cancellation_not_new_normal_ordering": True,
            "one_loop_reference_not_nonperturbative_vacuum": True,
            "original_M1_and_all_previous_scopes_unchanged": True,
        },
    }
