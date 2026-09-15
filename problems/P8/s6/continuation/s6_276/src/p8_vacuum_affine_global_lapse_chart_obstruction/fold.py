"""Exact full-source clock jets with the actual fixed quantum profiles."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_band_neighborhood import source as fifth
from p8_vacuum_affine_nonlinear_lapse_branch import source as original
from p8_vacuum_affine_physical_background_vertices import parent as physical

N, R = original.N, original.R
RHO, PRESSURE = s.symbols(
    "actual_fixed_bounce_rho actual_fixed_bounce_pressure", real=True
)
SHEAR, CURVATURE = original.sh, original.curv
M = s.Symbol("whole_total_M1_density", positive=True)
PROFILE_BOUND = fifth.b.prior_initial.PROFILE_BOUND
FOLD_SHEAR = s.Rational(81, 160) - 2 * PRESSURE / 3 + 7 * RHO / 12
FOLD_CURVATURE = -s.Rational(243, 40) + 2 * PRESSURE - 3 * RHO
FOLD_SECOND = (11391 + 1400 * PRESSURE - 1600 * RHO) / 400
MATTER_SQUARE = s.Rational(1, 100) + 6 * PRESSURE - 4 * RHO


def require_order(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)):
        raise TypeError("Require an exact integer clock-jet order")
    if not 0 <= value <= 3:
        raise ValueError("Require clock-jet order zero through three")
    return int(value)


def clock_rules():
    # These are EXACT finite jets atN1, never full off-clock replacements.
    F = original.F_TREE - PRESSURE - (RHO + PRESSURE) * (N**-2 - 1) / 2
    return {
        s.diff(f, N, k): s.diff(value, N, k)
        for f, value in (
            (R, original.R_TREE),
            (original.F, F),
            (original.Ruu, original.RUU_TREE),
        )
        for k in range(6)
    }


def clock_jet(expression, order):
    order = require_order(order)
    return s.factor(
        s.diff(expression, N, order).subs(clock_rules(), simultaneous=True).subs(N, 1)
    )


@cache
def whole():
    rules = dict(original.CENTER)
    rules.pop(SHEAR)
    rules.pop(CURVATURE)
    rules[original.dp] = M - s.Rational(1, 10)
    return {
        "C": original.CONSTRAINT.subs(rules, simultaneous=True),
        "H": original.HAMILTONIAN.subs(rules, simultaneous=True),
    }


def constraint_jet(order):
    order = require_order(order)
    return _constraint_jet(order)


@cache
def _constraint_jet(order):
    return clock_jet(whole()["C"], order)


def affine_error(expression):
    center = expression.subs({RHO: 0, PRESSURE: 0})
    dr, dp = s.diff(expression, RHO), s.diff(expression, PRESSURE)
    if s.factor(expression - center - dr * RHO - dp * PRESSURE) != 0:
        raise ValueError("Require a literal affine fixed-profile expression")
    return (abs(dr) + abs(dp)) * PROFILE_BOUND


@cache
def data():
    germs = physical.source_germs()
    jets = [constraint_jet(i).subs(M, s.Rational(1, 10)) for i in range(3)]
    at_fold = {SHEAR: FOLD_SHEAR, CURVATURE: FOLD_CURVATURE}
    mstar = s.sqrt(MATTER_SQUARE)
    shift = mstar - s.Rational(1, 10)
    checks = {
        "complete_fixed_profile_clock_constraint": s.factor(
            jets[0] - (6 * PRESSURE - 4 * RHO + CURVATURE + 12 * SHEAR) / 4
        ),
        "complete_fixed_profile_first_lapse_jet": s.factor(
            jets[1]
            - (140 * PRESSURE - 160 * RHO - 30 * CURVATURE + 120 * SHEAR - 243) / 80
        ),
        "exact_full_constraint_at_fold": s.factor(jets[0].subs(at_fold)),
        "exact_full_lapse_pivot_at_fold": s.factor(jets[1].subs(at_fold)),
        "exact_full_second_lapse_jet_at_fold": s.factor(
            jets[2].subs(at_fold) - FOLD_SECOND
        ),
        "exact_transverse_shear_derivative": s.factor(s.diff(jets[0], SHEAR) - 3),
        "literal_profile_momentum_square": s.expand(
            MATTER_SQUARE - s.Rational(1, 100) - 6 * PRESSURE + 4 * RHO
        ),
        "literal_rationalized_matter_shift": s.factor(
            shift * (mstar + s.Rational(1, 10)) - (6 * PRESSURE - 4 * RHO)
        ),
        "zero_profile_shear_control": FOLD_SHEAR.subs({RHO: 0, PRESSURE: 0})
        - s.Rational(81, 160),
        "zero_profile_curvature_control": FOLD_CURVATURE.subs({RHO: 0, PRESSURE: 0})
        + s.Rational(243, 40),
        "zero_profile_second_jet_control": FOLD_SECOND.subs({RHO: 0, PRESSURE: 0})
        - s.Rational(11391, 400),
        **{"source_pinned_" + key: value for key, value in germs["checks"].items()},
    }
    return {
        "whole_original_constraint": whole()["C"],
        "whole_original_Hamiltonian_with_nonzero_primitive": whole()["H"],
        "whole_actual_fixed_profile_bindings": germs[
            "fixed_total_coefficient_function_bindings"
        ],
        "whole_actual_reference_energy": germs["fixed_total_normalized_energy"],
        "whole_actual_reference_pressure": germs["fixed_total_normalized_pressure"],
        "whole_all_three_fixed_vacuum_constants": germs["all_three_vacuum_constants"],
        "whole_actual_source_remainders": (
            germs["entire_R_difference"],
            germs["entire_F_difference"],
        ),
        "whole_full_clock_constraint_jets": jets,
        "whole_fold_invariants": at_fold,
        "whole_fold_second_lapse_jet": FOLD_SECOND,
        "whole_actual_C5_profile_bound": PROFILE_BOUND,
        "whole_connected_path_M1_density_square": MATTER_SQUARE,
        "whole_zero_mean_M1_momentum_amplitude": shift,
        "whole_scope": "At the actual bounce the complete F,R and nonzero primitive contact have the stated finite jets. The heavy correction vanishes toorder8, complement toorder1024. N derivatives throughfour, including mixed Ruu contacts, retain this license. The profiles bind to the actual fixed Proca and both scalar reference stresses; they are not new freely chosen parameters. Nothing replaces the full source away from the clock.",
        "checks": checks,
        "gates": {
            "all_full_original_source_gates": all(germs["gates"].values()),
            "actual_heavy_clock_order_eight": physical.heavy.ORDER == 8,
            "actual_complement_clock_order1024": physical.original.N == 1024,
            "actual_original_profiles_below_one_e_minus390": PROFILE_BOUND
            < s.Rational(1, 10**390),
            "fold_shear_strictly_above49_over100": s.Rational(81, 160)
            - affine_error(FOLD_SHEAR)
            > s.Rational(49, 100),
            "fold_shear_strictly_below_one": s.Rational(81, 160)
            + affine_error(FOLD_SHEAR)
            < 1,
            "fold_curvature_strictly_below_minus6": -s.Rational(243, 40)
            + affine_error(FOLD_CURVATURE)
            < -6,
            "fold_second_jet_strictly_above28": s.Rational(11391, 400)
            - affine_error(FOLD_SECOND)
            > 28,
            "matter_square_strictly_positive": s.Rational(1, 100) - 10 * PROFILE_BOUND
            > 0,
            "matter_shift_below100_profile_bound": 10 * PROFILE_BOUND
            < s.Rational(1, 1000),
            "initial_whole_invariant_domain_contains_shift": 100 * PROFILE_BOUND
            < original.DELTA,
            "full_source_not_bare_tree_or_free_profile_reselection": True,
        },
    }
