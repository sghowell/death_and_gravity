"""Constraint, matter-variation and vacuum-domain boundaries of the decomposition."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_clock_transparent_map import norms
from p8_vacuum_flat_dirac_hadamard.symbols import rational


def normalized_first_jet_gap(kappa):
    k = rational(kappa)
    if k <= 1:
        raise ValueError("The compared canonical-clock normalization must exceed one")
    return 1 - 1 / s.sqrt(k)


@cache
def data():
    h = s.Symbol("h", positive=True)
    K, V = s.symbols("K V", real=True)
    B = -s.Rational(1, 2)
    fx = -1 / (2 * h)
    A3 = 1 / h
    A4 = -1 / h - 7 / (4 * h * h)
    A5 = 1 / (h * h)
    C = 4 * fx + A3
    D = A3 + A4 + A5
    L = s.Rational(2, 3) * B * K * K + C * K * V + D * V * V
    Hess = s.hessian(L, (K, V))
    drop = s.hessian(s.Rational(2, 3) * B * K * K + C * K * V + A3 * V * V, (K, V))
    b, N, w, p = s.symbols("log_a lapse chidot momentum", real=True)
    volume = s.exp(3 * b)
    matter = volume * w * w / (2 * N)
    wsol = N * p / volume
    early = s.simplify(matter.subs(w, wsol))
    routh = s.simplify(early - p * wsol)
    rho = lambda L: s.simplify(-s.diff(L, N).subs(N, 1) / volume)
    pressure = lambda L: s.simplify(s.diff(L, b).subs(N, 1) / (3 * volume))
    correct_rho = s.simplify(rho(matter).subs(w, p / volume))
    correct_P = s.simplify(pressure(matter).subs(w, p / volume))
    k = analytic.KAPPA
    gap = normalized_first_jet_gap(k)
    return {
        "full_clock_trace_lapse_velocity_hessian": Hess,
        "deleted_A4_A5_velocity_hessian": drop,
        "deleted_A4_A5_hessian_determinant": s.factor(drop.det()),
        "correct_original_M1_energy": correct_rho,
        "correct_original_M1_pressure": correct_P,
        "incorrect_early_matter_substitution_energy": rho(early),
        "incorrect_early_matter_substitution_pressure": pressure(early),
        "correct_fixed_momentum_Routhian": routh,
        "same_vacuum_matching_class": norms.data()["common_class"],
        "actual_vacuum_class_gradient_invariant_upper": 4 / k,
        "actual_normalized_clock_first_jet_distance_lower": gap,
        "constraint_warning": "A4,A5 have zero first background metric variation but nonzero lapse-velocity quadratic terms. With the actual clock coefficients their contribution is necessary for the Ia trace/lapse degeneracy. Deleting them changes the Hessian rank; this is not an allowed way to simplify the propagating theory.",
        "matter_warning": "Vary the metric while holding the scalar field fixed, before substituting its conserved-momentum solution. Direct substitution into the Lagrangian flips both energy and pressure. A fixed-momentum Routhian includes the Legendre term and restores the correct signs.",
        "domain_warning": "Every field in the frozen flat Schwartz/Fourier unit class has |partial_t Psi|<=1. The canonical clock has partial_t Psi=sqrt(kappa), so its normalized pointwise first-jet distance from this class is at least1-1/sqrt(kappa)>1/2, even on a finite interval. The class also has |X|<=4/kappa, unlike clock X=1. The small common-class action norm does not bound clock actions, equations, or variations outside this class.",
        "checks": {
            "actual_trace_lapse_C": s.factor(C + 1 / h),
            "actual_trace_lapse_D": s.factor(D + 3 / (4 * h * h)),
            "full_Ia_hessian_determinant": s.factor(Hess.det()),
            "trace_lapse_complete_square": s.factor(
                L + (K / s.sqrt(3) + s.sqrt(3) * V / (2 * h)) ** 2
            ),
            "deleted_A4_A5_determinant": s.factor(
                drop.det() + (4 * h + 3) / (3 * h * h)
            ),
            "matter_canonical_momentum": s.diff(matter, w) - volume * w / N,
            "correct_original_matter_energy": s.simplify(
                correct_rho - p * p / (2 * volume**2)
            ),
            "correct_original_matter_pressure": s.simplify(
                correct_P - p * p / (2 * volume**2)
            ),
            "early_substitution_energy_sign_flip": s.simplify(rho(early) + correct_rho),
            "early_substitution_pressure_sign_flip": s.simplify(
                pressure(early) + correct_P
            ),
            "Routhian_energy_restores_fixed_momentum_variation": s.simplify(
                rho(routh) - correct_rho
            ),
            "Routhian_pressure_restores_fixed_momentum_variation": s.simplify(
                pressure(routh) - correct_P
            ),
            "actual_canonical_clock_velocity": s.sqrt(k) - 10**400,
            "actual_first_jet_gap": gap - 1 + s.Rational(1, 10**400),
        },
        "gates": {
            "full_velocity_hessian_rank_one": Hess.rank() == 1,
            "deleting_zero_background_terms_changes_rank": drop.rank() == 2,
            "deleted_hessian_determinant_strictly_negative": bool(
                s.factor(drop.det()).is_negative
            ),
            "actual_vacuum_class_X_does_not_reach_clock": bool(4 / k < 1),
            "actual_normalized_first_jet_gap_above_one_half": bool(
                gap > s.Rational(1, 2)
            ),
            "correct_matter_variation_retained_before_source_substitution": True,
        },
    }
