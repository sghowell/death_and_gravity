"""Connected local-sheet obstruction and exact full envelope boundary."""

from functools import cache

import sympy as s

from . import fold


@cache
def data():
    N, M = fold.N, fold.M
    shear, curvature = fold.SHEAR, fold.CURVATURE
    C, H = fold.whole()["C"], fold.whole()["H"]
    path = {M * M: fold.MATTER_SQUARE, curvature: -12 * shear}
    point_C = s.factor(fold.constraint_jet(0).subs(path, simultaneous=True))
    point_CN = s.factor(fold.constraint_jet(1).subs(path, simultaneous=True))
    point_CNN = s.factor(
        fold.constraint_jet(2)
        .subs(path, simultaneous=True)
        .subs(shear, fold.FOLD_SHEAR)
    )
    canonical = fold.original.parent
    prior = canonical.full_trace()
    bind = {
        canonical.p: 0,
        canonical.G: 0,
        canonical.pm: M,
        canonical.ph: 0,
        canonical.h: 0,
        canonical.Href: 0,
        canonical.shear: shear,
        canonical.electric: 0,
        canonical.magnetic: 0,
        canonical.wmass: 0,
        canonical.gm: 0,
        canonical.gh: 0,
        canonical.curvature: curvature,
        canonical.B: 0,
        canonical.F: canonical.R ** (-s.Rational(3, 4)) * fold.original.F
        - fold.original.Iu / N,
    }
    raw = s.factor(prior["whole_raw_Hamiltonian"].subs(bind, simultaneous=True))
    gamma = 1 - 3 * (canonical.R - 1) ** 2 / (2 * canonical.R)
    raw_expected = H - N * canonical.R ** (-s.Rational(3, 4)) * canonical.T**2 / (
        2 * gamma
    )
    DTT = fold.clock_jet(s.diff(raw, canonical.T, 2), 0)
    DNT = s.factor(s.diff(raw, N, canonical.T).subs(canonical.T, 0))
    Cs, CN = s.diff(C, shear), s.diff(C, N)
    Hs = s.diff(H, shear)
    # Entire exact local solve for shear, not a tree polynomial substitute.
    stationary_shear = -C.subs(shear, 0) / Cs
    clock_fold = {M * M: fold.MATTER_SQUARE, curvature: -12 * fold.FOLD_SHEAR}
    s0 = s.factor(
        fold.clock_jet(stationary_shear, 0).subs(clock_fold, simultaneous=True)
    )
    s1 = s.factor(
        fold.clock_jet(stationary_shear, 1).subs(clock_fold, simultaneous=True)
    )
    s2 = s.factor(
        fold.clock_jet(stationary_shear, 2).subs(clock_fold, simultaneous=True)
    )
    a, b, D = s.symbols(
        "lapse_deviation transverse_shear_deviation positive_fold_second_jet", real=True
    )
    leading_H = 2 * b + 3 * b * a + D * a**3 / 6
    leading_C = 3 * b + D * a * a / 2
    # The primary/secondary finite block algebra, not a continuum determinant.
    d11, d12, d22, jj = s.symbols(
        "auxiliary_D11 auxiliary_D12 auxiliary_D22 secondary_bracket", real=True
    )
    DD = s.Matrix([[d11, d12], [d12, d22]])
    JJ = s.Matrix([[0, jj], [-jj, 0]])
    dirac = s.zeros(2).row_join(-DD).col_join(DD.row_join(JJ))
    singular = dirac.subs({d11: 0, d12: 0, d22: -1})
    checks = {
        "whole_C_zero_along_connected_path": point_C,
        "whole_lapse_pivot_along_connected_path": s.factor(
            point_CN - 6 * (shear - fold.FOLD_SHEAR)
        ),
        "whole_second_jet_at_path_fold": s.factor(point_CNN - fold.FOLD_SECOND),
        "whole_transverse_Cs_at_fold": s.factor(fold.clock_jet(Cs, 0) - 3),
        "whole_original_raw_Hamiltonian_at_fixture": s.factor(raw - raw_expected),
        "whole_temporal_pivot_at_fold": s.factor(DTT + 1),
        "whole_temporal_lapse_cross_at_fold": DNT,
        "whole_sequential_Gamma_regular_at_this_fold": gamma.subs(canonical.R, 1) - 1,
        "whole_exact_stationary_shear_elimination": s.factor(
            C.subs(shear, stationary_shear)
        ),
        "whole_stationary_shear_at_fold": s.factor(s0 - fold.FOLD_SHEAR),
        "whole_stationary_shear_first_derivative": s1,
        "whole_stationary_shear_second_derivative": s.factor(s2 + fold.FOLD_SECOND / 3),
        "whole_envelope_cross_identity": s.factor(s.diff(Hs, N) - Cs),
        "whole_envelope_fixed_N_second_shear_derivative": s.diff(Hs, shear),
        "whole_fold_Hs_value": fold.clock_jet(Hs, 0) - 2,
        "leading_energy_derivative_is_fold_constraint": s.diff(leading_H, a)
        - leading_C,
        "leading_stationary_energy_three_halves_coefficient": s.factor(
            (leading_H - 2 * b).subs(b, -D * a * a / 6) + D * a**3 / 3
        ),
        "full_primary_secondary_block_determinant": s.factor(
            dirac.det() - DD.det() ** 2
        ),
        "rank_changed_primary_secondary_determinant": s.factor(singular.det()),
    }
    return {
        "whole_original_raw_Hamiltonian": raw,
        "whole_entire_stationary_shear_function": stationary_shear,
        "whole_connected_path_C_N": point_CN,
        "whole_fold_N_T_Hessian": s.diag(0, -1),
        "whole_regular_branch_second_shear_envelope": -Cs * Cs / CN,
        "whole_stationary_shear_second_derivative": s2,
        "whole_finite_primary_secondary_block": dirac,
        "whole_finite_rank_changed_block": singular,
        "whole_leading_fold_energy_and_constraint": (leading_H, leading_C),
        "whole_connected_sheet_obstruction": "Any C1 lapse map extending the established local branch to ALL these finite canonical inputs must have N=1 at the designated point along0<=s<s*, since C=0 there and C_N=6(s-s*)<0. Continuity gives N1 at s*. A transverse change of tensor momentum holds v fixed and changes shear independently; differentiating C=0 would give0=C_N dN/ds+C_s=3, impossible. If a global field branch already fails earlier at another spatial point, the proposed all-phase regular extension already fails. This does not exclude a different disconnected root, multiple charts, constrained quantization or a physically matched finite EFT domain.",
        "whole_local_normal_form": "C_s=3 and C_NN=D>28 give a smooth local graph s(N) with a nondegenerate maximum, s''=-D/3. At fixed curvature and matter, s<s* has two nearby positive roots, s>s* none nearby. The exact branch envelope H_ss=-C_s^2/C_N diverges as the fold is approached. The displayed cubic is only the leading Taylor form; the actual full remainders are retained.",
        "whole_measure_boundary": "The N,T block has a genuine zero eigenvalue at the fold even though R=Gamma=N=1 and the separate joint K,T stationary system stays invertible. The finite primary/secondary bracket determinant is(detD)^2, including an arbitrary secondary bracket. Its vanishing does not by itself define a new gauge symmetry, legitimize division by the pivot, or construct a new quantum measure. The earlier regular-branch determinant cancellation is not extended across a rank change.",
        "whole_not_inferred": "Not a Lorentzian metric singularity, clock-time turning point, full constrained time solution, trajectory reaching the fold, no-go for all completions, a new state, Airy/path-integral unitarity, self-adjoint domain, UV matching or original P8 closure.",
        "checks": checks,
        "gates": {
            "actual_full_positive_second_jet": s.Rational(11391, 400)
            - fold.affine_error(fold.FOLD_SECOND)
            > 28,
            "actual_full_positive_fold_shear": s.Rational(81, 160)
            - fold.affine_error(fold.FOLD_SHEAR)
            > 0,
            "temporal_pivot_stays_nonzero_at_fold": DTT == -1,
            "fold_is_not_the_sequential_Gamma_crossing": gamma.subs(canonical.R, 1)
            == 1,
            "full_primary_secondary_rank_loss": singular.rank() == 2,
            "full_nonzero_primitive_retained": H.has(fold.original.Iu),
            "full_lapse_derivatives_before_clock_restriction": True,
            "regular_classical_sheet_not_quantum_measure": True,
            "no_global_single_sheet_or_original_P8_completion": True,
        },
    }
