"""New covariant isotropic trace mass: exact full quotient update and reduction."""

from functools import cache

import sympy as s
from p8_affine import connection as old
from p8_affine_retuned import geometry


@cache
def matrices():
    d = geometry.update()
    p = old.P
    eta = geometry.ETA
    N, M, W = d["N"], d["M_old"], d["W_old"]
    D = d["D_old"]
    Xi = eta - D.inv()
    lift = geometry.clean(W * D.inv())
    updated = geometry.clean(M + N.T * Xi * N)
    newW = lift * eta
    projector = geometry.clean(s.eye(60) - lift * N)
    return {
        "p": p,
        "eta": eta,
        "N": N,
        "M_old": M,
        "D_old": D,
        "Xi": Xi,
        "lift": lift,
        "M_new": updated,
        "W_new": newW,
        "projector": projector,
        "checks": {
            "updated_all60_Euler_lift": geometry.clean(updated * newW - N.T),
            "retained_response_exact_physical_Lorentz": geometry.clean(N * newW - eta),
            "full_trace_right_inverse": geometry.clean(N * lift - s.eye(4)),
            "complement_trace_zero": geometry.clean(N * projector),
            "full_complement_cross_term": geometry.clean(projector.T * updated * lift),
            "retained_mass_exactly_Lorentz": geometry.clean(
                lift.T * updated * lift - eta
            ),
            "full_projective_columns": geometry.clean(
                d["N_full"] * old.quadratic()["gauge"]
            ),
        },
    }


@cache
def data():
    m = matrices()
    p = m["p"]
    tt = 3 * (2 * p**3 - 1) / p
    ss = (8 * p + 5) / (8 * p * p)
    xi_t = 1 - 1 / tt
    xi_s = -1 - 1 / ss
    anisotropy = s.factor(xi_t + xi_s)
    R = s.symbols("positive_R", positive=True)
    gamma = 1 - 3 * (R - 1) ** 2 / (2 * R)
    X = s.symbols("X", real=True)
    a = -xi_s
    b_numerator = xi_t - a
    lower = s.Rational(1, 2)
    upper = s.Rational(6, 5)
    return {
        "candidate": "CD-REG-AFFINE-ISO",
        "new_mass_update_at_timelike_frame": s.diag(xi_t, xi_s, xi_s, xi_s),
        "covariant_mass_update_metric_coefficient": a,
        "covariant_mass_update_u_u_coefficient": b_numerator / X,
        "mass_update_action": "Add +1/2*(T-Tstar)_mu*Xi^{mu nu}*(T-Tstar)_nu and -zeta*F(T-B du)^2/4 to the source-pinned S6.109 affine action, all with overall kappa. Source centering includes its linear and constant terms.",
        "exact_reduced_action": "S_CD_target+S_M1+kappa*int sqrt(-g)[-zeta F(W)^2/4+(W-S)_mu*(W-S)^mu/2]",
        "algebraic_reduction_scope": "The 56 trace-complement directions remain algebraic; the four projective directions remain gauge. Eliminating them is exact, with no omitted derivative expansion. W is retained and is not integrated out.",
        "unchanged_stationary_connection_when_zeta_zero": True,
        "updated_quotient_determinant_ratio": s.factor(-tt * ss**3),
        "normalized_temporal_schur_factor": gamma,
        "temporal_schur_lower": s.Rational(1, 4),
        "same_domain": (lower, upper),
        "checks": {
            "full_update_determinant_lemma": s.factor(
                (s.eye(4) + m["Xi"] * m["D_old"]).det() + tt * ss**3
            ),
            "mass_update_clock_equals_old_11_over_9": s.factor(
                xi_t.subs(p, s.Rational(1, 2)) - s.Rational(11, 9)
            ),
            "spatial_mass_update_clock_equals_old": s.factor(
                xi_s.subs(p, s.Rational(1, 2)) + s.Rational(11, 9)
            ),
            "anisotropy_has_p_minus_half_factor": s.factor(
                anisotropy
                + p
                * (2 * p - 1)
                * (24 * p**3 + 12 * p**2 + 6 * p - 5)
                / (3 * (8 * p + 5) * (2 * p**3 - 1))
            ),
            "regular_anisotropy_conversion": s.factor(
                (2 * p - 1) - (4 * p * p - 1) / (2 * p + 1)
            ),
            "positive_temporal_schur_margin_identity": s.factor(
                gamma
                - s.Rational(1, 4)
                - 3 * (R - s.Rational(1, 2)) * (2 - R) / (2 * R)
            ),
        },
        "not_inferred": "A full physical constraint count away from the clock/vacuum, a stable cone on every background, a physical cutoff, interacting quantum matching, or UV admissibility.",
    }
