"""Exact complete finite-local factorization and its retained boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_curved_scalar_reference import geometry as g
from p8_vacuum_affine_curved_scalar_reference import local as original

t = g.eta
alpha = s.Function("local_alpha")(t)
gamma = s.Function("local_gamma")(t)
A0 = s.Matrix([[18, -6], [10, -4]])
H0 = s.Matrix([[24, -18], [-18, 12]])
J0 = s.Matrix([[0, 0], [10, -4]])
TV = s.Matrix([3, -1])


def B(vector):
    return g.L * vector.diff(t, 2) + g.Q * vector + g.C * vector.diff(t)


def Bstar(vector):
    return g.L.T * vector.diff(t, 2) + g.Q.T * vector - (g.C.T * vector).diff(t)


@cache
def data():
    old = original.data()
    YD, YG = s.Matrix([g.wd, g.cd]), s.Matrix([g.wg, g.cg])
    SD, _WD = B(YD)
    SG, _WG = B(YG)
    TD, TG = (TV.T * YD)[0], (TV.T * YG)[0]
    local = (
        alpha * old["mixed_rescaled_scalar_curvature"]
        + 6 * alpha * (TD * SG + TG * SD)
        + gamma * TD * TG
    )
    mass = g.bridge.MASS
    actual_alpha = -2 * g.U / 3 + 5 * mass**2 * g.a**2 / 3
    actual_gamma = (
        -2 * g.U**2 + 10 * mass**2 * g.a**2 * g.U + s.Rational(5, 2) * mass**4 * g.a**4
    )
    A = alpha * A0
    H = alpha * H0
    J = -s.diff(alpha, t) * J0
    V = gamma * TV * TV.T
    factored = (
        B(YD).T * A * YG
        + YD.T * A.T * B(YG)
        + YD.diff(t).T * H * YG.diff(t)
        + YD.diff(t).T * J * YG
        + YD.T * J.T * YG.diff(t)
        + YD.T * V * YG
    )[0]
    boundary = alpha * (
        s.diff(g.cd, t) * (10 * g.wg - 4 * g.cg)
        + (10 * g.wd - 4 * g.cd) * s.diff(g.cg, t)
    )
    euler = s.Matrix(
        [
            sum(
                (-1) ** j * s.diff(s.diff(local, s.diff(field, t, j)), t, j)
                for j in range(3)
            )
            for field in YD
        ]
    )
    operator = (
        Bstar(A * YG)
        + A.T * B(YG)
        - (H * YG.diff(t)).diff(t)
        - (J * YG).diff(t)
        + J.T * YG.diff(t)
        + V * YG
    )
    noJ = operator + (J * YG).diff(t) - J.T * YG.diff(t)
    checks = {
        "actual_complete_original_finite_local_remainder_bridge": s.simplify(
            local.subs({alpha: actual_alpha, gamma: actual_gamma})
            - old["complete_lower_order_local_remainder"]
        ),
        "exact_full_density_boundary_not_discarded": s.expand(
            local - factored - s.diff(boundary, t)
        ),
        "complete_variable_coefficient_Euler_Hessian_factorization": (
            euler - operator
        ).applyfunc(s.expand),
        "all_spatial_transfer_terms_in_curvature_factors": (euler - operator)
        .diff(g.k, 2)
        .applyfunc(s.expand),
        "actual_alpha_from_source_fixed_coefficients": s.simplify(
            actual_alpha
            - 12 * old["actual_fixed_coefficients"]["R_old_squared"] * g.U
            - old["actual_fixed_coefficients"]["R_old"] * g.a**2
        ),
        "actual_gamma_from_source_fixed_coefficients": s.simplify(
            actual_gamma
            - 36 * old["actual_fixed_coefficients"]["R_old_squared"] * g.U**2
            - 6 * old["actual_fixed_coefficients"]["R_old"] * g.a**2 * g.U
            - old["actual_fixed_coefficients"]["constant"] * g.a**4
        ),
        "A_matrix_Frobenius": sum(v * v for v in A0) - 476,
        "H_matrix_Frobenius": sum(v * v for v in H0) - 1368,
        "J_matrix_Frobenius": sum(v * v for v in J0) - 116,
        "trace_volume_rank_one_norm": TV.dot(TV) - 10,
    }
    return {
        "actual_scalar_coefficients": {"alpha": actual_alpha, "gamma": actual_gamma},
        "general_coefficient_factor_matrices": {"A": A, "H": H, "J": J, "V0": V},
        "exact_bilinear_boundary": boundary,
        "Euler_operator_identity": "Rloc=Bc* A+A^T Bc+D* H D+D* J+J^T D+V0, D*=-D, with every coefficient derivative retained.",
        "ordered_coordinate_conjugate": "Omega=Z Rloc Y=A Y+Z A^T+(Z D*)H(DY)+(Z D*)JY+Z J^T(DY)+Z V0 Y. Adjacent Bc/Bc* inverse pairs cancel; coefficients are not commuted.",
        "actual_bridge": "This is the complete S225 original finite local remainder, including second metric variation, background curvature, Einstein and volume terms. The flat factors already retain the original curvature-square finite terms.",
        "boundary_scope": "The density difference is exactly the displayed total derivative. The corresponding Euler-Hessian identity is an identity of causal distributions; all initial derivative atoms are retained rather than omitted.",
        "checks": checks,
        "gates": {
            "actual_remainder_not_zero": old["complete_lower_order_local_remainder"]
            != 0,
            "deleting_coefficient_derivative_J_changes_operator": (
                euler - noJ
            ).applyfunc(s.expand)
            != s.zeros(2, 1),
            "curvature_leg_matrix_norm_below22": 476 < 22**2,
            "kinetic_matrix_norm_below37": 1368 < 37**2,
            "coefficient_derivative_matrix_norm_below11": 116 < 11**2,
            "unchanged_original_finite_mass": mass == 1000,
        },
    }
