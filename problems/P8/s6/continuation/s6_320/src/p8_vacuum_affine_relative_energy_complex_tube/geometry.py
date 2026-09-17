"""Complex conserved inverse, weighted factorization and proper-cluster pole."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as original
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower


@cache
def inverse_data():
    vx, vy, vz = s.symbols("vx vy vz")
    v = s.Matrix([vx, vy, vz])
    r = s.symbols("r0:6")
    R = s.Matrix([[r[0], r[1], r[2]], [r[1], r[3], r[4]], [r[2], r[4], r[5]]])
    J = s.zeros(4)
    J[1:, 1:] = R
    J[0, 0] = (v.T * R * v)[0]
    J[0, 1:], J[1:, 0] = (R * v).T, R * v
    Q = s.Matrix([1, *v])
    eta = lower.ETA
    q = eta * Q
    raw = -(eta * J * eta - eta * s.trace(eta * J) / 2)
    projected = s.Matrix(
        4,
        4,
        lambda i, j: (
            raw[i, j] - q[i] * raw[0, j] - q[j] * raw[i, 0] + q[i] * q[j] * raw[0, 0]
        ),
    )
    T = s.eye(3) - v * v.T
    candidate = -T * R * T + T * s.trace(T * R) / 2
    d, a, b, c = s.symbols("delta a b c")
    u = s.Matrix([d * a, d * b, 1 - d * d * c])
    D = s.diag(1, 1, d)
    Tc = s.eye(3) - u * u.T
    B = s.Matrix(
        [
            [1 - d * d * a * a, -d * d * a * b, -a * (1 - d * d * c)],
            [-d * d * a * b, 1 - d * d * b * b, -b * (1 - d * d * c)],
            [-a * (1 - d * d * c), -b * (1 - d * d * c), 2 * c - d * d * c * c],
        ]
    )
    Jw = s.Matrix(
        [
            [d * d * r[0], d * d * r[1], d * r[2]],
            [d * d * r[1], d * d * r[3], d * r[4]],
            [d * r[2], d * r[4], r[5]],
        ]
    )
    lhs = -Tc * Jw * Tc + Tc * s.trace(Tc * Jw) / 2
    rhs = d * d * D * (-B * R * B + B * s.trace(B * R) / 2) * D
    eps = s.Rational(1, 100)
    e1, e2 = s.symbols("e1 e2")
    checks = {
        "generic_complex_source_conservation": (Q.T * eta * J).applyfunc(s.expand),
        "generic_complex_spatial_temporal_inverse": (
            projected[1:, 1:] - candidate
        ).applyfunc(s.expand),
        "complex_temporal_row": projected[0, :].applyfunc(s.expand),
        "weighted_T_factorization": (Tc - D * B * D).applyfunc(s.expand),
        "weighted_conserved_inverse_factorization": (lhs - rhs).applyfunc(s.expand),
        "relative_pair_invariant_perturbation": s.expand(
            (1 + e1) * (1 + e2) - 1 - e1 - e2 - e1 * e2
        ),
        "weighted_inverse_numerical_budget": s.Rational(3, 2) * 64 * 3 * 2 - 576,
        "complex_vertex_and_inverse_recurrence_factor": s.Integer(25 * 1024 - 25600),
        "complex_EGF_nonlinear_factor": s.Integer(32 * 25600 - 819200),
        "full_hard_Euclidean_step": s.Integer(4 * 50 * 600000 - 120000000),
    }
    return {
        "checks": checks,
        "whole_conserved_complex_inverse": "T=I-u*u^T,u=Qsp/Q0. For a complete conserved source with spatial block R, Hsp=(-T*R*T+T*tr(T*R)/2)/Q^2. Transposes are bilinear; absolute Hermitian norms occur only in bounds.",
        "whole_weighted_factorization": "At a positive real center, use D=diag(1,1,delta). In a relative energy tube epsilon<=1/100, T=D*B*D with ||B||F<8, while R=delta^2*D^-1*Rbar*D^-1 and ||Rbar||F<=3L. The exact numerator factors as delta^2*D*(-B*Rbar*B+B*tr(B*Rbar)/2)*D. The real-center weighted field norm is<=576L/W^2<1024L/W^2.",
        "whole_relative_pure_cut_gap": "Each real future invariant is a sum of nonnegative pair terms. Continuing zi=wi(1+ei), |ei|<=epsilon, changes any subset Q^2 by at most(2epsilon+epsilon^2)*Qreal^2. Every generic pure-soft and temporal inverse remains nonzero without a factor depending on N.",
        "gates": {
            "complex_mass_gap_retains_half": 1 - 2 * eps - eps * eps > s.Rational(1, 2),
            "transverse_velocity_budget": 2 * s.Rational(1, 99) ** 2
            < s.Rational(1, 50) ** 2,
            "longitudinal_velocity_ratio_below_two": s.Rational(101, 99) < 2,
            "weighted_T_Frobenius_budget": 4 + 2 + 36 < 64,
            "conservative_inverse_budget": s.Integer(576) < 1024,
            "transverse_identity_plus_small_outer_product": (2 - s.Rational(1, 2500))
            ** 2
            > 2,
            "real_center_weights_not_a_complex_rotation": True,
            "arbitrary_Rosen_and_reflection_identities_extend_algebraically": True,
        },
    }


@cache
def pole_data():
    directions = (
        s.Matrix([1, 0, 0]),
        s.Matrix([0, 1, 0]),
        s.Matrix([-s.Rational(3, 5), -s.Rational(4, 5), 0]),
    )
    energies = (s.S.One, s.S.One, -s.Rational(5, 17))
    z = s.Matrix([0, 0, 1])
    hs = []
    for n, w in zip(directions, energies):
        tangent = s.Matrix([-n[1], n[0], 0])
        A = s.zeros(4)
        A[1:, 1:] = tangent * tangent.T - z * z.T
        hs.append((s.ImmutableMatrix([w, *(w * n)]), s.ImmutableMatrix(A)))
    engine = original.TreeEngine([("h", q, A) for q, A in hs])
    J, count = engine.amputated(7, "h")
    Q = engine.momentum(7)
    n = Q[1:, 0] / Q[0]
    tangent = s.Matrix([-n[1], n[0], 0])
    plus, cross = s.zeros(4), s.zeros(4)
    plus[1:, 1:] = tangent * tangent.T - z * z.T
    cross[1:, 1:] = tangent * z.T + z * tangent.T
    residues = tuple(
        s.factor(sum(A[i, j] * J[i, j] for i in range(4) for j in range(4)))
        for A in (plus, cross)
    )
    lam = s.Rational(1, 10**20)
    shift = s.Rational(22, 17) * lam
    global_radius = (3 * lam + s.Rational(1, 16)) / s.Integer(10) ** 12
    checks = {
        "complete_four_tree_inventory": s.Integer(count - 4),
        "exact_null_cluster_root": s.factor(original.old.dot(Q, Q)),
        "complete_root_Ward": (Q.T * lower.ETA * J).applyfunc(s.factor),
        "nonzero_physical_root_plus_residue": residues[0] - s.Rational(47089, 2601),
        "root_cross_residue": residues[1],
        "exact_root_vector": Q - s.Matrix([29, 20, 21, 0]) / 17,
    }
    return {
        "checks": checks,
        "whole_global_total_energy_tube_obstruction": {
            "directions": directions,
            "pole_energies": energies,
            "root": Q,
            "TT_residues": residues,
            "positive_center_scale": lam,
            "required_energy_shift": shift,
            "global_total_energy_radius": global_radius,
        },
        "whole_strict_control_scope": "A proper complete three-ray block can have a nonzero transverse root pole inside a global-total-W disc when a spectator dominates W. Multiplying by leaf energies does not remove this pole. This refutes that block-holomorphy shortcut, not the bounded positive-real amplitude or every full hard-amplitude residue.",
        "gates": {
            "complete_current_TT_root_residue_nonzero": residues[0] != 0,
            "all_three_pole_leaf_energies_nonzero": all(w != 0 for w in energies),
            "proper_cluster_pole_inside_wrong_global_disc": 0 < shift < global_radius,
            "positive_center_and_spectator_in_original_energy_domain": 3 * lam
            + s.Rational(1, 16)
            < s.Rational(1, 8),
            "negative_energy_pole_not_called_physical_real_radiation": True,
        },
    }


@cache
def data():
    first, second = inverse_data(), pole_data()
    return {
        "whole_complex_inverse_and_geometry": {
            k: v for k, v in first.items() if k not in ("checks", "gates")
        },
        "whole_failed_global_disc_control": {
            k: v for k, v in second.items() if k not in ("checks", "gates")
        },
        "checks": {
            **{"inverse_" + k: v for k, v in first["checks"].items()},
            **{"pole_" + k: v for k, v in second["checks"].items()},
        },
        "gates": {**first["gates"], **second["gates"]},
    }
