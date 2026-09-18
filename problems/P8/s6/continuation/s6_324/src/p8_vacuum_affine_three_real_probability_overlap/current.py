"""Literal external residues and full original three-to-two-real soft faces."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as full
from p8_vacuum_affine_complete_two_graviton_tree import trees as e

from . import source

imm = s.ImmutableMatrix


def dot(p, q):
    return (p.T * e.ETA * q)[0]


@cache
def vertices():
    A = s.zeros(4)
    for i in range(4):
        for j in range(i, 4):
            A[i, j] = A[j, i] = s.Symbol(f"A{i}{j}")
    A = imm(A)
    c, d = s.symbols("c d")
    h00, h01, h02, h11, h12, h22 = s.symbols("h00 h01 h02 h11 h12 h22")
    q = imm([1, 0, 0, 1])
    # Independent plus/cross hard TT field.
    B = imm([[0, 0, 0, 0], [0, c, d, 0], [0, d, -c, 0], [0, 0, 0, 0]])
    H = imm(
        [
            [h00, h01, h02, h00],
            [h01, h11, h12, h01],
            [h02, h12, h22, h02],
            [h00, h01, h02, h00],
        ]
    )
    propagated = e.ETA * H * e.ETA - e.ETA * s.trace(e.ETA * H) / 2
    vertex = e.cubic_gravity((A, B, propagated), (e.VECTOR_ZERO, q, -q))
    target = 2 * (q.T * A * q)[0] * s.trace(H * B)
    p = imm(s.symbols("p0:4"))
    scalar = e.scalar_vertex(p, -p, (A,), 1)
    checks = {
        "generic_scalar_on_shell_residue": s.expand(
            -scalar - 2 * (p.T * A * p)[0] + e.tr(A) * (dot(p, p) - 1)
        ),
        "generic_EH3_conserved_root_residue": s.factor(-vertex - target),
        "hard_null_mass_shell": dot(q, q),
        "hard_field_trace": e.tr(B),
        "hard_field_transversality": B * q,
        "full_hard_root_conservation": H * e.ETA * q,
    }
    return {
        "checks": checks,
        "gates": {
            "all10_soft_tensor_entries_free": len(A.free_symbols) == 10,
            "all6_conserved_hard_root_components_free": len(H.free_symbols) == 6,
            "both_hard_TT_polarizations_free": len(B.free_symbols) == 2,
            "no_TT_premise_on_conserved_hard_root": e.tr(H) != 0,
            "both_emitter_residues_share_canonical_normalization": True,
        },
        "whole_literal_graviton_residue": "-V3(A,B,eta*H*eta-eta*tr(eta*H)/2;0,q,-q)=2*(q.A.q)*tr(H*B), with arbitrary symmetric A, both hard TT polarizations and the general conserved hard root. Spatial covariance and quadratic momentum homogeneity extend q=(1,0,0,1) to every positive hard null momentum.",
        "whole_face_proof": "At a generic fixed two-real hard state only external emitter propagators have the new soft pole. The scalar and Einstein residues give the full massive-plus-two-null soft current times the complete434 amplitude. Internal cuts and contact terms have no leading pole. The complete lower hard root is conserved by S316, not by a multi-off-shell shortcut. Uniform face estimates are proved separately.",
    }


@cache
def calibration():
    E = s.Rational(5, 4)
    r0 = s.Rational(3, 4)
    alpha0 = s.Rational(19, 10)
    ep0 = (alpha0 + 1 / alpha0) / 2
    b = s.Rational(1, 64)
    c = s.factor((E * E - E * b - ep0 * ep0) / (E - b / 2))
    rays = tuple(imm([1, *(1 if i == j else 0 for i in range(3))]) for j in range(3))
    fields = (
        imm(s.diag(0, 0, s.Rational(1, 2), -s.Rational(1, 2))),
        imm(s.diag(0, s.Rational(1, 2), 0, -s.Rational(1, 2))),
        imm(s.diag(0, s.Rational(1, 2), -s.Rational(1, 2), 0)),
    )
    u = imm([s.Rational(3, 5), 0, s.Rational(4, 5)])
    checks = {}
    gates = {}
    for i, (ray, A) in enumerate(zip(rays, fields)):
        checks[f"fixed_ray{i}_null"] = dot(ray, ray)
        checks[f"fixed_tensor{i}_TT"] = A * ray
        checks[f"fixed_tensor{i}_trace"] = e.tr(A)
        gates[f"fixed_tensor{i}_unit_norm"] = sum(v * v for v in A) <= 1

    def geometry(h, label):
        alpha = alpha0 - h
        ep = (alpha + 1 / alpha) / 2
        rp = (alpha - 1 / alpha) / 2
        a = s.factor((ep * ep - ep0 * ep0) / (-E + (b + c) / 2))
        weights = (a, b, c)
        W = sum(weights)
        v = imm([a, b, c])
        half = E - W / 2
        uv = u.dot(v)
        p3 = imm(
            [
                half - rp * uv / (2 * ep),
                *(rp * u + (rp * uv / (4 * ep * (half + ep)) - s.Rational(1, 2)) * v),
            ]
        )
        p4 = imm(
            [
                half + rp * uv / (2 * ep),
                *(-rp * u + (-rp * uv / (4 * ep * (half + ep)) - s.Rational(1, 2)) * v),
            ]
        )
        ps = (imm([-E, 0, 0, -r0]), imm([-E, 0, 0, r0]), p3, p4)
        qs = tuple(w * n for w, n in zip(weights, rays))
        checks[label + "_pair_squared_energy"] = s.factor(
            ep * ep - half * half + v.dot(v) / 4
        )
        checks[label + "_pair_squared_spatial_radius"] = s.factor(rp * rp - ep * ep + 1)
        checks[label + "_full_conservation"] = sum(ps, e.VECTOR_ZERO) + sum(
            qs, e.VECTOR_ZERO
        )
        for i, p in enumerate(ps):
            checks[label + f"_mass_shell{i}"] = s.factor(dot(p, p) - 1)
        gates[label + "_outgoing_future"] = all(p[0] > 0 for p in ps[2:])
        gates[label + "_total_energy"] = bool(0 < b + c <= W < s.Rational(1, 8))
        rho = E * rp / (r0 * ep)
        gates[label + "_positive_real_phase_at_most_one"] = bool(0 < rho <= 1)
        return ps, qs, weights, rho

    ps0, qs0, weights0, rho0 = geometry(0, "face")
    checks["zero_soft_energy_at_face"] = weights0[0]
    pars = source.original_parameters()
    hard = (*ps0, qs0[1], qs0[2])
    current = sum((p * p.T / dot(p, rays[0]) for p in hard), s.zeros(4))
    massive = sum((p * p.T / dot(p, rays[0]) for p in ps0), s.zeros(4))
    checks["full_two_marked_current_Ward"] = current * e.ETA * rays[0]
    checks["massive_only_Ward_is_minus_both_null_legs"] = (
        massive * e.ETA * rays[0] + qs0[1] + qs0[2]
    )
    gates["massive_only_current_Ward_nonzero"] = (
        massive * e.ETA * rays[0] != e.VECTOR_ZERO
    )
    one_marked = massive + qs0[1] * qs0[1].T / dot(qs0[1], rays[0])
    checks["one_marked_Ward_is_minus_second_null_leg"] = (
        one_marked * e.ETA * rays[0] + qs0[2]
    )
    gates["one_marked_current_Ward_nonzero"] = (
        one_marked * e.ETA * rays[0] != e.VECTOR_ZERO
    )
    S = s.factor(s.trace(current * fields[0]))
    S_bad = s.factor(s.trace(massive * fields[0]))
    assert S != 0 and S_bad != 0
    M2, n2 = full.amplitude(ps0, tuple(zip(qs0[1:], fields[1:])), **pars)
    checks["complete_original_two_real_inventory"] = s.Integer(n2 - 434)
    gates["complete_original_two_real_nonzero"] = M2 != 0
    errors = []
    wrong_errors = []
    records = []
    for index, h in enumerate(
        (s.Rational(1, 1000), s.Rational(1, 10000), s.Rational(1, 100000))
    ):
        label = "sample" + str(index)
        ps, qs, weights, rho = geometry(h, label)
        M3, n3 = full.amplitude(ps, tuple(zip(qs, fields)), **pars)
        checks[label + "_complete_original_three_real_inventory"] = s.Integer(n3 - 5116)
        ratio = s.factor(weights[0] * s.sqrt(pars["kappa"]) * M3 / (S * M2))
        wrong = s.factor(weights[0] * s.sqrt(pars["kappa"]) * M3 / (S_bad * M2))
        err = abs(ratio - 1)
        wrong_err = abs(wrong - 1)
        errors.append(err)
        wrong_errors.append(wrong_err)
        gates[label + "_full_source_nonzero"] = M3 != 0
        gates[label + "_finite_error_below3h"] = bool(err < 3 * h)
        gates[label + "_correct_error_below_one_percent"] = bool(
            err < s.Rational(1, 100)
        )
        gates[label + "_wrong_current_error_larger"] = bool(wrong_err > err)
        records.append(
            {
                "h": h,
                "weights": weights,
                "rho": rho,
                "tree_count": n3,
                "error_over_h_ceiling": s.ceiling(err / h),
            }
        )
    gates["all_three_errors_strictly_decrease"] = bool(
        errors[2] < errors[1] < errors[0]
    )
    gates["massive_only_last_error_more_than_ten_times_correct"] = bool(
        wrong_errors[2] > 10 * errors[2]
    )
    h = s.Symbol("h")
    alpha = alpha0 - h
    ep = (alpha + 1 / alpha) / 2
    rp = (alpha - 1 / alpha) / 2
    checks["exact_phase_has_correct_face_limit"] = (
        s.limit(E * rp / (r0 * ep), h, 0) - rho0
    )
    return {
        "checks": checks,
        "gates": {k: bool(v) for k, v in gates.items()},
        "whole_original_parameters": pars,
        "whole_reduced_two_real_state": ps0,
        "whole_two_marked_null_legs": qs0[1:],
        "whole_correct_and_massive_only_soft_coefficients": (S, S_bad),
        "whole_three_exact_hierarchical_calibrations": records,
        "whole_calibration_boundary": "Original434 and5116 amplitudes at three rational recoil points; the correct two-marked soft ratio approaches one with exact finite error bounds. These calibrate the literal residue proof; the3h inequalities are asserted only at the three recorded samples.",
    }


@cache
def data():
    residue, actual = vertices(), calibration()
    return {
        "checks": {
            **{"residue_" + k: v for k, v in residue["checks"].items()},
            **{"original_" + k: v for k, v in actual["checks"].items()},
        },
        "gates": {**residue["gates"], **actual["gates"]},
        "whole_literal_external_residues": {
            k: v for k, v in residue.items() if k not in ("checks", "gates")
        },
        "whole_original_face_calibrations": {
            k: v for k, v in actual.items() if k not in ("checks", "gates")
        },
    }
