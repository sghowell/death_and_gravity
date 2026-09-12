"""Normalized scalar channels and the three still-independent ordered kernels."""

from functools import cache
from itertools import permutations, product

import sympy as s

TRACE = s.eye(3) / s.sqrt(3)


def scalar_axis(e):
    return (3 * e * e.T - s.eye(3)) / s.sqrt(6)


def scalar_components(Q, e):
    axis = scalar_axis(e)
    return s.trace(TRACE * Q), s.trace(axis * Q)


def rotations():
    out = []
    for perm in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            R = s.zeros(3)
            for j in range(3):
                R[j, perm[j]] = signs[j]
            if R.det() == 1:
                out.append(R)
    return out


@cache
def data():
    q0, q1, q2, q3, q4, q5 = s.symbols("q0:6", real=True)
    Q = s.Matrix([[q0, q3, q4], [q3, q1, q5], [q4, q5, q2]])
    e = s.Matrix([0, 0, 1])
    axis = scalar_axis(e)
    R = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    avg = sum((R**j * Q * (R.T) ** j for j in range(4)), s.zeros(3)) / 4
    u, v = scalar_components(Q, e)
    rots = rotations()
    origin = sum((r * Q * r.T for r in rots), s.zeros(3)) / 24
    rtt, rt0, r0t, r00 = s.symbols("Rtt Rt0 R0t R00")
    D = s.Matrix([[2, 1, 3], [1, 4, -1], [3, -1, 5]])
    d = s.Matrix(scalar_components(D, e))
    g = s.Matrix([u, v])
    K = s.Matrix([[rtt, rt0], [r0t, r00]])
    Kknown = s.Matrix([[0, 0], [0, r00]])
    missing = d[0] * rtt * g[0] + d[0] * rt0 * g[1] + d[1] * r0t * g[0]
    checks = {
        "unit_normalized_trace": s.trace(TRACE * TRACE) - 1,
        "unit_normalized_tracefree_scalar": s.trace(axis * axis) - 1,
        "trace_and_scalar_orthogonal": s.trace(TRACE * axis),
        "axial_rotation_average_retains_only_two_scalars": avg - u * TRACE - v * axis,
        "full_origin_rotational_average_is_pure_trace": origin
        - s.trace(Q) * s.eye(3) / 3,
        "proper_origin_rotation_group_order": len(rots) - 24,
        "ordered_scalar_block_has_exactly_three_missing_entries": (
            d.T * (K - Kknown) * g
        )[0]
        - missing,
        "scalar_tracefree_axis_even_in_momentum": scalar_axis(-e) - axis,
    }
    return {
        "normalized_basis": "For P!=0, e=P/|P|, Etr=I/sqrt3 and E0=(3ee^T-I)/sqrt6. Both are Frobenius unit and orthogonal. Write Q=gamma+u Etr with gamma tracefree, and v=E0:gamma.",
        "rotation_reduction": "The unchanged isotropic three-polarization state and covariant prescription imply simultaneous spatial-rotation covariance. The trace is rotation-invariant, so its cross pairing with a tracefree tensor depends only on its axial average v E0. Vector and transverse tensor directions cannot couple to trace. A four-element axial subgroup already gives the exact rank-two scalar projection.",
        "three_ordered_kernels": [
            "R_trace_trace",
            "R_trace_scalar_tracefree",
            "R_scalar_tracefree_trace",
        ],
        "known_block": "S213 supplies the complete gamma-to-gamma matched spatial response, including the scalar-tracefree diagonal. The three displayed ordered scalar entries are additional; the two off-diagonal retarded entries are not identified.",
        "origin": "At P0 full SO3 covariance forces trace/tracefree cross entries to vanish, but does not determine the homogeneous trace/trace current. The bounded scalar projectors need not be analytic at P0; no Borel-state external-momentum analyticity is asserted.",
        "conditional_completion": "Construct, match and bound each missing scalar kernel in the same original preparation and prescription. Only then may the Ward formula reconstruct the full conditional Gaussian metric response. Symmetry alone supplies none of their numerical values or estimates.",
        "checks": checks,
        "gates": {
            "complete_tracefree_sector_unchanged": True,
            "three_ordered_scalar_entries_not_two": rt0 != r0t,
            "retarded_kernel_not_single_branch_symmetric_Hessian": True,
            "zero_transfer_cross_vanishes_not_trace_anchor": True,
            "unit_Frobenius_projectors_without_helicity_basis_factor": True,
            "no_unproved_scalar_UV_matching_or_norm": True,
        },
    }
