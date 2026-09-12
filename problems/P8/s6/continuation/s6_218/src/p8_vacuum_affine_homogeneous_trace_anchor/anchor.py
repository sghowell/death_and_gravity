"""Actual zero-transfer trace anchor, ordered cross zeros and Hilbert lift."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_covariant_current import current as shear
from p8_vacuum_affine_current_response import tail

from . import comparison, local

CURRENT = (s.Integer(10) ** 78, s.Integer(10) ** 95, 2 * s.Integer(10) ** 111)
FULL_HOMOGENEOUS = s.Integer(10) ** 95


def symmetric(prefix, n):
    entries = iter(s.symbols(prefix + ":" + str(n * (n + 1) // 2), real=True))
    out = s.zeros(n)
    for i in range(n):
        for j in range(i, n):
            out[i, j] = out[j, i] = next(entries)
    return out


@cache
def data():
    e, x = s.symbols("epsilon x", real=True)
    M0, M1, M2 = (symmetric(p, 4) for p in ("M0_", "M1_", "M2_"))
    C0, C1 = (symmetric(p, 4) for p in ("C0_", "C1_"))
    J = s.zeros(4)
    J[:2, 2:] = s.eye(2)
    J[2:, :2] = -s.eye(2)
    M, C = M0 + e * M1 + e**2 * M2 / 2, C0 + e * C1
    rhs = J * M * C + C * (J * M).T
    wanted = J * M0 * C1 + C1 * (J * M0).T + J * M1 * C0 + C0 * (J * M1).T
    readout = -s.trace((M1 + e * M2) * C) / 2
    Q, D = symmetric("Q", 3), symmetric("D", 3)
    Qt, Dt = Q - s.trace(Q) * s.eye(3) / 3, D - s.trace(D) * s.eye(3) / 3
    Q0, D0 = s.trace(Q) * s.eye(3) / 3, s.trace(D) * s.eye(3) / 3
    alpha, beta = s.symbols("alpha beta", real=True)
    equivariant = lambda q: alpha * q + beta * s.trace(q) * s.eye(3)
    total = tuple(tail.COMPARISON[j] + local.LOCAL for j in range(3))
    checks = {
        "full_finite_mode_tangent_with_both_source_terms": (
            s.diff(rhs, e).subs(e, 0) - wanted
        ).applyfunc(s.expand),
        "complete_trace_metric_second_vertex_contact": s.expand(
            s.diff(readout, e).subs(e, 0) + s.trace(M2 * C0 + M1 * C1) / 2
        ),
        "ordered_trace_source_tracefree_detector_anchor_zero": s.expand(
            s.trace(Dt * equivariant(Q0))
        ),
        "ordered_tracefree_source_trace_detector_anchor_zero": s.expand(
            s.trace(D0 * equivariant(Qt))
        ),
        "full_Frobenius_trace_split": s.expand(
            s.trace(Q.T * Q) - s.trace(Qt.T * Qt) - s.trace(Q0.T * Q0)
        ),
        "unit_trace_two_leg_conversion": (s.sqrt(3) / 2) ** 2 - s.Rational(3, 4),
        "same_shear_anchor_first_derivative_bound": shear.CURRENT[1] - FULL_HOMOGENEOUS,
        "complete_second_derivative_Taylor_factor": CURRENT[2]
        * s.integrate(1 - x, (x, 0, 1))
        - s.Integer(10) ** 111,
        "zero_amplitude_Taylor_edge": (e * e * s.Integer(10) ** 111).subs(e, 0),
        "canonical_unit_trace_anchor_factor": s.Rational(3, 4)
        * 4
        * CURRENT[1]
        / modes.KAPPA
        - 3 * s.Rational(1, 10) ** 705,
        "canonical_full_homogeneous_anchor_factor": 4 * FULL_HOMOGENEOUS / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 705,
    }
    return {
        "complete_trace_current_C2_bounds": CURRENT,
        "unrounded_trace_current_sums": total,
        "actual_trace_anchor": "At P0 the full S214/S216 constrained spatial vertex and second vertex are exactly the first two derivatives of this trace Hamiltonian. The covariance tangent has both J M_G C and C M_G J^t terms, with zero initial tangent from the common unchanged germ. Finite-mode ODE uniqueness, the proved full comparison bounds and original local matching identify the actual prescribed trace/trace response with the first derivative of this complete homogeneous current.",
        "finite_amplitude_boundary": "The displayed C2 family and epsilon^2*1e111 Taylor remainder concern the Gaussian determinant component. They justify the actual reference metric Hessian because S176's source and first source variation vanish there. They are not a uniform finite-amplitude remainder for the full affine sourced parent.",
        "ordered_cross_anchors": "At zero external momentum the unchanged state, regulator and prescribed covariant local terms are rotation invariant. The causal spatial tensor map commutes with O(3), so it has the form alpha Q+beta tr(Q)I with time-operator coefficients. Both ordered trace/tracefree cross anchors vanish separately. No equality of their nonzero-transfer retarded kernels is inferred.",
        "unit_trace_normalization": "For T0=I/sqrt(3), Q=T0 G0 equals(2/3)phi I with phi=sqrt(3)G0/2. The same conversion on the detector gives3/4 times the scalar phi Hessian. This is a nonzero finite rescaling, not deletion of the trace.",
        "full_homogeneous_Hilbert_lift": "Combine the unchanged tracefree anchor with the new normalized trace anchor, whose cross blocks vanish. The Frobenius source splits orthogonally. The initial germ gives sup_t|G^(j)|<=||G^(j+1)||L2 for j<=12. Apply the real bound separately to real/imaginary parts and square-sum orthogonal outputs. The full identity-multiplier homogeneous anchor obeys1e95||D||L2 Z130[G], Z130^2=sum_r0..13||partial_t^r G||L2^2. This is the explicit P0 kernel's spatial Hilbert lift, not evaluation of an arbitrary weak multiplier at a point or a full inhomogeneous current.",
        "canonical_scope": "Both external metric factors give3e-705 for the unit-trace anchor and4e-705 for the full homogeneous lift. These remain derivative-losing Gaussian weak bounds, not a reduced inverse, stability or nonlinear quantum background.",
        "checks": checks,
        "gates": {
            "full_comparison_hypotheses_checked": all(
                bool(v) for v in comparison.data()["gates"].values()
            ),
            "original_local_matching_and_volume_bounds_checked": all(
                bool(v) for v in local.data()["gates"].values()
            ),
            "all_complete_trace_current_C2_displays": all(
                total[j] < CURRENT[j] for j in range(3)
            ),
            "normalized_trace_bound_fits_full_homogeneous_display": s.Rational(3, 4)
            * CURRENT[1]
            < FULL_HOMOGENEOUS,
            "source_germ_not_independent_quantum_metric_reset": True,
            "both_ordered_cross_anchors_zero_not_nonzero_transfer_equality": True,
            "full_scalar_state_time_contact_and_transfer_remainders_remain": True,
            "no_full_sourced_parent_finite_amplitude_claim": True,
        },
    }
