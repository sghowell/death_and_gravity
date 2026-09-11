"""Bounded TT projection and nonzero compactly supported detector construction."""

from functools import cache

import sympy as s


@cache
def data():
    k = s.Matrix(s.symbols("k0:3", real=True))
    k2 = k.dot(k)
    Q = k2 * s.eye(3) - k * k.T
    A = s.Matrix(
        3, 3, lambda i, j: s.Symbol("A" + str(min(i, j)) + str(max(i, j)), real=True)
    )
    C = Q * A * Q - Q * s.trace(Q * A) / 2
    n = s.Matrix([0, 0, 1])
    P = s.eye(3) - n * n.T
    TT = lambda M: P * (M + M.T) * P / 2 - P * s.trace(P * M) / 2
    q = s.Symbol("q", real=True)
    control = C.subs({k[0]: 0, k[1]: 0, k[2]: q}).subs(
        {A[i, j]: int(i == 0 and j == 0) for i in range(3) for j in range(i, 3)}
    )
    expected = s.diag(q**4 / 2, -(q**4) / 2, 0)
    checks = {
        "polynomial_transverse_projector": s.simplify(k.T * Q),
        "polynomial_projector_square": s.simplify(Q * Q - k2 * Q),
        "polynomial_projector_trace": s.trace(Q) - 2 * k2,
        "full_local_compact_tensor_transverse": s.simplify(k.T * C),
        "full_local_compact_tensor_trace": s.expand(s.trace(C)),
        "full_local_compact_tensor_symmetric": s.simplify(C - C.T),
        "nonzero_compact_TT_symbol": s.simplify(control - expected),
        "orthogonal_TT_projection_idempotent": s.simplify(TT(TT(A)) - TT(A)),
        "orthogonal_TT_projection_trace": s.simplify(s.trace(TT(A))),
        "orthogonal_TT_projection_transverse": s.simplify(n.T * TT(A)),
    }
    return {
        "bounded_projector": "Pi_TT(A)=P A_sym P-(1/2)P tr(PA), P=I-n n^T, n=k/|k|; Frobenius operator norm1",
        "local_compact_operator_symbol": C,
        "local_compact_construction": "For any compact smooth symmetric potential A_ij(t,x), apply Q_ik Q_jl-(1/2)Q_ij Q_kl with Q_ij=-Delta delta_ij+partial_i partial_j. The resulting psi is compact, transverse and traceless.",
        "compatible_detector_construction": "q=L psi is also compact TT. A nonzero compact psi cannot solve Lpsi=0 by zero-data uniqueness, hence gives nonzero q.",
        "zero_momentum_boundary": "Use a bounded representative at k=0, a measure-zero point for the stated L2 spaces. No division by |k| occurs in the energy estimate or local compact construction.",
        "noncompact_projection_boundary": "General TT projection need not preserve compact spatial support. The stress covariance extends by continuity in its proven tensor norm, not by falsely asserting compact support.",
        "checks": checks,
    }
