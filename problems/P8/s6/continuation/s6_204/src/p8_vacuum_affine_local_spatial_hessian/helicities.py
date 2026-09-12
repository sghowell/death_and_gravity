"""Complete tracefree spatial space and smooth polynomial momentum maps."""

from functools import cache

import sympy as s


def tf(M):
    return s.ImmutableMatrix(M - s.eye(3) * s.trace(M) / 3)


def B(P, M):
    v = M * P
    return tf((P * v.T + v * P.T) / 2)


def C(P, M):
    return tf(P * P.T) * (P.T * M * P)[0]


@cache
def basis():
    return (
        s.ImmutableMatrix(s.diag(1, -1, 0) / s.sqrt(2)),
        s.ImmutableMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / s.sqrt(2),
        s.ImmutableMatrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / s.sqrt(2),
        s.ImmutableMatrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / s.sqrt(2),
        s.ImmutableMatrix(s.diag(-1, -1, 2) / s.sqrt(6)),
    )


@cache
def data():
    p = s.Symbol("p", real=True)
    q = p * p
    P = s.Matrix([0, 0, p])
    forms = basis()
    b = (0, 0, q / 2, q / 2, 2 * q / 3)
    c = (0, 0, 0, 0, 2 * q * q / 3)
    checks = {
        "complete_orthonormal_tracefree_basis": s.Matrix(
            5, 5, lambda i, j: s.trace(forms[i] * forms[j])
        )
        - s.eye(5),
        "all_five_traces_zero": s.Matrix([s.trace(M) for M in forms]),
    }
    for j, M in enumerate(forms):
        checks[f"full_B_eigenvector_{j}"] = B(P, M) - b[j] * M
        checks[f"full_C_eigenvector_{j}"] = C(P, M) - c[j] * M
    x = s.symbols("x0:5")
    M = sum((xj * T for xj, T in zip(x, forms)), s.zeros(3))
    checks["zero_momentum_B"] = B(s.zeros(3, 1), M)
    checks["zero_momentum_C"] = C(s.zeros(3, 1), M)
    return {
        "space": "All five symmetric tracefree spatial metric directions are retained: two tensor, two vector and one scalar shear relative to each nonzero external momentum. These are spatial helicity labels, not fully constraint-reduced propagating modes.",
        "smooth_maps": "For q=P.P, define B_P(M)=TF sym(P tensor(MP)) and C_P(M)=TF(P tensor P)(P.M.P). Both are exact polynomial maps, with no division by |P| or a projector light-cone denominator. Their zero-momentum limits are unambiguous.",
        "basis": forms,
        "B_eigenvalues": b,
        "C_eigenvalues": c,
        "checks": checks,
        "gates": {
            "five_complete_tracefree_spatial_directions": len(forms) == 5,
            "all_basis_matrices_symmetric": all(M == M.T for M in forms),
            "both_vector_directions_present": b[2] == b[3] and b[2] != 0,
            "scalar_shear_not_dropped": c[4] != 0,
        },
    }
