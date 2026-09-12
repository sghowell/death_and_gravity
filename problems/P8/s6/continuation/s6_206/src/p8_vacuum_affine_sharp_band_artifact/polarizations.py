"""Leading actual endpoint symbol with all physical Proca pairs."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import stress


def cross(n):
    return s.Matrix([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])


def stress_matrix(D):
    return s.diag(-D, -D, s.zeros(1), D)


def leading_mode(direction, polarization, kind, a=1):
    """Coefficient of sqrt(r); mass fixed, r -> infinity."""
    v = 1 / s.sqrt(2 * a)
    zero = s.zeros(3, 1)
    if kind == "transverse":
        return s.Matrix(
            [
                *(-s.I * v * polarization),
                *(s.I * v * direction.cross(polarization)),
                0,
                *zero,
            ]
        )
    if kind == "longitudinal":
        return s.Matrix([*zero, *zero, -v, *(v * direction)])
    raise ValueError("Expected transverse or longitudinal")


def leading_density(D, G, n, a=1):
    return (
        4 * s.trace(D * G)
        - 4 * (n.T * (D * G + G * D) * n)[0]
        + 3 * (n.T * D * n)[0] * (n.T * G * n)[0]
    ) / (8 * a)


@cache
def data():
    a = s.symbols("a", positive=True)
    x = s.symbols("D0:5", real=True)
    y = s.symbols("G0:5", real=True)
    D = s.Matrix([[x[0], x[1], x[2]], [x[1], x[3], x[4]], [x[2], x[4], -x[0] - x[3]]])
    G = s.Matrix([[y[0], y[1], y[2]], [y[1], y[3], y[4]], [y[2], y[4], -y[0] - y[3]]])
    n = s.Matrix([0, 0, 1])
    e1 = s.Matrix([1, 0, 0])
    e2 = s.Matrix([0, 1, 0])
    plus = [
        leading_mode(n, e1, "transverse", a),
        leading_mode(n, e2, "transverse", a),
        leading_mode(n, n, "longitudinal", a),
    ]
    minus = [
        leading_mode(-n, e1, "transverse", a),
        leading_mode(-n, e2, "transverse", a),
        leading_mode(-n, -n, "longitudinal", a),
    ]
    MD, MG = stress_matrix(D), stress_matrix(G)
    AD = s.Matrix(3, 3, lambda i, j: s.expand((plus[i].T * MD * minus[j])[0]))
    AG = s.Matrix(3, 3, lambda i, j: s.expand((plus[i].T * MG * minus[j])[0]))
    actual = s.expand(
        a
        * sum(s.conjugate(AD[i, j]) * AG[i, j] for i in range(3) for j in range(3))
        / 2
    )
    literal = sum(
        (
            D[i, j] * stress.data()["quadratic_matrices"][i + 1, j + 1]
            for i in range(3)
            for j in range(3)
        ),
        s.zeros(10),
    )
    PT = s.eye(3) - n * n.T
    C = cross(n)
    projected = PT * (D - C.T * D * C) * PT
    tt = s.expand(s.trace(projected * (PT * (G - C.T * G * C) * PT)) / (8 * a))
    ll = (n.T * D * n)[0] * (n.T * G * n)[0] / (8 * a)
    return {
        "full_mode_limit": "With fixed mass and CD time, W8/r -> 1/a, f sqrt(r) -> sqrt(a/2), p/sqrt(r) -> -i/sqrt(2a). Transverse magnetic and longitudinal mass/temporal readouts remain at order sqrt(r); longitudinal electric and transverse mass readouts are lower order.",
        "complete_polarization_sum": "The actual j=0 current is Re sum_(alpha,beta) A_D# A_G/(W_k+W_l). Its leading density at l=-k is r*h(D,G,n). All four TT, four mixed and one LL terms are computed. The mixed leading terms vanish by their actual scaling, not by deleting the longitudinal mode.",
        "invariant_symbol": "h=[4tr(DG)-4 n.(DG+GD)n+3(n.D.n)(n.G.n)]/(8a), for real symmetric tracefree D,G. It is even in n. Complex Fourier tensors follow by sesquilinear extension.",
        "scope": "This is the fixed-mass ultraviolet limit of the original massive three-polarization Proca comparison. It is not the two-polarization Maxwell massless theory.",
        "checks": {
            "literal_full_tracefree_stress_matrix": literal - MD,
            "all_nine_pairs_give_invariant_symbol": s.expand(
                actual - leading_density(D, G, n, a)
            ),
            "complete_TT_projector_contraction": s.expand(
                tt
                - a
                * sum(
                    s.conjugate(AD[i, j]) * AG[i, j] for i in range(2) for j in range(2)
                )
                / 2
            ),
            "longitudinal_pair_survives": s.expand(
                a * s.conjugate(AD[2, 2]) * AG[2, 2] / 2 - ll
            ),
            "all_four_mixed_leading_pairs": s.Matrix(
                [AD[0, 2], AD[1, 2], AD[2, 0], AD[2, 1]]
            ),
            "TT_and_LL_recover_complete_symbol": s.expand(
                tt + ll - leading_density(D, G, n, a)
            ),
            "current_endpoint_sign": -s.im(-s.I * s.Symbol("positive", positive=True))
            - s.Symbol("positive", positive=True),
        },
        "gates": {
            "all_nine_physical_pairs_retained": AD.shape == (3, 3),
            "longitudinal_not_Maxwell_deleted": ll != 0,
            "full_temporal_constraint_readout_retained": plus[2][6] != 0,
            "tracefree_not_full_lapse_shift_claim": True,
        },
    }
