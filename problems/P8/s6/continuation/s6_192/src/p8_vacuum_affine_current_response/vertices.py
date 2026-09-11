"""Physical metric vertices and their complete balanced-frame derivatives."""

from functools import cache
from math import comb

import sympy as s
from p8_vacuum_affine_matrix_response_tail import mixed
from p8_vacuum_affine_matrix_response_tail.covariance import multiply

G = (s.Integer(2), s.Integer(25), s.Integer(650))
FRAME = (s.Integer(1), s.Integer(5), s.Integer(110))


def vertex_jet(frame, omega, vertex):
    """Raw parameter derivatives of T^t M_D T / omega through order two."""
    if any(len(row) != 3 for row in (frame, omega, vertex)):
        raise ValueError("Require three complete raw amplitude jets")
    inv = [1 / omega[0]]
    for n in range(1, 3):
        inv.append(
            -sum(comb(n, j) * omega[j] * inv[n - j] for j in range(1, n + 1)) / omega[0]
        )
    force = []
    for n in range(3):
        value = s.zeros(frame[0].cols)
        for a in range(n + 1):
            for b in range(n - a + 1):
                c = n - a - b
                value += (
                    comb(n, a)
                    * comb(n - a, b)
                    * multiply(frame[a].T, vertex[b], frame[c])
                )
        force.append(value.applyfunc(s.expand))
    return [
        sum(
            (comb(n, j) * inv[j] * force[n - j] for j in range(n + 1)),
            s.zeros(frame[0].cols),
        ).applyfunc(s.expand)
        for n in range(3)
    ]


@cache
def constants():
    c = mixed.constants()
    root1, root2 = c["root"][0, 1], c["root"][0, 2]
    inv1, inv2 = c["inverse_root"][0, 1], c["inverse_root"][0, 2]
    w1, w2 = c["omega"][0, 1], c["omega"][0, 2]
    sqrt1 = w1 / 2
    sqrt2 = w2 / 2 + w1**2 / 4
    inverse_sqrt2 = w2 / 2 + 3 * w1**2 / 4
    first = max(root1, inv1) + sqrt1
    second = max(
        root2 + 2 * sqrt1 * root1 + inverse_sqrt2, inv2 + 2 * sqrt1 * inv1 + sqrt2
    )
    F1 = 2 * FRAME[1] * G[0] + 2
    F2 = 2 * FRAME[2] * G[0] + 2 * FRAME[1] ** 2 * G[0] + 4 * FRAME[1] * 2 + 2
    return {
        "root_parameter_bounds": (root1, root2),
        "inverse_root_parameter_bounds": (inv1, inv2),
        "balanced_frame_parameter_bounds_before_rounding": (first, second),
        "vertex_parameter_bounds_before_rounding": (
            F1 + w1 * G[0],
            F2 + 2 * w1 * F1 + (2 * w1**2 + w2) * G[0],
        ),
    }


@cache
def data():
    A = s.Matrix([[1, 2, 0], [2, -1, 1], [0, 1, 0]]) / 7
    D = s.Matrix([[0, 1, 2], [1, 2, -1], [2, -1, -2]]) / 9
    T = [
        s.diag(2, 3, 5, s.Rational(1, 2), s.Rational(1, 3), s.Rational(1, 5)),
        s.diag(A, D),
        s.diag(D, A),
    ]
    V = [s.diag(D, A), s.diag(A, D), s.diag(A + D, A - D)]
    w = [s.Integer(3), s.Rational(2, 5), s.Rational(-1, 7)]
    actual = vertex_jet(T, w, V)
    H1, H2 = multiply(T[0].inv(), T[1]), multiply(T[0].inv(), T[2])
    g = actual[0]
    md1 = multiply(T[0].T, V[1], T[0]) / w[0]
    md2 = multiply(T[0].T, V[2], T[0]) / w[0]
    f1 = multiply(H1.T, g) + multiply(g, H1) + md1
    f2 = (
        multiply(H2.T, g)
        + multiply(g, H2)
        + 2 * multiply(H1.T, g, H1)
        + 2 * multiply(H1.T, md1)
        + 2 * multiply(md1, H1)
        + md2
    )
    p1, p2 = w[1] / w[0], w[2] / w[0]
    eps = s.Symbol("epsilon", real=True)
    polT = sum((eps**j * T[j] / s.factorial(j) for j in range(3)), s.zeros(6))
    polV = sum((eps**j * V[j] / s.factorial(j) for j in range(3)), s.zeros(6))
    reciprocal = (
        1 / w[0]
        - eps * w[1] / w[0] ** 2
        + eps**2 * (2 * w[1] ** 2 / w[0] ** 3 - w[2] / w[0] ** 2) / 2
    )
    literal = (multiply(polT.T, polV, polT) * reciprocal).applyfunc(s.expand)
    checks = {
        "full_first_vertex_frame_frequency_derivative": actual[1] - f1 + p1 * g,
        "full_second_vertex_frame_frequency_derivative": actual[2]
        - f2
        + 2 * p1 * f1
        - (2 * p1**2 - p2) * g,
        "noncommuting_detector_and_direction": s.trace(A * D - D * A),
    }
    for n in range(3):
        checks[f"literal_full_vertex_product_jet_{n}"] = (
            literal.applyfunc(lambda x, n=n: s.factorial(n) * x.coeff(eps, n))
            - actual[n]
        )
        checks[f"real_block_diagonal_vertex_jet_{n}"] = actual[n][:3, 3:]
    c = constants()
    return {
        "physical_current": "J_D=-1/2 tr(M_D C_phys)=-omega/2 tr(G_D Sigma), M=diag(V,K), G_D=T^t M_D T/omega. This is the actual conjugate metric current for a fixed real symmetric tracefree shear-coordinate detector D, not a scalar surrogate.",
        "vertex_domain": "D is smooth and fixed, with ||D(t)||op<=1; gamma=epsilon Gamma, ||Gamma^(j)||op<=1 through j12, |epsilon|<=.01. The positive constraint term remains in the base Hamiltonian.",
        "ordered_exponential_bound": "For a0..2, energy-relative M_(D Gamma^a) norm <=exp(2 delta)<2, using the full ordered exponential Frechet simplex and the positive mass/electric/magnetic pieces. No commutation with the momentum constraint is assumed.",
        "frame": "T=diag(B/sqrt(omega),sqrt(omega) B^-1), B=K^(1/2), T^t M T=omega I. Raw amplitude derivatives use T^-1 T_e and T^-1 T_ee, not derivatives of a commuting scalar root.",
        "constants": c,
        "balanced_frame_parameter_displays": FRAME,
        "physical_vertex_parameter_displays": G,
        "checks": checks,
        "gates": {
            "root_and_inverse_raw_second_jets": c["root_parameter_bounds"] == (4, 68)
            and c["inverse_root_parameter_bounds"] == (4, 100),
            "frame_first_second_displays": all(
                c["balanced_frame_parameter_bounds_before_rounding"][j] < FRAME[j + 1]
                for j in range(2)
            ),
            "complete_vertex_first_second_displays": all(
                c["vertex_parameter_bounds_before_rounding"][j] < G[j + 1]
                for j in range(2)
            ),
            "ordered_exponential_relative_display": 1 / (1 - s.Rational(1, 50)) < 2,
            "noncommuting_matrix_fixture": A * D != D * A,
        },
    }
