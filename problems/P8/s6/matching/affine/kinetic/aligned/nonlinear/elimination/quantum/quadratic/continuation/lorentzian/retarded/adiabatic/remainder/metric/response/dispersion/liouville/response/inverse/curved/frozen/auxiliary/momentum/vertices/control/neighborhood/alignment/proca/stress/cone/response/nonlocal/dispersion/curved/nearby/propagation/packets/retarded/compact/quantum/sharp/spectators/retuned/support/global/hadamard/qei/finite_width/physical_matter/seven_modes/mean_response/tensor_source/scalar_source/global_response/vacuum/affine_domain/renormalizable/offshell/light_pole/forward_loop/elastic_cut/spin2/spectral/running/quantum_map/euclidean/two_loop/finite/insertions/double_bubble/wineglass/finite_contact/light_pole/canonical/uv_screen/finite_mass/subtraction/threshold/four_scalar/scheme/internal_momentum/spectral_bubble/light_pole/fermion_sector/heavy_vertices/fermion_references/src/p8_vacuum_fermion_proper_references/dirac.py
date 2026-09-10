"""Euclidean Dirac and inverse-resolvent normalization of proper subgraphs."""

from functools import cache

import sympy as sp
from p8_vacuum_two_loop_fermion_ledger import gauge as color


def gamma_matrices():
    zero = sp.zeros(2)
    identity = sp.eye(2)
    sigma = (
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.Matrix([[1, 0], [0, -1]]),
    )
    return tuple(
        zero.row_join(-sp.I * s).col_join((sp.I * s).row_join(zero)) for s in sigma
    ) + (zero.row_join(identity).col_join(identity.row_join(zero)),)


@cache
def data():
    gamma = gamma_matrices()
    q = sp.symbols("q0:4", real=True)
    m = sp.Symbol("m", positive=True)
    slash = sum((a * b for a, b in zip(q, gamma)), sp.zeros(4))
    numerator = m * sp.eye(4) - sp.I * slash
    denom = m * m + sum(a * a for a in q)
    checks = {}
    for i in range(4):
        for j in range(4):
            checks[f"Clifford_{i}_{j}"] = (
                gamma[i] * gamma[j] + gamma[j] * gamma[i] - 2 * (i == j) * sp.eye(4)
            )
    checks["free_Euclidean_Dirac_inverse"] = (
        m * sp.eye(4) + sp.I * slash
    ) * numerator - denom * sp.eye(4)
    contraction = sum((a * numerator * a for a in gamma), sp.zeros(4))
    checks["literal_four_dimensional_gauge_numerator"] = contraction - (
        2 * sp.I * slash + 4 * m * sp.eye(4)
    )
    matrices = color.data()["fundamental_generators"]
    checks["open_fermion_line_Casimir"] = sum(
        (a * a for a in matrices), sp.zeros(3)
    ) - sp.Rational(4, 3) * sp.eye(3)
    A = sp.Matrix([[2, 1], [1, 3]])
    V = sp.Matrix([[1, 2], [2, 0]])
    h = sp.Symbol("Gaussian_variance", real=True)
    S = A.inv()
    averaged = S + h * S * V * S * V * S
    inverse = A - h * V * S * V
    checks["Gaussian_resolvent_inverse_first_order"] = (
        (averaged * inverse).diff(h).subs(h, 0)
    )
    checks["scalar_vertex_sign"] = -(sp.Integer(1) ** 2) + 1
    checks["imaginary_gauge_vertex_sign"] = -(sp.I**2) - 1
    d = sp.Symbol("dimension")
    return {
        "Euclidean_gamma_matrices": gamma,
        "open_line_Casimir": sp.Rational(4, 3),
        "formal_dimensional_gauge_numerator": "gamma_mu(-i slash(q)+m)gamma_mu=(2-d)(-i slash(q))+d m",
        "scalar_inverse_correction": "-Y integral D_Phi(k) S(p-k)",
        "gauge_inverse_correction": "+a Cf integral D_A(k) gamma_mu S(p-k) gamma_mu in Feynman gauge",
        "formal_dimension": d,
        "scope": "The d-dimensional Clifford identity follows before setting d=4-2epsilon. Literal four-dimensional matrices only independently check its four-dimensional specialization; the regulator factors remain in the anchors.",
        "checks": checks,
    }
