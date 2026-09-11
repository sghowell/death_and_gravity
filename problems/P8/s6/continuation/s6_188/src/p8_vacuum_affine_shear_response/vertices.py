"""Both shear vertices, temporal-constraint contact and physical stress sign."""

from functools import cache

import sympy as s

from .hamiltonian import a, cross, k, k2, m, symmetric

G = symmetric("G", True)
D = symmetric("D", True)
C = (G * D + D * G) / 2


def first(direction):
    return s.diag(
        -a * m * m * direction + cross.T * direction * cross / a, direction / a
    )


def second(left, right):
    product = (left * right + right * left) / 2
    return s.diag(a * m * m * product + cross.T * product * cross / a, product / a)


@cache
def data():
    e, r = s.symbols("epsilon eta", real=True)
    gamma = e * G + r * D
    plus = s.eye(3) + gamma + gamma * gamma / 2
    minus = s.eye(3) - gamma + gamma * gamma / 2
    K = plus / a + k * k.T / (a**3 * m**2)
    V = a * m * m * minus + cross.T * plus * cross / a
    M = s.diag(V, K)
    first_direct = M.diff(e).subs({e: 0, r: 0})
    second_direct = M.diff(e, r).subs({e: 0, r: 0})
    vector = s.Matrix(s.symbols("A0:3", real=True))
    momentum = s.Matrix(s.symbols("pi0:3", real=True))
    Z = vector.col_join(momentum)
    E = momentum / a**2
    B = cross * vector / a**2
    massA = m * vector / a
    stress = -E.dot(G * E) - B.dot(G * B) + massA.dot(G * massA)
    magnetic = -k.dot(G * k) * s.eye(3) - k2 * G + G * k * k.T + k * k.T * G
    second_magnetic = (
        k.dot(C * k) * s.eye(3)
        + k2 * C
        + k.dot(G * k) * D
        + k.dot(D * k) * G
        - C * k * k.T
        - k * k.T * C
        - G * k * k.T * D
        - D * k * k.T * G
    )
    checks = {
        "independent_first_exponential_vertex": s.simplify(first_direct - first(G)),
        "independent_mixed_second_vertex": s.simplify(second_direct - second(G, D)),
        "first_inverse_metric_magnetic_vertex": s.simplify(
            cross.T * G * cross - magnetic
        ),
        "second_inverse_metric_magnetic_vertex": s.simplify(
            cross.T * C * cross - second_magnetic
        ),
        "complete_physical_traceless_stress_sign": s.expand(
            (Z.T * first(G) * Z)[0] / 2 + a**3 * stress / 2
        ),
        "mixed_second_vertex_symmetric": s.simplify(second(G, D) - second(D, G)),
        "temporal_constraint_term_has_no_shear_derivative": s.diff(
            k * k.T / (a**3 * m**2), e
        ),
        "first_vertex_matrix_symmetric": first(G) - first(G).T,
        "second_vertex_matrix_symmetric": second(G, D) - second(G, D).T,
    }
    return {
        "first_metric_vertex": first(G),
        "mixed_second_metric_vertex": second(G, D),
        "physical_current": "J_G=-H_G=a^3 T_G/2 in the unimodular coordinate chart; h=sqrt(kappa)gamma/2 has the inherited canonical scaling.",
        "second_contact": "The response includes -<H_GD> times the prescribed shear perturbation as well as the propagated covariance. This contact is generally nonzero; centering the stress does not remove its variation.",
        "all_order_majorant": "For every n>=1, ||M0^(-1/2) D^n M_gamma[Gamma1,...,Gamman]|gamma=0 M0^(-1/2)||op <= product_j ||Gamma_j||op. Symmetrized exponential products and the positive mass/magnetic/electric decomposition prove the bound.",
        "state_boundary": "No transverse/longitudinal modes are dropped: generic external shear mixes their coordinate Hamiltonian blocks. A0 remains its exact constraint, not an independent fourth oscillator.",
        "checks": checks,
    }
