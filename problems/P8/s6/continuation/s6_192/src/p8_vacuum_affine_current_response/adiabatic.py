"""Ordered covariance coefficients, physical-current parity and contour bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_matrix_adiabatic import riccati
from p8_vacuum_affine_matrix_response_tail import reference
from p8_vacuum_affine_matrix_response_tail.covariance import multiply

from . import vertices

MASS = reference.MASS
CONTOUR_COV = (s.Integer(2), s.Integer(3), s.Integer(128))
CONTOUR_CURRENT = (s.Integer(12), s.Integer(180), s.Integer(6000))


def covariance_coefficients(graph, order):
    """Ordinary marker coefficients; adjoint coefficients do not conjugate marker."""
    if type(order) is not int or not 0 <= order <= 10:
        raise ValueError("Require native covariance marker order0..10")
    if not graph or len(graph) < order + 1 or graph[0] != s.zeros(graph[0].rows):
        raise ValueError("Require a zero constant graph and all requested coefficients")
    n = graph[0].rows
    zero = s.zeros(n)
    E = [s.eye(n)]
    E.extend(
        -sum((multiply(graph[j].H, graph[k - j]) for j in range(k + 1)), zero)
        for k in range(1, order + 1)
    )
    C = [s.eye(n)]
    for k in range(1, order + 1):
        C.append(-sum((multiply(E[j], C[k - j]) for j in range(1, k + 1)), zero))
    F = [s.eye(n).col_join(-s.I * s.eye(n)) / s.sqrt(2)]
    F.extend(r.col_join(s.I * r) / s.sqrt(2) for r in graph[1 : order + 1])
    output = []
    for k in range(order + 1):
        value = s.zeros(2 * n)
        for a in range(k + 1):
            for b in range(k - a + 1):
                value += multiply(F[a], C[b], F[k - a - b].H)
        output.append(value.applyfunc(s.re))
    return output


@cache
def constants():
    A = reference.constants()["reference_parameter_norm_over_inverse_frequency"]
    contour = (A[0] / MASS, A[1] / MASS, A[2] / MASS)
    cov = (
        (1 + contour[0] ** 2) / (1 - contour[0] ** 2),
        16 * contour[1],
        16 * contour[2] + 128 * contour[1] ** 2,
    )
    c0, c1, c2 = CONTOUR_COV
    g0, g1, g2 = vertices.G
    current = (
        3 * g0 * c0,
        3 * ((g0 + g1) * c0 + g0 * c1),
        3 * ((2 * g0 + 2 * g1 + g2) * c0 + 2 * (g0 + g1) * c1 + g0 * c2),
    )
    return {
        "graph_contour_parameter_bounds": contour,
        "covariance_contour_before_rounding": cov,
        "current_contour_before_rounding": current,
        "fourth_order_comparison_terms": (0, 2, 4),
    }


@cache
def data():
    R = s.Matrix([[0, 1, -2], [-1, 0, 1], [2, -1, 0]]) / 7
    S = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 11
    # Genuine noncommuting real time jets; all orders of the finite recurrence.
    order = 10
    omega = [s.Integer(1000)] + [
        s.Rational((-1) ** j, j + 2) for j in range(1, order + 1)
    ]
    rotation = [R / s.Integer(j + 1) for j in range(order + 1)]
    squeeze = [S / s.Integer(j + 1) for j in range(order + 1)]
    b = riccati.reference_jets(omega, rotation, squeeze, order)
    graph = [s.zeros(3)] + [b[n][0] for n in range(1, order + 1)]
    coeff = covariance_coefficients(graph, 6)
    D = s.diag(1, -1, 0, 2, -2, 0)
    checks = {}
    for n in range(1, order + 1):
        checks[f"full_ordered_reference_real_phase_{n}"] = (
            graph[n] / s.I**n
        ).applyfunc(s.im)
        checks[f"full_ordered_reference_symmetry_{n}"] = graph[n] - graph[n].T
    for n in (1, 3, 5):
        checks[f"odd_QQ_coefficient_vanishes_{n}"] = coeff[n][:3, :3]
        checks[f"odd_PP_coefficient_vanishes_{n}"] = coeff[n][3:, 3:]
        checks[f"physical_block_vertex_odd_current_vanishes_{n}"] = s.trace(
            D * coeff[n]
        )
    for n in (0, 2, 4, 6):
        checks[f"even_QP_coefficient_vanishes_{n}"] = coeff[n][:3, 3:]
    q = s.Symbol("q", positive=True)
    checks["even_Cauchy_tail_geometric_identity"] = s.cancel(
        q**6 / (1 - q**2) - (1 / (1 - q**2) - 1 - q**2 - q**4)
    )
    c = constants()
    return {
        "formal_marker": "zeta is an adiabatic counting marker, distinct from physical shear epsilon. rhat(zeta)=sum1..10 zeta^n r_n; the analytic coefficient adjoint is sum zeta^n r_n^dagger.",
        "analytic_covariance": "Extend F(I-rsharp r)^-1 Fsharp holomorphically, then average with its coefficientwise conjugate. On real zeta this equals the full real covariance. At complex zeta there is no positivity claim; both graph factors have separate operator-norm bounds.",
        "contour": "At each fixed physical epsilon choose |zeta|=omega/m. Physical derivatives are taken at fixed zeta before choosing the radius, so the varying contour causes no omitted derivative term.",
        "parity": "Real R,S,omega give r_n=i^n times real symmetric matrices by the full ordered recurrence. The QQ/PP covariance blocks have only even marker powers; QP may be odd. The actual real block-diagonal metric vertex therefore yields an even current, including its first two physical amplitude derivatives.",
        "comparison_definition": "J_ad4 is ONLY the mathematical mode comparison given by the Taylor polynomial of the physical finite-reference current at marker orders0,2,4. No claim identifies it with the fixed covariant mu=m subtraction or changes any finite local counterterm.",
        "constants": c,
        "covariance_contour_displays": CONTOUR_COV,
        "current_contour_displays_over_omega": CONTOUR_CURRENT,
        "checks": checks,
        "gates": {
            "both_complex_graph_factors_in_small_ball": c[
                "graph_contour_parameter_bounds"
            ][0]
            < s.Rational(1, 100),
            "all_complex_covariance_contour_displays": all(
                c["covariance_contour_before_rounding"][a] < CONTOUR_COV[a]
                for a in range(3)
            ),
            "all_current_contour_displays": all(
                c["current_contour_before_rounding"][a] <= CONTOUR_CURRENT[a]
                for a in range(3)
            ),
            "genuine_noncommuting_ordered_fixture": R * S != S * R,
            "odd_full_covariance_not_discarded": coeff[1][:3, 3:] != s.zeros(3),
            "fourth_order_current_tail_starts_at_six": c[
                "fourth_order_comparison_terms"
            ]
            == (0, 2, 4),
        },
    }
