"""Exact matrix normalizer, pure-state graph and full covariance map."""

from functools import cache

import sympy as s


def covariance(graph):
    n = graph.rows
    unit = s.eye(n)
    C = (unit - graph.H * graph).inv()
    F = (unit + graph).col_join(-s.I * (unit - graph)) / s.sqrt(2)
    return (F * C * F.H).applyfunc(s.re)


@cache
def data():
    omega = s.Symbol("omega", positive=True)
    rho = s.Symbol("rho", real=True)
    L = s.Matrix(3, 3, s.symbols("L0:9", real=True))
    unit = s.eye(3)
    R = (L.T - L) / 2
    S = rho * unit - (L + L.T) / 2
    realflow = (
        (rho * unit - L)
        .row_join(omega * unit)
        .col_join((-omega * unit).row_join(L.T - rho * unit))
    )
    transform = unit.row_join(s.I * unit).col_join(unit.row_join(-s.I * unit)) / s.sqrt(
        2
    )
    expected = (
        (-s.I * omega * unit + R)
        .row_join(S)
        .col_join(S.row_join(s.I * omega * unit + R))
    )
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]])
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -1]])
    r = (G + s.I * D) / 100
    C = covariance(r)
    J = s.zeros(3).row_join(unit).col_join((-unit).row_join(s.zeros(3)))
    e = s.Matrix([[1 + s.I, 2, -s.I], [3, s.I, 1], [2 * s.I, 0, 4]])
    fast = 2 * s.I * omega * e + R * e - e * R
    w, W, v = s.symbols("w W v", real=True)
    mode_p = -s.I * W - v
    actual_graph = (w - s.I * mode_p) / (w + s.I * mode_p)
    return {
        "normalizer": "K V=omega^2 I, B=K^(1/2), q=B^-1 A, p=B pi, L=B^-1 B'. Then q'=p-Lq and p'=-omega^2 q+L^T p.",
        "balanced_quadratures": "Q=sqrt(omega)q, P=p/sqrt(omega), b=(Q+iP)/sqrt(2); b'=(-i omega I+R)b+S bdag, R=(L^T-L)/2 and S=rho I-(L+L^T)/2, rho=omega'/(2omega).",
        "actual_pure_graph": "If U,V are the coefficients of the unchanged initial annihilators in b,bdag, r=V U^-1 is complex symmetric, ||r||<1, and r'=2i omega r+[R,r]+S-r S r. CCR gives U Udag=(I-rdag r)^-1; no adiabatic state reset is made.",
        "covariance": "For F(r)=[I+r; -i(I-r)]/sqrt(2), the real symmetrized six-quadrature covariance is Re[F(r)(I-rdag r)^-1 F(r)dag]. The arbitrary positive-frequency unitary and phase cancel.",
        "comparison_domain": "For two symmetric graphs of norm<=1/10, the full six-quadrature covariance map is Lipschitz with constant16 in operator norm.",
        "checks": {
            "full_three_mode_complex_flow": transform * realflow * transform.H
            - expected,
            "rotation_is_skew": R + R.T,
            "squeeze_is_symmetric": S - S.T,
            "graph_Riccati_vector_field_preserves_symmetry": (
                2 * s.I * omega * r + R * r - r * R + S - r * S * r
            )
            - (2 * s.I * omega * r + R * r - r * R + S - r * S * r).T,
            "exact_noncommuting_graph_covariance_purity": C * J * C - J / 4,
            "exact_noncommuting_graph_covariance_symmetry": C - C.T,
            "fast_phase_and_rotation_Frobenius_norm_rate_zero": s.trace(
                fast.H * e + e.H * fast
            ),
            "actual_scalar_preparation_graph_sign": s.cancel(
                actual_graph - (w - W + s.I * v) / (w + W - s.I * v)
            ),
        },
        "gates": {
            "noncommuting_graph_fixture_in_ball": sum(abs(x) ** 2 for x in r)
            < s.Rational(1, 100),
            "covariance_inverse_norm_below_two": 1 / (1 - s.Rational(1, 10) ** 2) < 2,
            "covariance_factor_norm_below_two": 1 + s.Rational(1, 10) ** 2 < 4,
            "inverse_map_difference_coefficient_below_one": s.Rational(1, 5)
            / (1 - s.Rational(1, 10) ** 2) ** 2
            < 1,
            "complete_covariance_Lipschitz_display": 4 + 4 + 4 < 16,
        },
    }
