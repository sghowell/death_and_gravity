"""Full nonlocal vertex, both insertion positions and all nine products."""

from functools import cache
from itertools import combinations_with_replacement

import sympy as sp
from p8_vacuum_forward_loop import normalization as parent


@cache
def data():
    prior = parent.data()
    W = prior["full_quadratic_Hessian_insertion"]
    phi = sp.symbols("phi1 phi2", real=True)
    d11, d12, d22, p11, p12, p22 = sp.symbols("d11 d12 d22 p11 p12 p22", real=True)
    D = sp.Matrix([[d11, d12], [d12, d22]])
    P = sp.Matrix([[p11, p12], [p12, p22]])
    h = sp.Symbol("insertion_marker", real=True)
    original = -sp.trace((D + h * P) * W * (D + h * P) * W) / 4
    variation = -(sp.trace(P * W * D * W) + sp.trace(D * W * P * W)) / 4
    Q = {(i, j): W.diff(phi[i], phi[j]) for i in range(2) for j in range(2)}
    checks = {
        "literal_first_covariance_variation": sp.factor(
            sp.diff(original, h).subs(h, 0) - variation
        ),
        "cyclic_trace_half_factor": sp.factor(variation + sp.trace(P * W * D * W) / 2),
        "three_channels_two_positions_half_symmetry": 3 * 2 * sp.Rational(1, 2) - 3,
    }
    for indices in combinations_with_replacement(range(2), 4):
        i, j, k, l = indices
        expected = (
            -sum(
                sp.trace(P * Q[a, b] * D * Q[c, d])
                + sp.trace(D * Q[a, b] * P * Q[c, d])
                for a, b, c, d in ((i, j, k, l), (i, k, j, l), (i, l, j, k))
            )
            / 2
        )
        checks["four_external_derivatives_" + "".join(map(str, indices))] = sp.factor(
            sp.diff(variation, *(phi[n] for n in indices)) - expected
        )
    C, g = sp.symbols("external_vertex_C cubic_squared_g", real=True)
    h1, h2, k1, k2 = sp.symbols("left_h1 left_h2 right_h1 right_h2")
    terms = (
        ("bubble", C**2),
        ("triangle", C * g * h1),
        ("triangle", C * g * h2),
        ("triangle", C * g * k1),
        ("triangle", C * g * k2),
        ("box", g**2 * h1 * k1),
        ("box", g**2 * h1 * k2),
        ("box", g**2 * h2 * k1),
        ("box", g**2 * h2 * k2),
    )
    checks.update(
        {
            "complete_nonlocal_vertex_product": sp.expand(
                (C + g * (h1 + h2)) * (C + g * (k1 + k2))
                - sum(term for _, term in terms)
            ),
            "one_bubble_product": sum(kind == "bubble" for kind, _ in terms) - 1,
            "four_triangle_products": sum(kind == "triangle" for kind, _ in terms) - 4,
            "four_box_products": sum(kind == "box" for kind, _ in terms) - 4,
        }
    )
    return {
        "full_quadratic_Hessian": W,
        "quartic_vertex_tensor": prior["actual_quartic_vertex_tensor"],
        "insertion_positions": 2,
        "channels": 3,
        "symmetry_per_position": sp.Rational(1, 2),
        "vertex_product_terms": terms,
        "family": "scalar_Phi4_W2_F0 with P=-D F0,R D; all heavy contractions",
        "checks": checks,
    }
