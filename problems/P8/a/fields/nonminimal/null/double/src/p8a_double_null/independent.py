"""Fraction-only compact-polynomial integration and selected boost replay."""

from fractions import Fraction as F


def derivative(p):
    return [F(i) * p[i] for i in range(1, len(p))] or [F(0)]


def pairing(p, q):
    return sum(
        (
            a * b * F(2, i + j + 1)
            for i, a in enumerate(p)
            for j, b in enumerate(q)
            if (i + j) % 2 == 0
        ),
        F(0),
    )


def replay():
    p = [F(1), F(0), F(-2), F(0), F(1)]
    profiles = [p, derivative(p), derivative(derivative(p))]
    norm = F(315, 256)
    pair = [[norm * pairing(a, b) for b in profiles] for a in profiles]
    orders = [(2, 0), (1, 1), (0, 2)]
    gram = [[pair[i[0]][j[0]] * pair[i[1]][j[1]] for j in orders] for i in orders]
    xi, r = F(1, 6), F(4, 3)
    h1, h2 = pair[1][1], pair[2][2]
    cost = (
        F(2, 3) * (1 + 4 * xi) * r * h2
        + 4 * (1 - 2 * xi) * h1 * h1 / r
        + F(2, 3) * h2 / r**3
    )
    return {
        "shape_norm_squared": pair[0][0],
        "shape_first_norm_squared": h1,
        "shape_second_norm_squared": h2,
        "actual_shape_derivative_pairings": pair,
        "actual_product_derivative_Gram": gram,
        "selected_boost_ratio": r,
        "selected_quantum_cost": cost,
        "selected_quantum_coefficient_without_hbar_over_pi_squared": cost / 8,
        "strict_margin_below_74": F(74) - cost,
        "selected_state_coefficient": 8 * xi * h1,
        "boost_stationary_polynomial_at_1": F(35) - 24 - 63,
        "boost_stationary_polynomial_at_4_over_3": 35 * r**4 - 24 * r * r - 63,
    }
