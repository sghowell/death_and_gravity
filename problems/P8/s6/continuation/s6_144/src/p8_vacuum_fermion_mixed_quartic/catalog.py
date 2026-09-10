"""Literal operator variation, mixed weights and complex pair-denominator gap."""

from functools import cache
from itertools import combinations_with_replacement, permutations

import sympy as sp
from p8_vacuum_forward_loop import normalization as prior


@cache
def data():
    W = prior.data()["full_quadratic_Hessian_insertion"]
    phi = sp.symbols("phi1 phi2", real=True)
    v1111, v1112, v1122, v1222, v2222 = sp.symbols("v1111 v1112 v1122 v1222 v2222")
    potential = (
        v1111 * phi[0] ** 4
        + 4 * v1112 * phi[0] ** 3 * phi[1]
        + 6 * v1122 * phi[0] ** 2 * phi[1] ** 2
        + 4 * v1222 * phi[0] * phi[1] ** 3
        + v2222 * phi[1] ** 4
    ) / 24
    F = sp.hessian(potential, phi)
    d11, d12, d22, h = sp.symbols("d11 d12 d22 h")
    D = sp.Matrix([[d11, d12], [d12, d22]])
    original = -sp.trace(D * (W + h * F) * D * (W + h * F)) / 4
    mixed = -sp.trace(D * W * D * F) / 2
    Q = {(i, j): W.diff(phi[i], phi[j]) for i in range(2) for j in range(2)}
    V = {(i, j): F.diff(phi[i], phi[j]) for i in range(2) for j in range(2)}
    checks = {
        "literal_Hessian_variation": sp.factor(sp.diff(original, h).subs(h, 0) - mixed)
    }
    for indices in combinations_with_replacement(range(2), 4):
        i, j, k, l = indices
        expected = (
            -sum(
                sp.trace(D * Q[a, b] * D * V[c, d])
                + sp.trace(D * Q[c, d] * D * V[a, b])
                for a, b, c, d in ((i, j, k, l), (i, k, j, l), (i, l, j, k))
            )
            / 2
        )
        checks["four_derivative_" + "".join(map(str, indices))] = sp.factor(
            sp.diff(mixed, *(phi[n] for n in indices)) - expected
        )
    # Distinct internal labels A,B; two background legs 1,2. Fix A cyclically.
    words = tuple(("A",) + w for w in permutations(("B", "1", "2")))
    arcs = [(w.index("B"), 4 - w.index("B")) for w in words]
    x, a, b, t = sp.symbols("x real_z imag_z boundary_displacement", real=True)
    S = x + 1
    pair = (S - (a + sp.I * b) / 4) ** 2 + (a + sp.I * b) * sp.Symbol(
        "q_axis_squared", real=True
    )
    # Drop the nonnegative a*q_axis_squared, use b^2<=1-t^2 and a=2+t.
    low = (S - (2 + t) / 4) ** 2 - (1 - t * t) / 16
    min_boundary = (S - sp.Rational(3, 4)) ** 2
    checks.update(
        {
            "six_cyclic_boxes": len(words) - 6,
            "two_arc_types_each": sum(
                arcs.count(ab) == 2 for ab in ((1, 3), (2, 2), (3, 1))
            )
            - 3,
            "three_channels_two_assignments_half": 3 * 2 * sp.Rational(1, 2) - 3,
            "disc_gap_endpoint_difference": sp.factor(
                low - min_boundary - (1 - t) * (8 * S - 6 - 2 * t) / 16
            ),
            "disc_gap_radial_positive_slack": sp.factor(
                min_boundary - S * S / 16 - 3 * (S - 1) * (5 * S - 3) / 16
            ),
            "bubble_difference_numerator_cap": sp.Rational(3, 2)
            + sp.Rational(9, 16)
            - sp.Rational(33, 16),
            "heavy_Neumann_ratio_at_minimum_M": sp.Rational(18, 100)
            + sp.Rational(18**2, 10000)
            - sp.Rational(531, 2500),
            "heavy_ratio_below_half_slack": sp.Rational(1, 2)
            - sp.Rational(531, 2500)
            - sp.Rational(719, 2500),
        }
    )
    return {
        "family": "scalar_Phi4_W1_F2",
        "operator": mixed,
        "cyclic_box_words": words,
        "fermion_arc_counts": arcs,
        "channels": 3,
        "assignments_per_channel": 2,
        "coefficient_per_assignment": sp.Rational(1, 2),
        "vertex": "C(z)+g(h1+h2), C(z)=-L+g/(M-z)",
        "full_pair_denominator": pair,
        "complex_gap": "For |z-2|<=1 or z=0, |D_+ D_-|<=16/(1+q^2)^2; each internal heavy line <=2/(M+q^2) for M>=10000 and shift one-norm<18.",
        "checks": checks,
    }
