"""Full transformed interaction degrees and explicit one-loop topology."""

from functools import cache
from itertools import product

import sympy as sp


def topologies(external, loops):
    if type(external) is not int or external not in (2, 4, 6):
        raise ValueError("Only native external counts two, four or six are in scope")
    if type(loops) is not int or not 0 <= loops <= 2:
        raise ValueError("Only native loop orders zero through two are in scope")
    budget = 2 * loops + external - 2
    return tuple(
        counts
        for counts in product(range(budget // 2 + 1), repeat=5)
        if sum((valence - 2) * n for valence, n in zip((4, 6, 8, 10, 12), counts))
        == budget
    )


@cache
def data():
    x, y, r, s, t = sp.symbols("site_x site_y cubic_R_x cubic_R_y field_degree")
    k11, k12, k22 = sp.symbols("free_K11 free_K12 free_K22")
    h11, h12, h22 = sp.symbols(
        "full_heavy_inverse11 full_heavy_inverse12 full_heavy_inverse22"
    )
    L, g = sp.symbols("polynomial_quartic cubic_squared")
    a, b = sp.symbols("quadratic_counterterm quartic_counterterm")
    Q = (k11 * x * x + 2 * k12 * x * y + k22 * y * y) / 2
    V = (
        L * (x**4 + y**4) / 24
        - g * (h11 * x**4 + 2 * h12 * x * x * y * y + h22 * y**4) / 8
    )
    F = {x: t * x + t**3 * r, y: t * y + t**3 * s}
    full = sp.Poly(sp.expand((Q + V).subs(F, simultaneous=True)), t)
    degrees = {n: full.coeff_monomial(t**n) for n in (2, 4, 6, 8, 10, 12)}
    direction = lambda f: sp.diff(f, x) * r + sp.diff(f, y) * s
    expected = {
        2: Q,
        4: V + direction(Q),
        6: direction(V) + direction(direction(Q)) / 2,
        8: direction(direction(V)) / 2,
        10: direction(direction(direction(V))) / 6,
        12: direction(direction(direction(direction(V)))) / 24,
    }
    ct = a * (x * x + y * y) / 2 + b * (x**4 + y**4) / 24
    ct4 = sp.expand(ct.subs(F, simultaneous=True)).coeff(t, 4)
    external = sp.Symbol("external", integer=True)
    v4, v6, v8, v10, v12 = sp.symbols("V4 V6 V8 V10 V12", integer=True)
    count = (v4, v6, v8, v10, v12)
    degree_sum = sum(n * v for n, v in zip((4, 6, 8, 10, 12), count))
    loop_formula = (
        1
        + sum((n - 2) * v for n, v in zip((4, 6, 8, 10, 12), count)) / 2
        - external / 2
    )
    return {
        "two_site_free_action": Q,
        "two_site_full_nonlocal_quartic_action": V,
        "all_generated_field_degrees": degrees,
        "quartic_transformed_counterterm": ct4,
        "one_loop_four_external_topologies_V4_V6_V8_V10_V12": topologies(4, 1),
        "loop_count_formula": loop_formula,
        "heavy_kernel_boundary": "Every h_ij denotes the full heavy inverse, not a finite derivative truncation. Light-only F leaves exact Gaussian heavy integration commutative. This algebraic two-site polarization check supplements the continuous functional proof; it is not a momentum cutoff.",
        "counterterm_boundary": "Transform all parent counterterms. A one-loop quadratic counterterm produces a quartic tree insertion; its generated sixth-order tadpole would first contribute at two loops. Higher derivative terms are perturbative insertions, never resummed as a new propagator.",
        "checks": {
            **{
                "full_generated_degree_" + str(n): sp.expand(degrees[n] - expected[n])
                for n in degrees
            },
            "all_generated_terms_retained": sp.expand(
                full.as_expr() - sum(t**n * v for n, v in degrees.items())
            ),
            "quadratic_counterterm_generates_quartic": sp.expand(
                ct4 - b * (x**4 + y**4) / 24 - a * (x * r + y * s)
            ),
            "Euler_and_half_edge_loop_identity": sp.expand(
                (degree_sum - external) / 2 - sum(count) + 1 - loop_formula
            ),
            "sixth_order_free_action_piece_kept": sp.expand(
                direction(direction(Q)) / 2
                - (k11 * r * r + 2 * k12 * r * s + k22 * s * s) / 2
            ),
        },
    }
