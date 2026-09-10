"""Independent spanning-forest and matrix-tree polynomials for every refinement."""

from functools import cache
from itertools import combinations

import sympy as sp

from . import graphs


def data(kind, choices):
    graphs.refinement(kind, choices)
    return _data(kind, choices)


@cache
def _data(kind, choices):
    g = graphs.refinement(kind, choices)
    vertices, edges = g["vertices"], g["edges"]
    a = sp.symbols("alpha0:" + str(len(edges)), positive=True)
    ext = g["external_labels"]

    def partition(selected):
        parent = {v: v for v in vertices}

        def root(v):
            while parent[v] != v:
                v = parent[v]
            return v

        for j in selected:
            u, v = edges[j]
            left, right = root(u), root(v)
            if left == right:
                return None
            parent[left] = right
        return tuple(
            tuple(v for v in vertices if root(v) == r)
            for r in sorted({root(v) for v in vertices})
        )

    U = sp.Integer(0)
    cuts = {key: sp.Integer(0) for key in ("single", "s", "t", "u")}
    tree_count = forest_count = 0
    for size in (len(vertices) - 1, len(vertices) - 2):
        for selected in combinations(range(len(edges)), size):
            parts = partition(selected)
            if parts is None or len(parts) != len(vertices) - size:
                continue
            monomial = sp.prod(v for j, v in enumerate(a) if j not in selected)
            if len(parts) == 1:
                U += monomial
                tree_count += 1
                continue
            forest_count += 1
            labels = set().union(*(set(ext[v]) for v in parts[0]))
            if len(labels) in (0, 4):
                continue
            if len(labels) in (1, 3):
                key = "single"
            elif labels in ({0, 1}, {2, 3}):
                key = "s"
            elif labels in ({0, 2}, {1, 3}):
                key = "t"
            else:
                key = "u"
            cuts[key] += monomial
    U = sp.expand(U)
    cuts = {k: sp.expand(v) for k, v in cuts.items()}
    # Independent polynomial Laplacian cofactor; only afterward substitute 1/alpha.
    weights = sp.symbols("edge_weight0:" + str(len(edges)))
    Lap = sp.zeros(len(vertices))
    lookup = {v: i for i, v in enumerate(vertices)}
    for w, (u, v) in zip(weights, edges):
        i, j = lookup[u], lookup[v]
        if i != j:
            Lap[i, i] += w
            Lap[j, j] += w
            Lap[i, j] -= w
            Lap[j, i] -= w
    Kirchhoff = sp.expand(Lap[:-1, :-1].det(method="domain-ge"))
    matrix_U = sp.expand(
        sp.prod(a) * Kirchhoff.subs({w: 1 / x for w, x in zip(weights, a)})
    )
    M = sp.Symbol("heavy_mass_squared_at_least_sixteen", positive=True)
    v = sp.Symbol("complex_forward_increment")
    total = sum(a)
    H = sum(a[4:])
    F = sp.expand(
        U * (sum(a[:4]) + M * H)
        - cuts["single"]
        - (2 + v) * cuts["s"]
        - (2 - v) * cuts["u"]
    )
    coarse = sp.expand(
        U * (sum(a[:4]) + 16 * H)
        - cuts["single"]
        - 3 * (cuts["s"] + cuts["t"] + cuts["u"])
    )
    pair = cuts["s"] + cuts["t"] + cuts["u"]
    derivative_enclosure = sp.Poly(sp.expand(U * total - pair), *a)
    return {
        "graph": g,
        "parameters": a,
        "U": U,
        "cut_polynomials": cuts,
        "F": F,
        "coarse_mass_sixteen_F": coarse,
        "H": H,
        "M": M,
        "v": v,
        "matrix_tree_U": matrix_U,
        "spanning_tree_count": tree_count,
        "spanning_two_forest_count": forest_count,
        "crossing_derivative_polynomial": cuts["u"] - cuts["s"],
        "derivative_majorant_remainder": derivative_enclosure.as_expr(),
        "derivative_majorant_coefficients_nonnegative": all(
            c >= 0 for _, c in derivative_enclosure.terms()
        ),
        "checks": {
            "independent_matrix_tree_matches_spanning_tree_sum": sp.expand(
                matrix_U - U
            ),
            "actual_forward_crossing_derivative": sp.expand(
                sp.diff(F, v) - cuts["u"] + cuts["s"]
            ),
            "U_has_loop_degree_two": sp.Poly(U, *a).total_degree() - 2,
            "F_has_loop_plus_one_degree_three": sp.Poly(F, *a).total_degree() - 3,
            "positive_mass_lift_exact": sp.expand(F - F.subs(M, 16) - (M - 16) * U * H),
            "coarse_forward_disc_kinematic_gap": sp.expand(
                F.subs(M, 16)
                - coarse
                - (1 - v) * cuts["s"]
                - 3 * cuts["t"]
                - (1 + v) * cuts["u"]
            ),
        },
    }
