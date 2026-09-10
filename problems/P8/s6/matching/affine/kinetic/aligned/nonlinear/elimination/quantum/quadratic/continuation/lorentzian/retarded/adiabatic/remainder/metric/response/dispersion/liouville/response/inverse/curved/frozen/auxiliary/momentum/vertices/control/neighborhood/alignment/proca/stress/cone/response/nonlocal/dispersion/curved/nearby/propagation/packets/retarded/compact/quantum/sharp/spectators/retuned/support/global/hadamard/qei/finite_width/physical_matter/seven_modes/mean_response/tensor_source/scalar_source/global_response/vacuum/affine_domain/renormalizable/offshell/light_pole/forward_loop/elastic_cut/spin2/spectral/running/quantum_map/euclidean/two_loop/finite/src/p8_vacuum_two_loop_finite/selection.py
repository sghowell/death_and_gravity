"""Finite-subsector selection, with exhaustive first-polynomial cross-check."""

from functools import cache
from itertools import combinations

import sympy as sp
from p8_vacuum_two_loop_denom import forests, graphs, subgraphs


def require_finite(kind, choices):
    g = graphs.refinement(kind, choices)
    if subgraphs.data(kind, choices)["UV_subgraphs"]:
        raise ValueError(
            "This refinement needs UV subtractions and is outside the finite subsector"
        )
    return g


@cache
def cases():
    return tuple(
        case for case in graphs.cases() if not subgraphs.data(*case)["UV_subgraphs"]
    )


def cotrees(kind, choices):
    require_finite(kind, choices)
    return _cotrees(kind, choices)


@cache
def _cotrees(kind, choices):
    g = graphs.refinement(kind, choices)
    edges = g["edges"]
    return tuple(
        pair
        for pair in combinations(range(len(edges)), 2)
        if graphs.connected(
            g["vertices"], tuple(e for j, e in enumerate(edges) if j not in pair)
        )
    )


@cache
def data():
    rows = []
    checks = {}
    for kind, choices in cases():
        d = forests.data(kind, choices)
        basis = cotrees(kind, choices)
        rebuilt = sum(sp.prod(d["parameters"][j] for j in pair) for pair in basis)
        label = kind + "_" + "".join(map(str, choices))
        checks["independent_cotree_" + label] = sp.expand(rebuilt - d["U"])
        rows.append(
            {
                "kind": kind,
                "choices": choices,
                "heavy_edge_count": d["graph"]["heavy_edge_count"],
                "cotrees": basis,
                "source_U_tree_count": d["spanning_tree_count"],
            }
        )
    return {
        "finite_refinement_count": len(rows),
        "excluded_subtraction_dependent_count": len(graphs.cases()) - len(rows),
        "finite_refinements": tuple(rows),
        "scope": "Precisely the UV-finite bare four-point refinements. Counterterm diagrams, divergent refinements and complete two-loop LSZ are not included.",
        "checks": checks,
    }
