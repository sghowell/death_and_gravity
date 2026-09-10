"""Actual local ultraviolet subgraphs and vertex-aware restricted forests."""

from functools import cache
from itertools import combinations

from . import graphs


def data(kind, choices):
    graphs.refinement(kind, choices)
    return _data(kind, choices)


@cache
def _data(kind, choices):
    g = graphs.refinement(kind, choices)
    edges = g["edges"]
    rows = []
    for size in range(1, len(edges) + 1):
        for selected in combinations(range(len(edges)), size):
            subedges = tuple(edges[j] for j in selected)
            vertices = tuple(sorted({v for edge in subedges for v in edge}))
            loops = size - len(vertices) + 1
            degree = 4 * loops - 2 * size
            if (
                loops < 1
                or degree < 0
                or not graphs.one_particle_irreducible(vertices, subedges)
            ):
                continue
            light_edges = sum(j < 4 for j in selected)
            heavy_edges = size - light_edges
            light_external = (
                sum(4 if g["vertex_valences"][v] == 4 else 2 for v in vertices)
                - 2 * light_edges
            )
            heavy_external = (
                sum(g["vertex_valences"][v] == 3 for v in vertices) - 2 * heavy_edges
            )
            rows.append(
                {
                    "edges": selected,
                    "vertices": vertices,
                    "loops": loops,
                    "superficial_momentum_degree": degree,
                    "light_external_legs": light_external,
                    "heavy_external_legs": heavy_external,
                    "is_whole_refined_graph": size == len(edges),
                }
            )
    forests = []
    for size in range(len(rows) + 1):
        for selected in combinations(range(len(rows)), size):
            compatible = True
            for i, j in combinations(selected, 2):
                a, b = rows[i], rows[j]
                disjoint = set(a["vertices"]).isdisjoint(b["vertices"])
                nested = set(a["edges"]) < set(b["edges"]) or set(b["edges"]) < set(
                    a["edges"]
                )
                if not (disjoint or nested):
                    compatible = False
                    break
            if compatible:
                forests.append(selected)
    overlap = []
    for i, j in combinations(range(len(rows)), 2):
        a, b = rows[i], rows[j]
        if set(a["edges"]).isdisjoint(b["edges"]) and not set(a["vertices"]).isdisjoint(
            b["vertices"]
        ):
            overlap.append((i, j))
    return {
        "kind": kind,
        "choices": choices,
        "UV_subgraphs": tuple(rows),
        "restricted_forest_index_sets": tuple(forests),
        "edge_disjoint_but_vertex_overlapping_pairs": tuple(overlap),
        "scope": "Combinatorial divergent subgraphs and compatible forest sets only. No subtraction integral, finite on-shell/MS conversion or two-loop LSZ coefficient has been evaluated. Counterterm operations on nested elements must be applied inside out.",
    }


@cache
def summary():
    rows = []
    for kind, choices in graphs.cases():
        d = data(kind, choices)
        rows.append(
            {
                "kind": kind,
                "choices": choices,
                "UV_subgraph_count": len(d["UV_subgraphs"]),
                "restricted_forest_count": len(d["restricted_forest_index_sets"]),
                "UV_external_types_and_degrees": tuple(
                    sorted(
                        {
                            (
                                r["light_external_legs"],
                                r["heavy_external_legs"],
                                r["superficial_momentum_degree"],
                                r["loops"],
                            )
                            for r in d["UV_subgraphs"]
                        }
                    )
                ),
                "edge_disjoint_vertex_overlap_count": len(
                    d["edge_disjoint_but_vertex_overlapping_pairs"]
                ),
            }
        )
    return tuple(rows)
