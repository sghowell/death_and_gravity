"""Three cyclic words, all proper cycles and overlapping forest ownership."""

from functools import cache
from itertools import combinations

import sympy as sp


def word_rows():
    rows = []
    for pair in combinations((1, 2, 3), 2):
        word = ["P1", None, None, None]
        for j in pair:
            word[j] = "B"
        word[next(j for j in (1, 2, 3) if word[j] is None)] = "P2"
        distance = pair[1] - pair[0]
        arcs = (distance, 4 - distance)
        rows.append(
            {
                "word": tuple(word),
                "boson_vertices": pair,
                "fermion_arcs": arcs,
                "kind": "vertex" if arcs == (2, 2) else "self_energy",
                "proper_cycle_degrees": (0, 2 - arcs[0], 2 - arcs[1]),
            }
        )
    return rows


def graph_cycles(row):
    a, b = row["boson_vertices"]
    edges = [(i, (i + 1) % 4, "f") for i in range(4)] + [(a, b, "b")]
    cycles = []
    for n in range(2, 6):
        for ids in combinations(range(5), n):
            degree = {}
            for i in ids:
                u, v, _ = edges[i]
                degree[u] = degree.get(u, 0) + 1
                degree[v] = degree.get(v, 0) + 1
            if not all(d == 2 for d in degree.values()):
                continue
            seen = {next(iter(degree))}
            while True:
                new = seen | {
                    v
                    for i in ids
                    for u, v in (edges[i][:2], edges[i][:2][::-1])
                    if u in seen
                }
                if new == seen:
                    break
                seen = new
            if seen != set(degree):
                continue
            ferm = sum(edges[i][2] == "f" for i in ids)
            bos = n - ferm
            cycles.append(
                {
                    "edges": ids,
                    "fermion_lines": ferm,
                    "boson_lines": bos,
                    "UV_degree": 4 - ferm - 2 * bos,
                    "local_after_chord_closure": bos == 0,
                }
            )
    return cycles


@cache
def data():
    rows = word_rows()
    checks = {
        "cyclic_determinant_and_contraction_factor": sp.factorial(3) / 2 - 3,
        "two_self_energy_words": sum(r["kind"] == "self_energy" for r in rows) - 2,
        "one_vertex_word": sum(r["kind"] == "vertex" for r in rows) - 1,
        "overall_quadratic_UV_degree": 8 - 4 - 2 - 2,
    }
    for i, row in enumerate(rows):
        cycles = graph_cycles(row)
        checks[f"three_proper_cycles_{i}"] = len(cycles) - 3
        checks[f"one_local_fermion_cycle_{i}"] = (
            sum(c["local_after_chord_closure"] for c in cycles) - 1
        )
        for j, (a, b) in enumerate(combinations(cycles, 2)):
            checks[f"proper_cycles_overlap_{i}_{j}"] = (
                int(bool(set(a["edges"]) & set(b["edges"]))) - 1
            )
        expected = [-1, 0, 1] if row["kind"] == "self_energy" else [0, 0, 0]
        checks[f"proper_degrees_{i}"] = sp.Matrix(
            sorted(c["UV_degree"] for c in cycles)
        ) - sp.Matrix(expected)
    return {
        "words": rows,
        "all_proper_cycles": [graph_cycles(r) for r in rows],
        "scope": "Three words per scalar/gauge quadratic primitive. The whole four-fermion cycle gives a local Phi^2 term after chord closure; its degree-zero subtraction vanishes under the nonlocal on-shell projection. All pairs of proper cycles overlap, so no product of their counterterms is a forest.",
        "checks": checks,
    }
