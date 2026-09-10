"""Exhaustive two-loop full-quartic light skeletons and heavy refinements."""

from functools import cache
from itertools import product
from math import factorial

import sympy as sp

KINDS = ("double_bubble", "wineglass", "tadpole_insertion")
PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def connected(vertices, edges):
    seen = {vertices[0]}
    for _ in vertices:
        for a, b in edges:
            if a in seen or b in seen:
                seen.update((a, b))
    return set(vertices) == seen


def one_particle_irreducible(vertices, edges):
    return connected(vertices, edges) and all(
        connected(vertices, edges[:j] + edges[j + 1 :])
        for j, (a, b) in enumerate(edges)
        if a != b
    )


@cache
def skeletons():
    result = []
    for counts in product(range(5), repeat=6):
        if sum(counts) != 4:
            continue
        edges = tuple(pair for pair, n in zip(PAIRS, counts) for _ in range(n))
        if not one_particle_irreducible((0, 1, 2), edges):
            continue
        degree = [0, 0, 0]
        for a, b in edges:
            degree[a] += 1
            degree[b] += 1
        if any(v > 4 for v in degree):
            continue
        external = tuple(4 - v for v in degree)
        if sum(external) != 4:
            continue
        kind = (
            "tadpole_insertion"
            if sum(counts[:3])
            else "wineglass"
            if sorted(external) == [1, 1, 2]
            else "double_bubble"
        )
        assignments = factorial(4) // sp.prod(factorial(n) for n in external)
        denom = (
            factorial(3)
            * sp.prod(factorial(n) for n in counts[3:])
            * sp.prod(2**n * factorial(n) for n in counts[:3])
        )
        result.append(
            {
                "kind": kind,
                "multiplicities": counts,
                "edges": edges,
                "external_counts": external,
                "external_label_assignments": int(assignments),
                "weight_per_external_assignment": sp.Rational(1, denom),
                "summed_labelled_vertex_Wick_weight": sp.Rational(assignments, denom),
            }
        )
    return tuple(result)


def refinement(kind, choices):
    if type(kind) is not str or kind not in KINDS:
        raise ValueError("Require one of the three actual two-loop skeleton families")
    if type(choices) is not tuple or len(choices) != 3:
        raise TypeError("Require a native three-vertex choice tuple")
    if any(type(n) is not int or not 0 <= n <= 3 for n in choices):
        raise ValueError("Every vertex is local zero or one of three heavy pairings")
    return _refinement(kind, choices)


@cache
def _refinement(kind, choices):
    row = next(r for r in skeletons() if r["kind"] == kind)
    edges = row["edges"]
    slots = {v: [] for v in range(3)}
    for j, (a, b) in enumerate(edges):
        slots[a].append(("edge", j, 0))
        slots[b].append(("edge", j, 1))
    label = 0
    for v, n in enumerate(row["external_counts"]):
        for _ in range(n):
            slots[v].append(("external", label, 0))
            label += 1
    attach = {}
    external = {}
    heavy = []
    valences = {}
    for v, choice in enumerate(choices):
        ordered = sorted(slots[v])
        if len(ordered) != 4:
            raise ValueError("The actual quartic half-edge count changed")
        left = ordered if choice == 0 else [ordered[0], ordered[choice]]
        right = [] if choice == 0 else [s for s in ordered if s not in left]
        for side, selected in enumerate((left, right)):
            if not selected:
                continue
            node = 2 * v + side
            external[node] = set()
            valences[node] = 4 if choice == 0 else 3
            for slot in selected:
                attach[(v, slot)] = node
                if slot[0] == "external":
                    external[node].add(slot[1])
        if choice:
            heavy.append((2 * v, 2 * v + 1))
    refined = tuple(
        (attach[(a, ("edge", j, 0))], attach[(b, ("edge", j, 1))])
        for j, (a, b) in enumerate(edges)
    ) + tuple(heavy)
    return {
        "kind": kind,
        "choices": choices,
        "vertices": tuple(sorted(external)),
        "edges": refined,
        "light_edge_count": 4,
        "heavy_edge_count": len(heavy),
        "external_labels": {v: tuple(sorted(s)) for v, s in external.items()},
        "vertex_valences": valences,
        "loop_count": len(refined) - len(external) + 1,
        "scope": "The heavy inverse is expanded only into its exact three exchange channels, never as a momentum derivative series. Light self-loop half-edges stay distinct during vertex splitting.",
    }


def cases():
    return tuple(
        (kind, choices) for kind in KINDS for choices in product(range(4), repeat=3)
    )


@cache
def data():
    rows = skeletons()
    phi, L, K = sp.symbols(
        "diagnostic_background diagnostic_quartic positive_free_kernel", positive=True
    )
    background_inverse = K + L * phi**2 / 2
    gamma2 = 3 * L / (24 * background_inverse**2) - 6 * (L * phi) ** 2 / (
        2 * 36 * background_inverse**3
    )
    four_point = sp.diff(gamma2, phi, 4).subs(phi, 0)
    return {
        "labelled_skeletons": rows,
        "labelled_skeleton_count": len(rows),
        "external_assignment_count": sum(r["external_label_assignments"] for r in rows),
        "full_heavy_refinement_count": len(cases()),
        "family_Wick_weights": {
            kind: sum(
                r["summed_labelled_vertex_Wick_weight"]
                for r in rows
                if r["kind"] == kind
            )
            for kind in KINDS
        },
        "independent_zero_dimensional_two_loop_effective_action": gamma2,
        "independent_zero_dimensional_fourth_derivative": four_point,
        "scope": "Positive combinatorial Wick weights and complete bare two-loop light-1PI topology, not signed renormalized amplitudes or a higher-loop error bound.",
        "checks": {
            "two_loop_four_external_requires_three_quartic_vertices": sp.Rational(
                2 * 2 + 4 - 2, 2
            )
            - 3,
            "four_internal_light_edges_from_half_edge_count": sp.Rational(3 * 4 - 4, 2)
            - 4,
            "total_external_assignment_count": sum(
                r["external_label_assignments"] for r in rows
            )
            - 72,
            "double_bubble_combinatorial_weight": sum(
                r["summed_labelled_vertex_Wick_weight"]
                for r in rows
                if r["kind"] == "double_bubble"
            )
            - sp.Rational(3, 4),
            "wineglass_combinatorial_weight": sum(
                r["summed_labelled_vertex_Wick_weight"]
                for r in rows
                if r["kind"] == "wineglass"
            )
            - 3,
            "tadpole_insertion_combinatorial_weight": sum(
                r["summed_labelled_vertex_Wick_weight"]
                for r in rows
                if r["kind"] == "tadpole_insertion"
            )
            - sp.Rational(3, 2),
            "independent_background_Gaussian_four_point_Wick_weight": sp.factor(
                four_point * K**4 / L**3
                - sum(r["summed_labelled_vertex_Wick_weight"] for r in rows)
            ),
        },
    }
