"""Exhaustive two-loop two-point skeletons, heavy refinements and UV forests."""

from functools import cache
from itertools import combinations, product
from math import factorial, prod

import sympy as sp
from p8_vacuum_two_loop_denom.graphs import one_particle_irreducible

PAIRS = ((0, 0), (1, 1), (0, 1))


@cache
def skeletons():
    rows = []
    for counts in product(range(4), repeat=3):
        if sum(counts) != 3:
            continue
        edges = tuple(e for e, n in zip(PAIRS, counts) for _ in range(n))
        if not one_particle_irreducible((0, 1), edges):
            continue
        degrees = [0, 0]
        for a, b in edges:
            degrees[a] += 1
            degrees[b] += 1
        if max(degrees) > 4:
            continue
        ext = tuple(4 - d for d in degrees)
        assignments = factorial(2) // prod(factorial(n) for n in ext)
        denominator = (
            factorial(2)
            * factorial(counts[2])
            * prod(2**n * factorial(n) for n in counts[:2])
        )
        rows.append(
            {
                "kind": "sunset" if counts[2] == 3 else "nested_tadpole",
                "edges": edges,
                "counts": counts,
                "external_counts": ext,
                "assignments": assignments,
                "weight": sp.Rational(assignments, denominator),
            }
        )
    return rows


def refine(kind, choices):
    if type(kind) is not str or kind not in ("sunset", "nested_tadpole"):
        raise ValueError("Require an actual two-point skeleton family")
    if type(choices) is not tuple or len(choices) != 2:
        raise TypeError("Require a native two-vertex choice tuple")
    if any(type(c) is not int or not 0 <= c <= 3 for c in choices):
        raise ValueError("Require local zero or one of three heavy pairings")
    return _refine(kind, choices)


@cache
def _refine(kind, choices):
    row = next(r for r in skeletons() if r["kind"] == kind)
    slots = {v: [] for v in range(2)}
    for j, (a, b) in enumerate(row["edges"]):
        slots[a].append(("edge", j, 0))
        slots[b].append(("edge", j, 1))
    label = 0
    for v, n in enumerate(row["external_counts"]):
        for _ in range(n):
            slots[v].append(("external", label, 0))
            label += 1
    attach, ext, valences, heavy = {}, {}, {}, []
    for v, choice in enumerate(choices):
        ordered = sorted(slots[v])
        left = ordered if choice == 0 else [ordered[0], ordered[choice]]
        right = [] if choice == 0 else [s for s in ordered if s not in left]
        for side, selected in enumerate((left, right)):
            if not selected:
                continue
            node = 2 * v + side
            ext[node] = []
            valences[node] = 4 if choice == 0 else 3
            for slot in selected:
                attach[(v, slot)] = node
                if slot[0] == "external":
                    ext[node].append(slot[1])
        if choice:
            heavy.append((2 * v, 2 * v + 1))
    edges = tuple(
        (attach[(a, ("edge", j, 0))], attach[(b, ("edge", j, 1))])
        for j, (a, b) in enumerate(row["edges"])
    ) + tuple(heavy)
    return {
        "kind": kind,
        "choices": choices,
        "vertices": tuple(sorted(ext)),
        "edges": edges,
        "external": ext,
        "valences": valences,
    }


def uv(g):
    rows = []
    for n in range(1, len(g["edges"]) + 1):
        for selected in combinations(range(len(g["edges"])), n):
            edges = tuple(g["edges"][i] for i in selected)
            vertices = tuple(sorted({v for e in edges for v in e}))
            loops = n - len(vertices) + 1
            degree = 4 * loops - 2 * n
            if loops < 1 or degree < 0 or not one_particle_irreducible(vertices, edges):
                continue
            light = sum(i < 3 for i in selected)
            h = n - light
            ephi = sum(4 if g["valences"][v] == 4 else 2 for v in vertices) - 2 * light
            eh = sum(g["valences"][v] == 3 for v in vertices) - 2 * h
            rows.append((selected, loops, degree, ephi, eh))
    return rows


def paths(g, a, b):
    result = []

    def visit(v, seen, path):
        if v == b:
            result.append(path)
            return
        for j, (u, w) in enumerate(g["edges"]):
            if u == v and w not in seen:
                visit(w, seen | {w}, path + ((j, 1),))
            if w == v and u not in seen:
                visit(u, seen | {u}, path + ((j, -1),))

    visit(a, {a}, ())
    return result


def cases():
    return tuple(
        (kind, c)
        for kind in ("sunset", "nested_tadpole")
        for c in product(range(4), repeat=2)
    )


def forests(g):
    rows = uv(g)
    vertices = [{v for i in r[0] for v in g["edges"][i]} for r in rows]
    result = []
    for n in range(len(rows) + 1):
        for selected in combinations(range(len(rows)), n):
            if all(
                set(rows[i][0]) <= set(rows[j][0])
                or set(rows[j][0]) <= set(rows[i][0])
                or vertices[i].isdisjoint(vertices[j])
                for i, j in combinations(selected, 2)
            ):
                result.append(selected)
    return tuple(result)


def group(kind, choices):
    a, b = choices
    if kind == "sunset":
        return (
            "local"
            if a == b == 0
            else "mixed"
            if not a or not b
            else "same_heavy"
            if a == b
            else "finite"
        )
    return (
        "outer_constant"
        if a in (0, 1)
        else "inner_constant"
        if b in (0, 1)
        else "inner_mixed"
    )


@cache
def data():
    rows = skeletons()
    phi, L, K = sp.symbols("diagnostic_phi diagnostic_L positive_K", positive=True)
    background = K + L * phi**2 / 2
    gamma2 = L / (8 * background**2) - (L * phi) ** 2 / (12 * background**3)
    derivative = sp.diff(gamma2, phi, 2).subs(phi, 0)
    summaries = []
    checks = {
        "two_point_requires_two_quartic_vertices": sp.Rational(2 * 2 + 2 - 2, 2) - 2,
        "three_internal_light_edges": sp.Rational(2 * 4 - 2, 2) - 3,
        "three_labelled_skeletons": len(rows) - 3,
        "sunset_Wick_weight": sum(r["weight"] for r in rows if r["kind"] == "sunset")
        - sp.Rational(1, 6),
        "nested_Wick_weight": sum(
            r["weight"] for r in rows if r["kind"] == "nested_tadpole"
        )
        - sp.Rational(1, 4),
        "independent_background_second_derivative": sp.factor(
            -derivative * K**3 / L**2 - sum(r["weight"] for r in rows)
        ),
    }
    expected = {
        "sunset_local": 1,
        "sunset_mixed": 6,
        "sunset_same_heavy": 3,
        "sunset_finite": 6,
        "nested_tadpole_outer_constant": 8,
        "nested_tadpole_inner_constant": 4,
        "nested_tadpole_inner_mixed": 4,
    }
    counts = {k: 0 for k in expected}
    for kind, choice in cases():
        g = refine(kind, choice)
        key = kind + "_" + str(choice[0]) + str(choice[1])
        label = kind + "_" + group(kind, choice)
        counts[label] += 1
        cores = uv(g)
        summaries.append(
            {
                "graph": g,
                "group": label,
                "UV_cores": cores,
                "restricted_forests": forests(g),
            }
        )
        checks[key + "_loop_number"] = len(g["edges"]) - len(g["vertices"]) + 1 - 2
        checks[key + "_all_local_UV_types"] = (
            int(
                all(
                    (r[3], r[4]) in ((2, 0), (0, 1), (0, 2), (2, 1), (4, 0))
                    for r in cores
                )
            )
            - 1
        )
        proper = [r for r in cores if r[1] == 1]
        if kind == "sunset":
            target = (
                (3, 4, 0)
                if choice == (0, 0)
                else (
                    (1, 2, 1)
                    if 0 in choice
                    else (1, 0, 2)
                    if choice[0] == choice[1]
                    else (0, 0, 0)
                )
            )
            checks[key + "_proper_core_count"] = len(proper) - target[0]
            checks[key + "_proper_core_external_type"] = (
                int(all((r[3], r[4]) == target[1:] for r in proper)) - 1
            )
        if label == "nested_tadpole_inner_mixed":
            checks[key + "_sole_inner_mixed_core"] = (
                int(len(cores) == 1 and cores[0][1:] == (1, 0, 2, 0)) - 1
            )
    checks.update({"partition_" + k: counts[k] - v for k, v in expected.items()})
    return {
        "labelled_skeletons": rows,
        "raw_refinement_count": len(cases()),
        "UV_and_forest_inventory": summaries,
        "group_counts": counts,
        "independent_zero_dimensional_effective_action": gamma2,
        "independent_second_derivative": derivative,
        "scope": "Exactly the 32 reduced Gaussian-heavy two-loop light two-point refinements. UV cores and forests are retained even when the final outer on-shell operation annihilates their momentum-independent contribution.",
        "checks": checks,
    }


def bad_cases():
    return (
        [
            ("bad_family_" + str(i), refine, (v, (0, 0)))
            for i, v in enumerate((True, 0, None, "", "wineglass"))
        ]
        + [
            ("bad_choice_tuple_" + str(i), refine, ("sunset", v))
            for i, v in enumerate((None, [0, 0], (0,), (0, 0, 0)))
        ]
        + [
            ("bad_vertex_" + str(i), refine, ("sunset", (0, v)))
            for i, v in enumerate((True, False, 1.0, sp.Integer(1), -1, 4, None, "1"))
        ]
    )
