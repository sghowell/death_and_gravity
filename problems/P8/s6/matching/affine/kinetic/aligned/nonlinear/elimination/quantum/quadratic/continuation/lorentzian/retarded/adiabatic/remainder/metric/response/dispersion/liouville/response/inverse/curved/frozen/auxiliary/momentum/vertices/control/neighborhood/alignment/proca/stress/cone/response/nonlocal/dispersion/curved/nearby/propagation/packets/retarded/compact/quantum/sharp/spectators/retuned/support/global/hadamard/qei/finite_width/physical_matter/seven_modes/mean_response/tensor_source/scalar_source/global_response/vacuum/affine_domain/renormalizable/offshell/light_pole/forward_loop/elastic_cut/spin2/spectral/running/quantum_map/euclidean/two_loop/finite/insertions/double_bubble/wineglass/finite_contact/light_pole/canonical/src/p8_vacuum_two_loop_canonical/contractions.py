"""Independent marked-core contraction and counterterm insertion graph census."""

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations, pairwise, permutations, product

import sympy as sp
from p8_vacuum_two_loop_denom import graphs as parent
from p8_vacuum_two_loop_denom import subgraphs


def assignments(counts):
    rows = set()
    for labels in permutations(range(4)):
        pos = 0
        row = []
        for n in counts:
            row.append(tuple(sorted(labels[pos : pos + n])))
            pos += n
        rows.add(tuple(row))
    return tuple(sorted(rows))


def refine(row, external, choices):
    slots = {v: [] for v in range(len(choices))}
    for i, (u, v) in enumerate(row["edges"]):
        slots[u].append(("edge", i, 0))
        slots[v].append(("edge", i, 1))
    for v, labels in enumerate(external):
        slots[v].extend(("external", j, 0) for j in labels)
    attach, colors, ext, heavy = {}, {}, {}, []
    for v, c in enumerate(choices):
        ordered = sorted(slots[v])
        assert len(ordered) == 4
        left = ordered if c == 0 else [ordered[0], ordered[c]]
        right = [] if c == 0 else [j for j in ordered if j not in left]
        for side, selected in enumerate((left, right)):
            if not selected:
                continue
            n = 2 * v + side
            colors[n] = "L" if c == 0 else "G"
            ext[n] = tuple(
                sorted(slot[1] for slot in selected if slot[0] == "external")
            )
            for slot in selected:
                attach[(v, slot)] = n
        if c:
            heavy.append((2 * v, 2 * v + 1, "H"))
    edges = tuple(
        (attach[(u, ("edge", i, 0))], attach[(v, ("edge", i, 1))], "P")
        for i, (u, v) in enumerate(row["edges"])
    ) + tuple(heavy)
    return colors, ext, edges


def key(graph, powers):
    colors, ext, edges = graph
    groups = defaultdict(list)
    for v in colors:
        groups[(colors[v], ext[v])].append(v)
    group_names = sorted(groups)
    fixed_colors = tuple(color for color in group_names for _ in groups[color])
    least = None
    for blocks in product(*(permutations(groups[c]) for c in group_names)):
        ordered = tuple(v for block in blocks for v in block)
        lookup = {v: i for i, v in enumerate(ordered)}
        encoded = tuple(
            sorted(
                (min(lookup[u], lookup[v]), max(lookup[u], lookup[v]), t)
                for u, v, t in edges
            )
        )
        if least is None or encoded < least:
            least = encoded
    return powers, fixed_colors, least


def cores(graph):
    colors, _ext, edges = graph
    found = []
    for i, j in combinations(range(len(edges)), 2):
        a, b = edges[i], edges[j]
        if a[2] != "P" or b[2] != "P" or a[0] == a[1] or set(a[:2]) != set(b[:2]):
            continue
        vs = tuple(sorted(a[:2]))
        nL = sum(colors[v] == "L" for v in vs)
        found.append(((i, j), vs, ("C2", "C3", "C4")[nL]))
    return tuple(found)


def contract(graph, selected):
    colors, ext, edges = graph
    new_colors, new_ext = dict(colors), dict(ext)
    mapping = {v: v for v in colors}
    removed = set()
    for n, (es, vs, color) in enumerate(selected, max(colors) + 1):
        new_colors[n] = color
        new_ext[n] = tuple(sorted(label for v in vs for label in ext[v]))
        for v in vs:
            mapping[v] = n
            del new_colors[v]
            del new_ext[v]
        removed.update(es)
    return (
        new_colors,
        new_ext,
        tuple(
            (mapping[u], mapping[v], t)
            for i, (u, v, t) in enumerate(edges)
            if i not in removed
        ),
    )


def mark_vertex(graph, v, color):
    colors, ext, edges = graph
    changed = dict(colors)
    changed[v] = color
    return changed, dict(ext), edges


def insert_mass(graph, edge, count=1):
    colors, ext, edges = graph
    u, v, kind = edges[edge]
    assert kind == "H"
    new_colors, new_ext = dict(colors), dict(ext)
    ns = tuple(range(max(colors) + 1, max(colors) + 1 + count))
    for n in ns:
        new_colors[n] = "C2"
        new_ext[n] = ()
    chain = (u,) + ns + (v,)
    new_edges = (
        edges[:edge]
        + edges[edge + 1 :]
        + tuple((a, b, "H") for a, b in pairwise(chain))
    )
    return new_colors, new_ext, new_edges


@cache
def data():
    left, left_pairs = defaultdict(F), defaultdict(F)
    left_signed, pair_signed = defaultdict(F), defaultdict(F)
    checks = {}
    by_family = Counter()
    record_count = single_count = pair_count = 0
    hist = Counter()
    assert sum(len(assignments(r["external_counts"])) for r in parent.skeletons()) == 72
    for row in parent.skeletons():
        w = row["weight_per_external_assignment"]
        weight = F(int(w.p), int(w.q))
        for ex in assignments(row["external_counts"]):
            for choices in product(range(4), repeat=3):
                graph = refine(row, ex, choices)
                powers = (
                    sum(c == "L" for c in graph[0].values()),
                    sum(c == "G" for c in graph[0].values()),
                )
                selected = cores(graph)
                for core in selected:
                    left[key(contract(graph, (core,)), powers)] += weight
                    left_signed[key(contract(graph, (core,)), powers)] += (-1) ** (
                        powers[0] + 1
                    ) * weight
                    single_count += 1
                    hist[core[2]] += 1
                for a, b in combinations(selected, 2):
                    if set(a[1]).isdisjoint(b[1]):
                        left_pairs[key(contract(graph, (a, b)), powers)] += weight
                        pair_signed[key(contract(graph, (a, b)), powers)] += (
                            -1
                        ) ** powers[0] * weight
                        pair_count += 1
                record_count += 1
                by_family[row["kind"]] += 1

    for kind, choices in parent.cases():
        row = next(r for r in parent.skeletons() if r["kind"] == kind)
        graph = refine(row, assignments(row["external_counts"])[0], choices)
        expected = {
            r["edges"]
            for r in subgraphs.data(kind, choices)["UV_subgraphs"]
            if r["loops"] == 1
            and (r["light_external_legs"], r["heavy_external_legs"])
            in ((4, 0), (2, 1), (0, 2))
        }
        checks[
            "same_frozen_proper_UV_core_selection_"
            + kind
            + "_"
            + "".join(map(str, choices))
        ] = int({r[0] for r in cores(graph)} == expected) - 1

    right = defaultdict(F)
    right_signed = defaultdict(F)
    row = {"edges": ((0, 1), (0, 1)), "external_counts": (2, 2)}
    for ex in assignments((2, 2)):
        for choices in product(range(4), repeat=2):
            graph = refine(row, ex, choices)
            nL, nG = (
                sum(c == "L" for c in graph[0].values()),
                sum(c == "G" for c in graph[0].values()),
            )
            for v, color in graph[0].items():
                coefficient = F(3, 2) if color == "L" else F(1, 2)
                right[
                    key(
                        mark_vertex(graph, v, "C4" if color == "L" else "C3"),
                        (nL + 1, nG),
                    )
                ] += F(1, 4) * coefficient
                right_signed[
                    key(
                        mark_vertex(graph, v, "C4" if color == "L" else "C3"),
                        (nL + 1, nG),
                    )
                ] += (-1) ** nL * F(1, 4) * coefficient
            for j, (_, _, kind) in enumerate(graph[2]):
                if kind == "H":
                    right[key(insert_mass(graph, j), (nL, nG + 2))] += F(1, 8)
                    right_signed[key(insert_mass(graph, j), (nL, nG + 2))] += (-1) ** (
                        nL + 1
                    ) * F(1, 8)
    differences = [
        (k, left[k], right[k]) for k in set(left) | set(right) if left[k] != right[k]
    ]

    right_pairs = defaultdict(F)
    right_pairs_signed = defaultdict(F)
    for ex in assignments((2, 2)):
        tree = ({0: "G", 1: "G"}, {0: ex[0], 1: ex[1]}, ((0, 1, "H"),))
        right_pairs[key(mark_vertex(mark_vertex(tree, 0, "C3"), 1, "C3"), (2, 2))] += F(
            1, 8
        )
        right_pairs_signed[
            key(mark_vertex(mark_vertex(tree, 0, "C3"), 1, "C3"), (2, 2))
        ] += F(1, 8)
        for v in (0, 1):
            right_pairs[key(insert_mass(mark_vertex(tree, v, "C3"), 0), (1, 4))] += F(
                1, 8
            )
            right_pairs_signed[
                key(insert_mass(mark_vertex(tree, v, "C3"), 0), (1, 4))
            ] -= F(1, 8)
        right_pairs[key(insert_mass(tree, 0, 2), (0, 6))] += F(1, 8)
        right_pairs_signed[key(insert_mass(tree, 0, 2), (0, 6))] += F(1, 8)
    differences2 = [
        (k, left_pairs[k], right_pairs[k])
        for k in set(left_pairs) | set(right_pairs)
        if left_pairs[k] != right_pairs[k]
    ]
    assert not differences and not differences2

    rows = []
    pair_rows = []
    for prefix, first, second, first_signed, second_signed, target in (
        ("single", left, right, left_signed, right_signed, rows),
        (
            "disjoint",
            left_pairs,
            right_pairs,
            pair_signed,
            right_pairs_signed,
            pair_rows,
        ),
    ):
        for i, k in enumerate(sorted(set(first) | set(second))):
            checks[prefix + "_absolute_Wick_weight_" + str(i)] = sp.Rational(
                first[k] - second[k]
            )
            checks[prefix + "_signed_Feynman_weight_" + str(i)] = sp.Rational(
                first_signed[k] - second_signed[k]
            )
            target.append(
                {
                    "index": i,
                    "coupling_powers_L_G": k[0],
                    "vertex_types_and_external_labels": k[1],
                    "typed_cograph_edges": k[2],
                    "contracted_raw_Wick_weight": sp.Rational(first[k]),
                    "independent_counterterm_Wick_weight": sp.Rational(second[k]),
                    "contracted_signed_coefficient": sp.Rational(first_signed[k]),
                    "independent_signed_coefficient": sp.Rational(second_signed[k]),
                }
            )
            checks[prefix + "_cograph_loop_count_" + str(i)] = (
                len(k[2]) - len(k[1]) + 1 - (1 if prefix == "single" else 0)
            )
    checks.update(
        {
            "all_labelled_two_loop_records": record_count - 4608,
            "all_marked_proper_vertex_or_heavy_mass_cores": single_count - 1152,
            "all_disjoint_proper_core_pairs": pair_count - 72,
            "all_one_core_cograph_classes": len(left) - 144,
            "all_disjoint_cograph_classes": len(left_pairs) - 12,
            "Phi4_core_count": hist["C4"] - 288,
            "HPhi2_core_count": hist["C3"] - 576,
            "H2_core_count": hist["C2"] - 288,
            "double_bubble_records": by_family["double_bubble"] - 1152,
            "wineglass_records": by_family["wineglass"] - 2304,
            "tadpole_records": by_family["tadpole_insertion"] - 1152,
        }
    )
    return {
        "raw_external_labelled_graph_records": record_count,
        "single_proper_core_contractions": single_count,
        "disjoint_proper_core_pairs": pair_count,
        "core_type_counts": dict(hist),
        "raw_records_by_family": dict(by_family),
        "all_single_core_cograph_weight_matches": rows,
        "all_disjoint_tree_cograph_weight_matches": pair_rows,
        "scope": "Exact full external-leg labels, local/cubic/counterterm vertex colors, light/heavy propagator types and coupling powers are preserved. Every marked contraction is matched to an independently generated parameter-counterterm insertion graph, including signs. One factor I0/(16pi^2) per contracted logarithmic light core is suppressed in the tables. Light quadratic and tadpole cores use the separate fixed on-shell insertion calculation.",
        "checks": checks,
    }
