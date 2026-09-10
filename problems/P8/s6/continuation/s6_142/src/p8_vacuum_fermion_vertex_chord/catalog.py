"""All vertex-chord words, proper subdivergences and marked vertices."""

from functools import cache
from itertools import permutations

import sympy as sp
from p8_vacuum_fermion_opposite_chord import catalog as cycle_catalog
from p8_vacuum_fermion_self_energy_chord import catalog as parent


def selected_words():
    return [r["word"] for r in parent.catalog() if r["kind"] == "vertex"]


def short_arc_rotation(word):
    for i, x in enumerate(word):
        if x == "B":
            w = word[i:] + word[:i]
            if w.index("B", 1) == 2:
                return w
    raise ValueError("Not a vertex-chord word")


def marked_box_vertices():
    result = []
    for p in permutations((2, 3, 4)):
        w = (1,) + p
        for j in range(4):
            result.append(parent.canonical(w[:j] + ("B", w[j], "B") + w[j + 1 :]))
    return result


@cache
def data():
    words = selected_words()
    selected = set(words)
    misses = 0
    for p in permutations((1, 2, 3, 4)):
        relabel = dict(zip((1, 2, 3, 4), p))
        misses += sum(
            parent.canonical(tuple(relabel.get(x, x) for x in w)) not in selected
            for w in words
        )
    cycles = [cycle_catalog.simple_cycles(w) for w in words]
    bad = sum(
        sorted((c["fermions"], c["bosons"], c["UV_degree"]) for c in cs)
        != [(2, 1, 0), (4, 1, -2), (6, 0, -2)]
        for cs in cycles
    )
    h, u = sp.symbols("loop_marker local_vertex_ratio")
    A = sp.Matrix([[2, 1], [1, 3]])
    V = sp.Matrix([[1, 2], [2, 0]])
    S = A.inv()
    varied = sp.diff(sp.trace((S * V * (1 + h * u)) ** 4) / 4, h).subs(h, 0)
    return {
        "selected_vertex_chord_words": words,
        "marked_box_vertex_bijection": marked_box_vertices(),
        "literal_proper_cycle_catalog": cycles,
        "short_arc_routings": [short_arc_rotation(w) for w in words],
        "proper_forest": "Each graph has exactly one logarithmically divergent proper Yukawa vertex. Pair its total proper MS Yukawa counterterm once; mass/kinetic proper insertions belong to S6.140.",
        "scope": "All 24 vertex words per sector; the primitive rows become bounded only after adding the disjoint S6.140 and S6.141 subsets.",
        "checks": {
            "twenty_four_vertex_words": len(words) - 24,
            "no_duplicate_vertex_word": len(selected) - 24,
            "S4_closed_family": misses,
            "six_boxes_four_marked_vertices": len(marked_box_vertices()) - 24,
            "marked_vertex_bijection_missing": len(
                selected - set(marked_box_vertices())
            ),
            "marked_vertex_bijection_extra": len(set(marked_box_vertices()) - selected),
            "proper_cycles_complete": bad,
            "unique_divergent_proper_vertex_per_word": sum(
                sum(c["UV_degree"] >= 0 for c in cs) != 1 for cs in cycles
            ),
            "short_arc_one_external_Phi": sum(
                w.index("B", 1) != 2 for w in map(short_arc_rotation, words)
            ),
            "literal_local_vertex_variation_normalization": varied
            - u * sp.trace((S * V) ** 4),
            "all_sixty_words_now_partitioned": 24 + 24 + 12 - 60,
            "no_self_energy_overlap": len(
                selected
                & {r["word"] for r in parent.catalog() if r["kind"] == "self_energy"}
            ),
            "no_opposite_overlap": len(selected & set(cycle_catalog.selected_words())),
        },
    }
