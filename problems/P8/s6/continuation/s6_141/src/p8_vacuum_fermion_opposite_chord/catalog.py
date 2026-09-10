"""Complete opposite-chord word ownership and proper-cycle UV degrees."""

from functools import cache
from itertools import combinations, permutations

import sympy as sp
from p8_vacuum_fermion_self_energy_chord import catalog as parent


def selected_words():
    return [r["word"] for r in parent.catalog() if r["cyclic_distance"] == 3]


def simple_cycles(word):
    ends = [i for i, x in enumerate(word) if x == "B"]
    edges = [(i, (i + 1) % 6, "F") for i in range(6)] + [(ends[0], ends[1], "B")]
    cycles = []
    for size in range(2, 8):
        for indices in combinations(range(7), size):
            degrees, neighbours = {}, {}
            for i in indices:
                u, v, _ = edges[i]
                for a, b in ((u, v), (v, u)):
                    degrees[a] = degrees.get(a, 0) + 1
                    neighbours.setdefault(a, set()).add(b)
            if any(d != 2 for d in degrees.values()):
                continue
            reached, pending = set(), [next(iter(degrees))]
            while pending:
                v = pending.pop()
                if v not in reached:
                    reached.add(v)
                    pending.extend(neighbours[v] - reached)
            if reached != set(degrees):
                continue
            nf = sum(edges[i][2] == "F" for i in indices)
            nb = size - nf
            cycles.append(
                {
                    "edges": indices,
                    "fermions": nf,
                    "bosons": nb,
                    "UV_degree": 4 - nf - 2 * nb,
                }
            )
    return cycles


@cache
def data():
    words = selected_words()
    selected = set(words)
    misses = 0
    for perm in permutations((1, 2, 3, 4)):
        relabel = dict(zip((1, 2, 3, 4), perm))
        misses += sum(
            parent.canonical(tuple(relabel.get(x, x) for x in w)) not in selected
            for w in words
        )
    cycles = [simple_cycles(w) for w in words]
    bad_cycles = sum(
        sorted((c["fermions"], c["bosons"], c["UV_degree"]) for c in cs)
        != [(3, 1, -1), (3, 1, -1), (6, 0, -2)]
        for cs in cycles
    )
    return {
        "selected_opposite_chord_words": words,
        "literal_proper_cycle_catalog": cycles,
        "word_normalization": "The n=6 determinant coefficient 1/6 and the Gaussian 1/2 times two B derivatives leave weight one for each cyclic word. Both loop orientations are retained.",
        "proper_subgraphs": "Exactly three simple cycles: two three-fermion/one-boson cycles of UV degree -1 and the six-fermion cycle of degree -2. No divergent proper 1PI subgraph.",
        "whole_graph": "Two loops, six fermion propagators and one boson: degree zero; only an overall quartic contact subtraction remains.",
        "scope": "Twelve opposite-endpoint words per sector, disjoint from the 24 self-energy and 24 vertex words.",
        "checks": {
            "twelve_selected_words": len(words) - 12,
            "no_duplicate_selected_word": len(selected) - 12,
            "S4_closed_subset": misses,
            "twelve_complete_proper_cycle_catalogs": len(cycles) - 12,
            "no_missing_or_extra_proper_cycle": bad_cycles,
            "sixth_determinant_all_word_factor": sp.binomial(6, 2) * sp.factorial(4) / 6
            - 60,
            "Gaussian_second_derivative_half_cancels": sp.Rational(1, 2) * 2 - 1,
            "opposite_endpoint_selections_times_external_orders": 2 * sp.factorial(3)
            - 12,
            "full_graph_two_loops": 7 - 6 + 1 - 2,
            "full_graph_logarithmic_degree": 4 * 2 - 6 - 2,
            "self_energy_and_opposite_disjoint": len(
                selected
                & {r["word"] for r in parent.catalog() if r["kind"] == "self_energy"}
            ),
            "remaining_twenty_four_vertex_words": 60 - 24 - len(words) - 24,
        },
    }
