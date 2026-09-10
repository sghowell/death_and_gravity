"""Literal six-vertex cyclic words and the adjacent self-energy subset."""

from functools import cache
from itertools import combinations, permutations

import sympy as sp


def canonical(word):
    i = word.index(1)
    return tuple(word[i:] + word[:i])


def catalog():
    rows = []
    for ends in combinations(range(1, 6), 2):
        for perm in permutations((2, 3, 4)):
            letters = iter(perm)
            word = (1,) + tuple(
                "B" if i in ends else next(letters) for i in range(1, 6)
            )
            distance = min(ends[1] - ends[0], 6 - ends[1] + ends[0])
            rows.append(
                {
                    "word": word,
                    "endpoints": ends,
                    "cyclic_distance": distance,
                    "kind": {1: "self_energy", 2: "vertex", 3: "finite_proper"}[
                        distance
                    ],
                }
            )
    return rows


def box_placements():
    return [
        canonical((1,) + perm[:gap] + ("B", "B") + perm[gap:])
        for perm in permutations((2, 3, 4))
        for gap in range(4)
    ]


@cache
def data():
    rows = catalog()
    selected = {r["word"] for r in rows if r["kind"] == "self_energy"}
    all_words = {r["word"] for r in rows}
    misses = 0
    phi, h = sp.symbols("background_phi insertion_h")
    A = sp.Matrix([[2, 1], [1, 3]])
    V = sp.Matrix([[1, 2], [2, 0]])
    Sigma = sp.Matrix([[3, 1], [1, 2]])
    S = A.inv()
    direct = -sp.trace((A + phi * V).inv() * Sigma)
    expected = -sp.trace((S * V) ** 4 * S * Sigma)
    shifted = S - h * S * Sigma * S
    varied = sp.diff(sp.trace((shifted * V) ** 4) / 4, h).subs(h, 0)
    for perm in permutations((1, 2, 3, 4)):
        mapping = dict(zip((1, 2, 3, 4), perm))
        for word in selected:
            mapped = canonical(tuple(mapping.get(x, x) for x in word))
            misses += mapped not in selected
    return {
        "all_primitive_words": rows,
        "selected_self_energy_words": sorted(selected, key=str),
        "first_box_marked_propagator_words": box_placements(),
        "sector_dictionary": "Scalar N Y^3; gauge N Cf a Y^2 with N=6,Cf=4/3. The proper self-energy kernel carries the internal coupling and open-line Casimir.",
        "scope": "Only 24 adjacent-endpoint words per scalar/gauge sector are bounded here. The 24 vertex and 12 other primitive words in each sector remain separate.",
        "checks": {
            "literal_determinant_insertion_fourth_coefficient": sp.diff(
                direct, phi, 4
            ).subs(phi, 0)
            / sp.factorial(4)
            - expected,
            "one_loop_box_covariance_variation_no_extra_half": varied - expected,
            "all_sixty_cyclic_words": len(rows) - 60,
            "no_duplicate_cyclic_words": len(all_words) - 60,
            "twenty_four_self_energy_words": len(selected) - 24,
            "twenty_four_vertex_words": sum(r["kind"] == "vertex" for r in rows) - 24,
            "twelve_remaining_primitive_words": sum(
                r["kind"] == "finite_proper" for r in rows
            )
            - 12,
            "six_boxes_four_internal_positions": len(box_placements()) - 6 * 4,
            "box_marking_bijection_missing": len(selected - set(box_placements())),
            "box_marking_bijection_extra": len(set(box_placements()) - selected),
            "self_energy_subset_closed_under_all_external_permutations": misses,
        },
    }
