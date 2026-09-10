"""Exhaustive disjoint selection from the immutable full heavy graph inventory."""

from functools import cache

import sympy as sp
from p8_vacuum_two_loop_denom import graphs, subgraphs


def predicate(choices):
    return choices[2] in (0, 1) and not (choices[0] in (2, 3) and choices[1] in (2, 3))


def cases():
    return tuple(
        choices
        for kind, choices in graphs.cases()
        if kind == "double_bubble" and predicate(choices)
    )


def require_group(choices):
    g = graphs.refinement("double_bubble", choices)
    if not predicate(choices):
        raise ValueError("Require a subtraction-dependent factorizing double bubble")
    return g


@cache
def data():
    checks, rows = {}, []
    finite_overlap = []
    for kind, choices in graphs.cases():
        if kind != "double_bubble":
            continue
        d = subgraphs.data(kind, choices)
        key = "".join(map(str, choices))
        checks["UV_selection_" + key] = sp.Integer(bool(d["UV_subgraphs"])) - int(
            predicate(choices)
        )
        if choices[2] in (0, 1) and not predicate(choices):
            finite_overlap.append(choices)
        if not predicate(choices):
            continue
        count = len(d["UV_subgraphs"])
        mode = {1: "single", 2: "disjoint", 3: "overlap"}[count]
        rows.append(
            {
                "choices": choices,
                "mode": mode,
                "UV_subgraphs": d["UV_subgraphs"],
                "restricted_forest_index_sets": d["restricted_forest_index_sets"],
            }
        )
        checks["all_UV_degrees_logarithmic_" + key] = sum(
            abs(r["superficial_momentum_degree"]) for r in d["UV_subgraphs"]
        )
    return {
        "selected_rows": tuple(rows),
        "selected_refinement_count": len(rows),
        "forest_class_counts": {
            mode: sum(r["mode"] == mode for r in rows)
            for mode in ("single", "disjoint", "overlap")
        },
        "already_finite_separable_refinements_excluded": tuple(finite_overlap),
        "already_finite_separable_count": len(finite_overlap),
        "remaining_wineglass_refinements": tuple(
            choices
            for kind, choices in graphs.cases()
            if kind == "wineglass" and subgraphs.data(kind, choices)["UV_subgraphs"]
        ),
        "local_UV_external_types": tuple(
            sorted(
                {
                    (s["light_external_legs"], s["heavy_external_legs"])
                    for r in rows
                    for s in r["UV_subgraphs"]
                }
            )
        ),
        "scope": "Exactly 24 subtraction-dependent bare refinements. The eight finite separable refinements and 32 nonfactorizing finite double bubbles are already in S6.121 and are not counted again.",
        "checks": checks,
    }
