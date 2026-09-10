"""Exact remaining sixteen refinements and their logarithmic local cores."""

from functools import cache
from itertools import permutations

import sympy as sp
from p8_vacuum_two_loop_denom import graphs, subgraphs


def cases():
    return tuple(
        c
        for kind, c in graphs.cases()
        if kind == "wineglass" and c[1] in (0, 3) and c[2] in (0, 3)
    )


def require_group(choices):
    g = graphs.refinement("wineglass", choices)
    if choices not in cases():
        raise ValueError(
            "Require one of the sixteen subtraction-dependent wineglass refinements"
        )
    return g


@cache
def data():
    checks, rows = {}, []
    for kind, c in graphs.cases():
        if kind != "wineglass":
            continue
        d = subgraphs.data(kind, c)
        key = "".join(map(str, c))
        checks["exact_remaining_selection_" + key] = sp.Integer(
            bool(d["UV_subgraphs"])
        ) - int(c in cases())
        if c not in cases():
            continue
        checks["all_local_cores_logarithmic_" + key] = sum(
            abs(r["superficial_momentum_degree"]) for r in d["UV_subgraphs"]
        )
        overall = c[0] in (0, 1) and c[1:] == (0, 0)
        checks["exact_overall_core_selection_" + key] = (
            len(d["UV_subgraphs"]) - 1 - int(overall)
        )
        rows.append(
            {
                "choices": c,
                "has_overall_core": overall,
                "UV_subgraphs": d["UV_subgraphs"],
                "restricted_forest_index_sets": d["restricted_forest_index_sets"],
            }
        )
    weight = sum(
        r["summed_labelled_vertex_Wick_weight"]
        for r in graphs.skeletons()
        if r["kind"] == "wineglass"
    )
    checks["same_total_external_assignment_Wick_weight"] = weight - 3
    channel_weights = {name: sp.Integer(0) for name in ("s", "t", "u")}
    for skeleton in graphs.skeletons():
        if skeleton["kind"] != "wineglass":
            continue
        counts = skeleton["external_counts"]
        assignments = set()
        for labels in permutations(range(4)):
            offset, assignment = 0, []
            for count in counts:
                assignment.append(tuple(sorted(labels[offset : offset + count])))
                offset += count
            assignments.add(tuple(assignment))
        for assignment in assignments:
            pair = assignment[counts.index(2)]
            channel = (
                "s"
                if pair in ((0, 1), (2, 3))
                else "t"
                if pair in ((0, 2), (1, 3))
                else "u"
            )
            channel_weights[channel] += skeleton["weight_per_external_assignment"]
    for channel, value in channel_weights.items():
        checks["all_labelled_" + channel + "_channel_Wick_weight"] = value - 1
    return {
        "selected_rows": tuple(rows),
        "selected_refinement_count": len(rows),
        "overall_subtraction_refinement_count": sum(
            r["has_overall_core"] for r in rows
        ),
        "finite_after_inner_subtraction_count": sum(
            not r["has_overall_core"] for r in rows
        ),
        "total_family_Wick_weight": weight,
        "explicit_external_assignment_channel_weights": channel_weights,
        "one_loop_local_external_types": tuple(
            sorted(
                {
                    (s["light_external_legs"], s["heavy_external_legs"])
                    for r in rows
                    for s in r["UV_subgraphs"]
                    if s["loops"] == 1
                }
            )
        ),
        "two_loop_local_external_types": tuple(
            sorted(
                {
                    (s["light_external_legs"], s["heavy_external_legs"])
                    for r in rows
                    for s in r["UV_subgraphs"]
                    if s["loops"] == 2
                }
            )
        ),
        "scope": "Exactly the sixteen raw refinements excluded from the preceding three disjoint groups. Fourteen become integrable after their one-loop core subtraction; two also require an overall logarithmic local subtraction.",
        "checks": checks,
    }
