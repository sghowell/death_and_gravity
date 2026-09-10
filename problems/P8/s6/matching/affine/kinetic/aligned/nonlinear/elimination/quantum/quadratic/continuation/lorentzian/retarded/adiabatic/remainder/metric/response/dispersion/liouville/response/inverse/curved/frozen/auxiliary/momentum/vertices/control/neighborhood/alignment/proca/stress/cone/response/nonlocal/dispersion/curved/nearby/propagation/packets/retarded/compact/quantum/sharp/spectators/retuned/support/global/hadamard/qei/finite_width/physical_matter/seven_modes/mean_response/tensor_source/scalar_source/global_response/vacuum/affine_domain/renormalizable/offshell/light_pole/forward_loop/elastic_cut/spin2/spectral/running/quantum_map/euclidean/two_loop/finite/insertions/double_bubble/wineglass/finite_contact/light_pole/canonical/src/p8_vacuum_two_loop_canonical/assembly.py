"""Exhaustive disjoint assembly and perturbative-order ownership of every term."""

from functools import cache
from itertools import combinations

from p8_vacuum_two_loop_denom import graphs, subgraphs
from p8_vacuum_two_loop_double_bubble import selection as double
from p8_vacuum_two_loop_finite import selection as finite
from p8_vacuum_two_loop_insertions import subtraction as inserted
from p8_vacuum_two_loop_wineglass import selection as wine


@cache
def data():
    raw = set(graphs.cases())
    groups = {
        "finite": set(finite.cases()),
        "insertion": {c for c in raw if c[0] == "tadpole_insertion"},
        "double_bubble": {("double_bubble", c) for c in double.cases()},
        "wineglass": {("wineglass", c) for c in wine.cases()},
    }
    checks = {
        "complete_bare_two_loop_refinement_union": len(
            raw.symmetric_difference(set().union(*groups.values()))
        ),
        "finite_group_count": len(groups["finite"]) - 88,
        "insertion_group_count": len(groups["insertion"]) - 64,
        "subtracted_double_bubble_group_count": len(groups["double_bubble"]) - 24,
        "subtracted_wineglass_group_count": len(groups["wineglass"]) - 16,
        "same_actual_inner_OS_insertion_group_count": len(groups["insertion"])
        - inserted.data()["grouped_refinement_count"],
    }
    for a, b in combinations(groups, 2):
        checks[a + "_" + b + "_disjoint"] = len(groups[a] & groups[b])
    ownership = []
    two_loop_types = set()
    for kind, c in sorted(raw):
        owners = [name for name, selected in groups.items() if (kind, c) in selected]
        checks["exactly_one_owner_" + kind + "_" + "".join(map(str, c))] = (
            len(owners) - 1
        )
        ownership.append({"kind": kind, "choices": c, "integrated_group": owners[0]})
        for core in subgraphs.data(kind, c)["UV_subgraphs"]:
            if core["loops"] == 2:
                two_loop_types.add(
                    (
                        core["light_external_legs"],
                        core["heavy_external_legs"],
                        core["superficial_momentum_degree"],
                    )
                )
    checks["all_two_loop_proper_or_overall_UV_types_are_local_reference_types"] = (
        int(two_loop_types == {(4, 0, 0), (2, 1, 0), (0, 2, 0)}) - 1
    )
    return {
        "exact_raw_refinement_ownership": ownership,
        "raw_group_sizes": {name: len(selected) for name, selected in groups.items()},
        "two_loop_UV_reference_types": tuple(sorted(two_loop_types)),
        "complete_order_two_term_ownership": [
            {
                "term": "three bare full quartic vertices in two light loops",
                "owner": "the four disjoint raw-graph groups, 192 refinements",
            },
            {
                "term": "one-loop interaction counterterms inside a one-loop graph",
                "owner": "all 1152 marked one-loop vertex/heavy-mass cores; finite quartic contact separately",
            },
            {
                "term": "one-loop light quadratic counterterms inside a one-loop graph",
                "owner": "the inner on-shell grouping of the 64-refinement insertion family",
            },
            {
                "term": "two disjoint first-order tree counterterms",
                "owner": "72 actual disjoint-core contractions; cubic-cubic, cubic-mass and mass-mass products",
            },
            {
                "term": "new two-loop local interaction references",
                "owner": "sum of the four source-pinned linear family references, with G rather than g as the canonical cubic vertex",
            },
            {
                "term": "once-fixed finite one-loop quartic insertion",
                "owner": "fixed-g,M quartic variation of the full one-loop amplitude and inherited local reference",
            },
            {
                "term": "new fixed two-loop quartic potential contact",
                "owner": "one momentum-independent local contact; its forward second coefficient is zero",
            },
            {
                "term": "two-loop light mass/residue and external self-energy chains",
                "owner": "the complete inherited two-loop on-shell pole calculation and canonical LSZ cancellation",
            },
            {
                "term": "light-reducible tree exchanges",
                "owner": "no internal light tree exchange because the symmetric light three-point kernel vanishes; heavy-reducible exchanges remain in the Gaussian-heavy quartic graphs",
            },
        ],
        "scope": "The complete canonical order-two ownership ledger, not merely a bare-graph count. All regulator-dependent references are retained until cancellation. No extra finite heavy mass/cubic prescription is introduced, no momentum-dependent local contact is used to tune b2, and no source-aware derivative map or all-orders claim is inferred.",
        "checks": checks,
    }
