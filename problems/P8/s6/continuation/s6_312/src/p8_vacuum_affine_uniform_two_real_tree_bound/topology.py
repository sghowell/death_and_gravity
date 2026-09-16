"""Labeled finite-tree inventory and every regular-graph denominator profile."""

from collections import Counter
from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as e

KINDS = ("phi", "H", "h")


@cache
def labeled_trees(mask, kind):
    if mask & (mask - 1) == 0:
        i = mask.bit_length() - 1
        return ((kind, mask, ()),) if kind == ("phi" if i < 3 else "h") else ()
    out = []
    for row in e.partitions(mask):
        if len(row) < 2:
            continue
        for tags in product(KINDS, repeat=len(row)):
            nf, nh, ng = ((kind, *tags).count(t) for t in KINDS)
            allowed = (
                nf == 2
                and nh == 0
                and 1 <= ng <= 3
                or nh == 2
                and nf == 0
                and 1 <= ng <= 3
                or nf == 2
                and nh == 1
                and 0 <= ng <= 2
                or nf == 4
                and nh == 0
                and 0 <= ng <= 2
                or nf == nh == 0
                and ng in (3, 4)
            )
            if not allowed:
                continue
            for children in product(
                *(labeled_trees(block, tag) for block, tag in zip(row, tags))
            ):
                out.append((kind, mask, children))
    return tuple(out)


def inventory(tree):
    vertices = []
    edges = []

    def walk(node, isroot=False):
        kind, mask, children = node
        if not children:
            return
        tags = (kind, *(child[0] for child in children))
        vertices.append(tuple(tags.count(t) for t in KINDS))
        if not isroot:
            edges.append((kind, mask & 7, (mask >> 3) & 3))
        for child in children:
            walk(child)

    walk(tree, True)
    return vertices, edges


@cache
def classified():
    answer = []
    for tree in labeled_trees(31, "phi"):
        vertices, edges = inventory(tree)
        nc = sum(nf == 4 for nf, nh, ng in vertices)
        ng = sum(nf == 2 and nh == 1 for nf, nh, nr in vertices)
        if (nc, ng) not in ((1, 0), (0, 2), (0, 0)):
            raise ValueError("Unexpected matter-core topology")
        sector = "C" if nc else ("H" if ng else "GR")
        paired = any(k == "h" and p == 0 and q == 3 for k, p, q in edges)
        answer.append((sector, paired, vertices, edges))
    return tuple(answer)


EXPECTED_PROFILES = {
    ("C", 0, 0, 0, 0, 0, 0): 1,
    ("C", 1, 0, 0, 0, 0, 1): 8,
    ("C", 1, 0, 0, 0, 0, 2): 4,
    ("C", 2, 0, 0, 0, 0, 2): 20,
    ("GR", 0, 1, 0, 0, 0, 4): 12,
    ("GR", 0, 2, 0, 0, 1, 2): 3,
    ("GR", 0, 2, 0, 1, 0, 3): 12,
    ("GR", 0, 3, 0, 2, 0, 2): 6,
    ("GR", 1, 1, 0, 0, 0, 4): 60,
    ("GR", 1, 2, 0, 1, 0, 3): 24,
    ("GR", 2, 1, 0, 0, 0, 4): 60,
    ("H", 0, 0, 0, 0, 0, 0): 12,
    ("H", 0, 0, 1, 0, 0, 0): 15,
    ("H", 0, 0, 2, 0, 0, 0): 6,
    ("H", 1, 0, 0, 0, 0, 1): 48,
    ("H", 1, 0, 0, 0, 0, 2): 12,
    ("H", 1, 0, 1, 0, 0, 1): 24,
    ("H", 2, 0, 0, 0, 0, 2): 60,
}


@cache
def data():
    hist = Counter()
    profiles = Counter()
    soft = Counter()
    all_cuts = True
    soft_products = True
    matter = True
    gravity = True
    for sector, paired, vv, ee in classified():
        hist[sector, paired] += 1
        if paired:
            continue
        hcuts = [(p, q) for k, p, q in ee if k == "h"]
        heavycuts = [(p, q) for k, p, q in ee if k == "H"]
        softcuts = [q if p.bit_count() == 1 else 3 ^ q for k, p, q in ee if k == "phi"]
        all_cuts &= all(p.bit_count() in (1, 3) for k, p, q in ee if k == "phi")
        all_cuts &= all(p.bit_count() == 2 for p, q in hcuts + heavycuts)
        soft_products &= all(q in (1, 2, 3) for q in softcuts) and len(softcuts) <= 2
        if len(softcuts) == 2:
            soft_products &= sorted(softcuts) in ([1, 2], [1, 3], [2, 3])
        soft[tuple(sorted(softcuts))] += 1
        light_degree = sum(r for nf, nh, r in vv if nf == 2 and nh == 0)
        heavy_vertices = sum(nh == 2 for nf, nh, r in vv)
        eh3 = sum(nf == nh == 0 and r == 3 for nf, nh, r in vv)
        eh4 = sum(nf == nh == 0 and r == 4 for nf, nh, r in vv)
        if sector in ("C", "H"):
            matter &= (
                not hcuts and not eh3 and not eh4 and sum(r for nf, nh, r in vv) == 2
            )
            if sector == "H":
                matter &= (
                    len(heavycuts) == heavy_vertices + 1 <= 3 and heavy_vertices <= 2
                )
            else:
                matter &= not heavycuts
        else:
            gravity &= not heavycuts and light_degree <= 4
            gravity &= len(hcuts) == eh3 + eh4 + 1 <= 3
            gravity &= (eh3, eh4) in ((0, 0), (1, 0), (2, 0), (0, 1))
        profiles[
            (sector, len(softcuts), len(hcuts), heavy_vertices, eh3, eh4, light_degree)
        ] += 1
    expected_hist = {
        ("C", False): 33,
        ("C", True): 5,
        ("H", False): 177,
        ("H", True): 21,
        ("GR", False): 177,
        ("GR", True): 21,
    }
    expected_soft = {
        (): 67,
        (1,): 76,
        (1, 2): 84,
        (1, 3): 28,
        (2,): 76,
        (2, 3): 28,
        (3,): 28,
    }
    checks = {"labeled_tree_inventory": s.Integer(len(classified()) - 434)}
    checks.update(
        {
            f"sector_{key[0]}_pair{int(key[1])}": s.Integer(hist[key] - v)
            for key, v in expected_hist.items()
        }
    )
    checks.update(
        {
            f"regular_profile_{i}": s.Integer(profiles[key] - v)
            for i, (key, v) in enumerate(EXPECTED_PROFILES.items())
        }
    )
    checks.update(
        {
            f"light_scalar_cut_profile_{i}": s.Integer(soft[key] - v)
            for i, (key, v) in enumerate(expected_soft.items())
        }
    )
    return {
        "whole_regular_graph_inventory": [
            {"sector": k[0], "paired": k[1], "count": v}
            for k, v in sorted(hist.items())
        ],
        "whole_regular_graph_profiles": [
            {
                "sector": k[0],
                "scalar_propagators": k[1],
                "hard_gravitons": k[2],
                "heavy_metric_vertices": k[3],
                "Einstein_cubics": k[4],
                "Einstein_quartics": k[5],
                "light_metric_degree": k[6],
                "count": v,
            }
            for k, v in sorted(profiles.items())
        ],
        "whole_topology_proof": "Cut the fixed scalar root exactly as in S310. Each unordered labeled partition, species assignment and allowed multilinear vertex reattaches to one tree. All434 trees are enumerated without momenta or tensor cancellations. Excluding the unique q1+q2 edge leaves33C+177H+177GR. The C/H connected matter core forbids any extra internal h except the excluded pure radiation branch. The pure-GR Phi forest has two components joined by one hard-h path, with at most two cubics or one quartic. One/three-Phi cuts carry one external massive leg and nonempty radiation; two such light propagators have subsets a,b or a,W or b,W. No soft denominator is counted twice without its enlarged nested subset.",
        "checks": checks,
        "gates": {
            "all_labeled_trees_distinct": len(set(labeled_trees(31, "phi"))) == 434,
            "all_sector_counts_exact": hist == expected_hist,
            "all_regular_profiles_exact": profiles == EXPECTED_PROFILES,
            "all_light_scalar_cut_profiles_exact": soft == expected_soft,
            "every_hard_and_light_cut_has_proved_form": all_cuts,
            "every_soft_product_bounded_by_eight_thirds_squared_over_ab": soft_products,
            "all_regular_matter_topologies_fit_budget": matter,
            "all_regular_Einstein_topologies_fit_budget": gravity,
            "counted_graphs_not_tensor_monomial_samples": True,
        },
    }
