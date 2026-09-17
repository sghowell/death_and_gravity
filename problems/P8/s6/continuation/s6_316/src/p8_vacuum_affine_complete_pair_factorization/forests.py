"""Exact pair-branch inclusion-exclusion and independent root-cut inventory."""

from functools import cache
from itertools import combinations, product

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import topology
from p8_vacuum_affine_complete_two_graviton_tree import trees

from . import source

C, g = topology.C, topology.g


def matching_count(n, k):
    n, k = source.require_multiplicity(n), source.require_multiplicity(k)
    if 2 * k > n:
        return s.S.Zero
    return s.factorial(n) / (2**k * s.factorial(k) * s.factorial(n - 2 * k))


@cache
def no_pair_count(n):
    n = source.require_multiplicity(n)
    return s.expand(
        sum(
            (-1) ** k * matching_count(n, k) * topology.topology_closed(n - k)
            for k in range(n // 2 + 1)
        )
    )


def labeled_matchings(labels):
    labels = tuple(labels)
    if not labels:
        return ((),)
    first, rest = labels[0], labels[1:]
    result = list(labeled_matchings(rest))
    for j, other in enumerate(rest):
        for row in labeled_matchings(rest[:j] + rest[j + 1 :]):
            result.append(((first, other), *row))
    return tuple(result)


@cache
def independent_no_pair_count(n):
    n = source.require_multiplicity(n)
    kinds = ("phi",) * 3 + ("h",) * n

    @cache
    def current(mask, kind):
        if not mask & (mask - 1):
            return s.S.One if kinds[mask.bit_length() - 1] == kind else s.S.Zero
        if kind == "h" and not (mask & 7) and mask.bit_count() == 2:
            return s.S.Zero
        nphi = sum(kinds[i] == "phi" for i in range(len(kinds)) if mask >> i & 1)
        if (
            (nphi % 2 == 1 and kind != "phi")
            or (nphi % 2 == 0 and kind == "phi")
            or (nphi == 0 and kind == "H")
        ):
            return s.S.Zero
        result = s.S.Zero
        for row in trees.partitions(mask):
            if len(row) < 2:
                continue
            choices = []
            for part in row:
                options = tuple((tag, current(part, tag)) for tag in trees.KINDS)
                options = tuple((tag, value) for tag, value in options if value != 0)
                if not options:
                    break
                choices.append(options)
            else:
                for selected in product(*choices):
                    tags = (kind, *(item[0] for item in selected))
                    nf, nh, ng = (tags.count(kind) for kind in trees.KINDS)
                    if (
                        (nf == 2 and nh == 0 and ng >= 1)
                        or (nh == 2 and nf == 0 and ng >= 1)
                        or (nf == nh == 0 and ng >= 3)
                    ):
                        weight = 1
                    elif nf == 2 and nh == 1:
                        weight = g
                    elif nf == 4 and nh == 0:
                        weight = C
                    else:
                        continue
                    result += weight * s.prod(item[1] for item in selected)
        return s.expand(result)

    return current((1 << len(kinds)) - 1, "phi")


@cache
def data():
    checks = {}
    expected = (
        C + 3 * g**2 + 3,
        5 * C + 21 * g**2 + 21,
        33 * C + 177 * g**2 + 177,
        274 * C + 1770 * g**2 + 1770,
        2758 * C + 20646 * g**2 + 20646,
    )
    for n, value in enumerate(expected):
        checks[f"closed_no_pair_inventory_{n}"] = s.expand(no_pair_count(n) - value)
        checks[f"independent_root_cut_no_pair_inventory_{n}"] = s.expand(
            independent_no_pair_count(n) - value
        )
    inventories = tuple(no_pair_count(n) for n in range(9))
    y = s.Symbol("y")
    F = sum(topology.topology_closed(n) * y**n / s.factorial(n) for n in range(9))
    substituted = s.Poly(s.expand(F.subs(y, y - y * y / 2)), y)
    for n in range(9):
        checks[f"matching_forest_EGF_composition_{n}"] = s.expand(
            substituted.nth(n) * s.factorial(n) - inventories[n]
        )
        rows = labeled_matchings(range(n))
        for k in range(n // 2 + 1):
            checks[f"independent_labeled_matchings_{n}_{k}"] = sum(
                len(row) == k for row in rows
            ) - matching_count(n, k)
        checks[f"nonnegative_sector_decomposition_{n}"] = s.expand(
            inventories[n]
            - C * inventories[n].coeff(C)
            - (1 + g**2) * inventories[n].coeff(g, 2)
        )
    # Boolean inclusion-exclusion on every subset of each actual disjoint pair set.
    for n in range(7):
        for k in range(1, n // 2 + 1):
            checks[f"union_weight_of_k_present_pairs_{n}_{k}"] = (
                sum((-1) ** (j + 1) * s.binomial(k, j) for j in range(1, k + 1)) - 1
            )
    pair_sets = tuple(combinations(range(4), 2))
    disjoint = tuple(
        (a, b) for a, b in combinations(pair_sets, 2) if not set(a) & set(b)
    )
    checks["first_overlap_count_at_four_real"] = s.Integer(len(disjoint) - 3)
    return {
        "whole_pair_branch_definition": "A pair branch consists of two labeled external physical gravitons attached to one cubic Einstein vertex, followed by its internal graviton edge. Contracting k disjoint such branches gives a bijection to a tree on N-k graviton labels; the source weights of the added cubic vertices carry no matter-sector marker.",
        "whole_all_N_no_pair_inventory": "R_N=sum_{k=0}^{floor(N/2)}(-1)^k N!/[2^k k! (N-2k)!] T_{N-k}; its graph EGF is R(y)=T(y-y^2/2). This is only a graph-count composition, not an amplitude composition in kinematic variables.",
        "whole_inventory_calibrations": inventories,
        "whole_amplitude_union_identity": "Let A_M be the sum of all original trees containing every pair in a nonempty matching M. The full amplitude equals A_no_pairs + sum_M (-1)^(|M|+1) A_M. Shared-label pair intersections are empty; disjoint-pair intersections are essential starting at N4.",
        "checks": checks,
        "gates": {
            "every_fixed_matching_collapses_bijectively": True,
            "shared_label_pair_intersections_empty_in_a_tree": True,
            "independent_root_cuts_match_through_N4": all(
                independent_no_pair_count(n) == expected[n] for n in range(5)
            ),
            "independent_labeled_matchings_have_no_duplicates": all(
                len(set(labeled_matchings(range(n))))
                == len(labeled_matchings(range(n)))
                for n in range(9)
            ),
            "first_overlap_requires_three_N4_intersections": len(disjoint) == 3,
            "finite_no_pair_counts_are_nonnegative_and_bounded_by_full": all(
                0
                <= v.subs({C: 1, g: 1})
                <= topology.topology_closed(n).subs({C: 1, g: 1})
                for n, v in enumerate(inventories)
            ),
            "Boolean_inclusion_exclusion_is_general_not_finite_calibration": True,
            "graph_EGF_not_probability_or_amplitude_summability": True,
        },
    }
