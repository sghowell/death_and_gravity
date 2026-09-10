"""Exact Hepp-sector anisotropic-box upper polynomials, with full enumeration."""

import hashlib
from fractions import Fraction
from functools import cache
from itertools import permutations
from math import factorial

from . import selection


def single(kind, choices, order):
    graph = selection.require_finite(kind, choices)
    size = len(graph["edges"])
    if (
        type(order) is not tuple
        or len(order) != size
        or any(type(j) is not int for j in order)
        or set(order) != set(range(size))
    ):
        raise ValueError("Require a native permutation of all refinement edges")
    return _single(kind, choices, order)


def _single(kind, choices, order):
    size = len(order)
    rank = {edge: i for i, edge in enumerate(order)}
    ranks = tuple(
        tuple(sorted(rank[e] for e in pair))
        for pair in selection.cotrees(kind, choices)
    )
    basis = min(ranks)
    # Componentwise dominance is independently verified, not assumed from lex order.
    if not all(all(a <= b for a, b in zip(basis, other)) for other in ranks):
        raise ValueError("The proposed cotree does not dominate this ordered sector")
    powers = tuple(-2 if i in basis else 0 for i in range(size))
    tails = tuple(sum(powers[i] + 1 for i in range(j, size)) for j in range(size))
    if min(tails) <= 0:
        raise ValueError("Nonpositive ultraviolet sector exponent")
    initial = next(i for i, e in enumerate(order) if e >= 4)
    coefficients = [Fraction(0)] * 3
    pieces = []
    for k in range(initial + 1):
        prefix = tuple(sum(powers[i] + 1 for i in range(j + 1)) for j in range(k))
        if any(a > 0 for a in prefix) or prefix.count(0) > 2:
            raise ValueError("Unsupported power growth above the heavy cutoff")
        zero = prefix.count(0)
        denominator = factorial(zero)
        for a in tails[k:]:
            denominator *= a
        for a in prefix:
            if a < 0:
                denominator *= -a
        coefficient = Fraction(1, denominator)
        coefficients[zero] += coefficient
        if size - 4 - tails[k] != sum(powers[i] + 1 for i in range(k)):
            raise ValueError("The heavy-mass scaling failed to cancel")
        pieces.append(
            {
                "large_light_count": k,
                "prefix_exponents": prefix,
                "small_tail_exponents": tails[k:],
                "log_power": zero,
                "rational_coefficient": coefficient,
            }
        )
    return {
        "order": order,
        "dominating_cotree_ranks": basis,
        "monomial_powers": powers,
        "strictly_positive_small_tail_exponents": tails,
        "pieces": tuple(pieces),
        "box_log_polynomial_coefficients": tuple(coefficients),
    }


def data(kind, choices):
    selection.require_finite(kind, choices)
    return _data(kind, choices)


@cache
def _data(kind, choices):
    graph = selection.require_finite(kind, choices)
    size = len(graph["edges"])
    coefficients = [Fraction(0)] * 3
    digest = hashlib.sha256()
    sector_count = piece_count = 0
    minimum_tail = size
    max_log = 0
    patterns = set()
    for order in permutations(range(size)):
        d = _single(kind, choices, order)
        for j, c in enumerate(d["box_log_polynomial_coefficients"]):
            coefficients[j] += c
        for row in d["pieces"]:
            patterns.add(d["monomial_powers"][: row["large_light_count"]])
            max_log = max(max_log, row["log_power"])
        minimum_tail = min(
            minimum_tail, min(d["strictly_positive_small_tail_exponents"])
        )
        digest.update((repr(d) + "\n").encode())
        sector_count += 1
        piece_count += len(d["pieces"])
    return {
        "kind": kind,
        "choices": choices,
        "heavy_edge_count": size - 4,
        "edge_order_sector_count": sector_count,
        "heavy_cutoff_piece_count": piece_count,
        "box_log_polynomial_coefficients": tuple(coefficients),
        "minimum_small_tail_exponent": minimum_tail,
        "maximum_log_power": max_log,
        "large_ordered_power_patterns": tuple(sorted(patterns)),
        "complete_sector_inventory_sha256": digest.hexdigest(),
        "every_monomial_dominance_tail_prefix_and_mass_cancellation_checked": True,
        "scope": "For B>=1, B^h times the positive anisotropic-box U^-2 integral is bounded by this rational polynomial in log B. This is a bound, not an exact value of the graph integral.",
    }


@cache
def summaries():
    return tuple(data(*case) for case in selection.cases())
