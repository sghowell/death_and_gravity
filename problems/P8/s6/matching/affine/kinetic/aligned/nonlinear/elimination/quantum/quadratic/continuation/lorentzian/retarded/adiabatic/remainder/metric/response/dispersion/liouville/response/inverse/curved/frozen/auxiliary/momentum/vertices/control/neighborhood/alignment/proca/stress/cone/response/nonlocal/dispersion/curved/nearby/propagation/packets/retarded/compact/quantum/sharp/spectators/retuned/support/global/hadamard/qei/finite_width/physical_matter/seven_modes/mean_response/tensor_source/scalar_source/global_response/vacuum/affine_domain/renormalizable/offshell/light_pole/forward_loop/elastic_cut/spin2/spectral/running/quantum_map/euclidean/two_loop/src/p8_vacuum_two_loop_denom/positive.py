"""Deterministic exact positive-polynomial witnesses, not numerical sampling."""

import hashlib
from functools import cache

import sympy as sp

from . import forests, graphs


def data(kind, choices):
    graphs.refinement(kind, choices)
    return _data(kind, choices)


@cache
def _data(kind, choices):
    d = forests.data(kind, choices)
    a = d["parameters"]
    gap = sp.expand(d["coarse_mass_sixteen_F"] - d["U"] * sum(a) / 4)
    base = sp.Integer(0)
    base_squares = []
    checks = {}
    if kind == "tadpole_insertion":
        loop, edge, side1, side2 = a[:4]
        A = edge + side1 + side2
        base = sp.Rational(3, 4) * loop * ((edge - side1 - side2) ** 2 + loop * A)
        base_squares.append((sp.Rational(3, 4) * loop, edge - side1 - side2))
        checks["unrefined_tadpole_base_gap"] = sp.expand(
            gap.subs({v: 0 for v in a[4:]}) - base
        )
        if choices[2] in (2, 3):
            h = a[-1]
            base += (
                sp.Rational(3, 4) * h * (edge - side1 - side2 - loop) ** 2
                + sp.Rational(63, 4) * h * loop * A
                + sp.Rational(63, 4) * h * h * (A + loop)
            )
            base_squares.append((sp.Rational(3, 4) * h, edge - side1 - side2 - loop))
            checks["mixed_self_energy_base_gap"] = sp.expand(
                gap.subs({v: 0 for v in a[4:-1]}) - base
            )
    remainder = sp.Poly(sp.expand(gap - base), *a)
    coefficients = dict(remainder.terms())
    negative = [
        (powers, -coefficient / 2)
        for powers, coefficient in coefficients.items()
        if coefficient < 0
    ]
    if any(sorted(powers, reverse=True)[:3] != [1, 1, 1] for powers, _ in negative):
        raise ValueError("Unexpected negative monomial outside the exact witness basis")
    capacities = {powers: c for powers, c in coefficients.items() if c > 0}
    squares = []
    for powers, need in negative:
        indices = [i for i, n in enumerate(powers) if n]
        candidates = []
        for k in indices:
            i, j = [v for v in indices if v != k]
            p1, p2 = list(powers), list(powers)
            p1[j] -= 1
            p1[i] += 1
            p2[i] -= 1
            p2[j] += 1
            candidates.append((k, i, j, tuple(p1), tuple(p2)))
        candidates.sort(
            key=lambda row: min(capacities.get(row[3], 0), capacities.get(row[4], 0)),
            reverse=True,
        )
        for k, i, j, p1, p2 in candidates:
            amount = min(need, capacities.get(p1, 0), capacities.get(p2, 0))
            if amount > 0:
                squares.append((amount, k, i, j))
                capacities[p1] -= amount
                capacities[p2] -= amount
                need -= amount
        if need:
            raise ValueError("No exact deterministic positive-polynomial witness")
    square_polynomial = sum(w * a[k] * (a[i] - a[j]) ** 2 for w, k, i, j in squares)
    residual = sp.Poly(sp.expand(gap - base - square_polynomial), *a)
    if any(c < 0 for _, c in residual.terms()):
        raise ValueError("An exact positive remainder coefficient was negative")
    positive_base_remainder = sp.Poly(
        sp.expand(base - sum(w * linear**2 for w, linear in base_squares)), *a
    )
    if any(c < 0 for _, c in positive_base_remainder.terms()):
        raise ValueError("The explicit tadpole base was not positive")
    checks["complete_exact_positive_polynomial_identity"] = sp.expand(
        gap - base - square_polynomial - residual.as_expr()
    )

    def digest(expr):
        return hashlib.sha256(str(sp.Poly(expr, *a).terms()).encode()).hexdigest()

    return {
        "kind": kind,
        "choices": choices,
        "gap": gap,
        "explicit_tadpole_base": base,
        "base_square_terms": tuple(base_squares),
        "pair_square_terms": tuple(squares),
        "nonnegative_monomial_remainder": residual.as_expr(),
        "nonnegative_base_monomial_remainder": positive_base_remainder.as_expr(),
        "summary": {
            "kind": kind,
            "choices": choices,
            "heavy_edge_count": d["graph"]["heavy_edge_count"],
            "spanning_tree_count": d["spanning_tree_count"],
            "spanning_two_forest_count": d["spanning_two_forest_count"],
            "explicit_base_square_count": len(base_squares),
            "pair_square_terms": tuple(squares),
            "positive_remainder_monomial_count": sum(
                bool(c > 0) for _, c in residual.terms()
            ),
            "positive_base_monomial_count": sum(
                bool(c > 0) for _, c in positive_base_remainder.terms()
            ),
            "U_polynomial_sha256": digest(d["U"]),
            "actual_F_polynomial_sha256": digest(d["F"]),
            "gap_polynomial_sha256": digest(gap),
            "all_square_weights_and_remainder_coefficients_nonnegative": True,
        },
        "checks": checks,
    }


@cache
def summaries():
    return tuple(data(kind, choices)["summary"] for kind, choices in graphs.cases())
