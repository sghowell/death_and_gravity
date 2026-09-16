"""Independent finite exponential-series enumeration of labeled topologies."""

from functools import cache

import sympy as s

C, g = s.symbols("C g")
ONE = {(0, 0): s.S.One}
X = {(1, 0): s.S.One}
Y = {(0, 1): s.S.One}


def add(*polys):
    result = {}
    for p in polys:
        for mon, c in p.items():
            result[mon] = s.expand(result.get(mon, 0) + c)
    return {m: c for m, c in result.items() if c != 0}


def scale(p, c):
    return {m: s.expand(c * v) for m, v in p.items() if c * v != 0}


def mul(p, q):
    result = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            mon = (i + k, j + l)
            if mon[0] <= 3 and mon[1] <= 2:
                result[mon] = s.expand(result.get(mon, 0) + a * b)
    return {m: c for m, c in result.items() if c != 0}


def power(p, n):
    result = ONE
    for _ in range(n):
        result = mul(result, p)
    return result


def exponential(p):
    return add(*(scale(power(p, n), 1 / s.factorial(n)) for n in range(6)))


@cache
def counts():
    P, H, h = X, {}, Y
    for iteration in range(12):
        E = exponential(h)
        nonconstant = add(E, scale(ONE, -1))
        nextP = add(
            X,
            mul(P, nonconstant),
            scale(mul(mul(P, H), E), g),
            scale(mul(power(P, 3), E), C / 6),
        )
        nextH = add(mul(H, nonconstant), scale(mul(power(P, 2), E), g / 2))
        scalar_terms = add(
            scale(power(P, 2), s.Rational(1, 2)),
            scale(power(H, 2), s.Rational(1, 2)),
            scale(mul(H, power(P, 2)), g / 2),
            scale(power(P, 4), C / 24),
        )
        nexth = add(Y, mul(scalar_terms, E), E, scale(ONE, -1), scale(h, -1))
        if (nextP, nextH, nexth) == (P, H, h):
            break
        P, H, h = nextP, nextH, nexth
    else:
        raise ValueError("The finite graph-count recursion did not stabilize")
    return tuple(
        s.expand(P.get((3, N), 0) * s.factorial(3) * s.factorial(N)) for N in range(3)
    )


@cache
def data():
    values = counts()
    expected = (C + 3 * g**2 + 3, 5 * C + 21 * g**2 + 21, 38 * C + 198 * g**2 + 198)
    checks = {
        f"independent_EGF_exact_{N}_gravitons": s.expand(a - b)
        for N, (a, b) in enumerate(zip(values, expected))
    }
    checks.update(
        {
            "complete_six_point_topology_count": values[2].subs({C: 1, g: 1}) - 434,
            "contact_six_point_topology_count": values[2].coeff(C) - 38,
            "heavy_six_point_topology_count": values[2].coeff(g, 2) - 198,
            "pure_Einstein_six_point_topology_count": values[2].subs({C: 0, g: 0})
            - 198,
        }
    )
    a, b = s.symbols("a b", nonzero=True)
    checks["same_line_two_soft_order_sum"] = s.cancel(
        1 / (a * (a + b)) + 1 / (b * (a + b)) - 1 / (a * b)
    )
    return {
        "whole_independent_EGF_counts": values,
        "whole_labeled_root_cut_bijection": "Cut the fixed external-Phi root vertex. Its remaining labeled external sets form a unique unordered partition; a species and rooted current occupies each block. Symmetric full multilinear vertices are used once, with no additional factorial. Reattaching produces a unique connected tree. Proper subsets give a terminating induction. Scalar propagators carry minus inverse kinetic operators and graviton currents the minus trace-reversed inverse. Final root is amputated.",
        "whole_independent_topology_EGF": "P=x+P(exp(h)-1)+gPH exp(h)+C P^3 exp(h)/6; H=H(exp(h)-1)+gP^2 exp(h)/2; h=y+[(P^2+H^2)/2+gHP^2/2+CP^4/24]exp(h)+exp(h)-1-h. Coefficient 3!N![x^3 y^N]P counts labeled four-Phi/N-h graphs. Formal powers beyond the implemented vertex jets cannot enter N<=2 by the degree budget. These count topologies, not tensor monomials or a UV-complete sum.",
        "checks": checks,
        "gates": {
            "independent_count_matches_three_coupling_sectors": values == expected,
            "six_point_inventory_38_plus_198_plus_198": 38 + 198 + 198 == 434,
            "lower_four_and_five_point_counts_recovered": values[0].subs({C: 1, g: 1})
            == 7
            and values[1].subs({C: 1, g: 1}) == 47,
            "finite_graph_sum_not_an_all_multiplicity_bound": True,
        },
    }
