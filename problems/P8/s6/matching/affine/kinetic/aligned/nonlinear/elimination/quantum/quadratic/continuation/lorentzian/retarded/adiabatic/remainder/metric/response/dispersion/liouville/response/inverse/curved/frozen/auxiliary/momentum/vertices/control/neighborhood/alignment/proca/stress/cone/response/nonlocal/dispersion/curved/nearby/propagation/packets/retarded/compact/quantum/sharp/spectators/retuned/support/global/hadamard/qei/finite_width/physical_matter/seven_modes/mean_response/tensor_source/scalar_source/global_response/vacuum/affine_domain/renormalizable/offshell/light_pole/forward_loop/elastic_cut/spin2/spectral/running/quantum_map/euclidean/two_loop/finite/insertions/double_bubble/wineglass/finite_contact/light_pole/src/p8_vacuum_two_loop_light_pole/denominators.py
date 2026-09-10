"""Independent two-forest polynomials and exact feasible-flow denominator squares."""

from functools import cache
from itertools import combinations

import sympy as sp
from p8_vacuum_two_loop_denom.graphs import connected

from . import graphs


def components(vertices, edges):
    unseen = set(vertices)
    result = []
    while unseen:
        block = {min(unseen)}
        for _ in vertices:
            for u, v in edges:
                if u in block or v in block:
                    block.update((u, v))
        unseen.difference_update(block)
        result.append(block)
    return result


def point(kind, choices):
    graphs.refine(kind, choices)
    return _point(kind, choices)


@cache
def _point(kind, choices):
    g = graphs.refine(kind, choices)
    vs, es = g["vertices"], g["edges"]
    a = sp.symbols("a0:" + str(len(es)), positive=True)
    w = sp.symbols("w0:" + str(len(es)))
    lookup = {v: i for i, v in enumerate(vs)}
    lap = sp.zeros(len(vs))
    for weight, (u, v) in zip(w, es):
        i, j = lookup[u], lookup[v]
        if i != j:
            lap[i, i] += weight
            lap[j, j] += weight
            lap[i, j] -= weight
            lap[j, i] -= weight
    source = next(v for v in vs if 0 in g["external"][v])
    sink = next(v for v in vs if 1 in g["external"][v])
    b = sp.zeros(len(vs), 1)
    b[lookup[source]] += 1
    b[lookup[sink]] -= 1
    red = lap[:-1, :-1]
    replacement = {ww: 1 / aa for ww, aa in zip(w, a)}
    product_a = sp.prod(a)
    U = sp.expand(product_a * red.det().subs(replacement))
    potentials = red.adjugate() * b[:-1, :]
    potential_num = [sp.expand(product_a * t.subs(replacement)) for t in potentials] + [
        sp.S.Zero
    ]
    J = tuple(
        sp.cancel((potential_num[lookup[u]] - potential_num[lookup[v]]) / aa)
        for aa, (u, v) in zip(a, es)
    )
    P_matrix = sp.expand(sum(b[i] * potential_num[i] for i in range(len(vs))))
    cotrees = tuple(
        pair
        for pair in combinations(range(len(es)), 2)
        if connected(vs, tuple(e for i, e in enumerate(es) if i not in pair))
    )
    U_tree = sp.expand(sum(a[i] * a[j] for i, j in cotrees))
    two_forests = []
    for selected in combinations(range(len(es)), len(vs) - 2):
        blocks = components(vs, tuple(es[i] for i in selected))
        if len(blocks) == 2 and all(
            not (source in block and sink in block) for block in blocks
        ):
            two_forests.append(selected)
    P = sp.expand(
        sum(
            sp.prod(aa for i, aa in enumerate(a) if i not in selected)
            for selected in two_forests
        )
    )
    flow = [sp.S.Zero] * len(es)
    selected_paths = ()
    if source != sink:
        found = graphs.paths(g, source, sink)
        selected_paths = next(
            (p, q)
            for p, q in combinations(found, 2)
            if {i for i, _ in p}.isdisjoint(i for i, _ in q)
        )
        for path in selected_paths:
            for i, sign in path:
                flow[i] += sp.Rational(sign, 2)
    divJ = {v: sp.S.Zero for v in vs}
    divf = {v: sp.S.Zero for v in vs}
    for j, (u, v) in enumerate(es):
        divJ[u] += J[j]
        divJ[v] -= J[j]
        divf[u] += flow[j]
        divf[v] -= flow[j]
    unused = sum(aa for aa, f in zip(a, flow) if f == 0)
    squares = sum(aa * (U * f - j) ** 2 for aa, f, j in zip(a, flow, J))
    sos = squares + U**2 * unused / 4
    M, s = sp.symbols("positive_heavy_mass_squared Minkowski_invariant")
    H = sum(a[3:])
    F = U * (sum(a[:3]) + M * H) - s * P
    checks = {
        "tree_and_Laplacian_U": sp.expand(U_tree - U),
        "two_forest_and_adjugate_P": sp.expand(P - P_matrix),
        "current_energy": sp.expand(sum(aa * j**2 for aa, j in zip(a, J)) - U * P),
        "feasible_flow_pairing": sp.expand(
            sum(aa * f * j for aa, f, j in zip(a, flow, J)) - P
        ),
        "all_flow_magnitudes_zero_or_half": int(
            all(abs(f) in (0, sp.Rational(1, 2)) for f in flow)
        )
        - 1,
        "literal_resistance_sum_of_squares": sp.expand(U**2 * sum(a) / 4 - U * P - sos),
        "actual_mass_first_sheet_lower_gap": sp.expand(
            U * (F.subs(s, 3) - U * sum(a) / 4 - (M - 1) * U * H) - 3 * sos
        ),
        "positive_U_coefficients": int(all(c > 0 for c in sp.Poly(U, *a).coeffs())) - 1,
        "nonnegative_P_coefficients": int(all(c >= 0 for c in sp.Poly(P, *a).coeffs()))
        - 1,
    }
    checks.update(
        {
            "current_conservation_" + str(v): sp.expand(divJ[v] - b[lookup[v]] * U)
            for v in vs
        }
    )
    checks.update(
        {"feasible_flow_conservation_" + str(v): divf[v] - b[lookup[v]] for v in vs}
    )
    if kind == "nested_tadpole":
        checks["external_P_zero_exactly_for_outer_constant"] = (
            int((P == 0) == (choices[0] in (0, 1))) - 1
        )
    return {
        "kind": kind,
        "choices": choices,
        "parameters": a,
        "cotrees": cotrees,
        "separating_two_forests": tuple(two_forests),
        "U": U,
        "P": P,
        "F": F,
        "disjoint_external_paths": selected_paths,
        "half_unit_feasible_flow": tuple(flow),
        "unit_current_numerator": J,
        "resistance_square_polynomial": sos,
        "scope": "Exact U^2 sum(a)/4 - U P is a sum of nonnegative weighted squares. For |s-1|<=2 and M>=1 this gives Re(F)>=U sum(a)/4+(M-1)U sum(a_heavy). It is a denominator bound, not by itself a UV convergence claim.",
        "checks": checks,
    }


@cache
def data():
    rows = [point(kind, choices) for kind, choices in graphs.cases()]
    return {
        "all_32_denominator_certificates": [
            {k: v for k, v in row.items() if k != "checks"} for row in rows
        ],
        "checks": {
            row["kind"] + "_" + "".join(map(str, row["choices"])) + "_" + k: v
            for row in rows
            for k, v in row["checks"].items()
        },
    }
