"""Actual restricted forest operations and local full-model counterterms."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import subtraction as previous

from . import selection


@cache
def data():
    I, I0, T = sp.symbols(
        "regulated_light_bubble same_regulated_reference finite_triangle"
    )
    L, g = sp.symbols("positive_quartic positive_cubic_squared", positive=True)
    h = sp.Symbol("external_heavy_inverse")
    a2, a3, b2, b3 = sp.symbols(
        "left_triangle2 left_triangle3 right_triangle2 right_triangle3"
    )
    C = -L + g * h
    left = (-L * I, g * h * I, a2, a3)
    right = (-L * I, g * h * I, b2, b3)
    center = (-L, g * h)
    raw = sum(center[c[2]] * left[c[0]] * right[c[1]] for c in selection.cases()) / 4
    target = (C**3 * I**2 + C**2 * I * (a2 + a3 + b2 + b3)) / 4
    rows, checks = (
        [],
        {
            "raw_selected_refinements_equal_factorization": sp.expand(raw - target),
            "eight_finite_refinements_removed_before_bounding": sp.expand(
                C * ((C * I + a2 + a3) * (C * I + b2 + b3) - (a2 + a3) * (b2 + b3)) / 4
                - target
            ),
        },
    )
    for row in selection.data()["selected_rows"]:
        mode = row["mode"]
        uv = row["UV_subgraphs"]
        terms = []
        for forest in row["restricted_forest_index_sets"]:
            if mode == "single":
                term = I * T if not forest else -I0 * T
            elif mode == "disjoint":
                term = (-I0) ** len(forest) * I ** (2 - len(forest))
            else:
                whole = next(i for i, s in enumerate(uv) if len(s["edges"]) == 4)
                if whole in forest:
                    term = (-1) ** len(forest) * I0**2
                else:
                    term = I**2 if not forest else -I0 * I
            terms.append({"forest": forest, "term": term})
        total = sp.expand(sum(t["term"] for t in terms))
        expected = (I - I0) * T if mode == "single" else (I - I0) ** 2
        key = "".join(map(str, row["choices"]))
        checks["every_actual_restricted_forest_sum_" + key] = sp.expand(
            total - expected
        )
        rows.append(
            {
                "choices": row["choices"],
                "mode": mode,
                "terms": tuple(terms),
                "sum": total,
            }
        )
    Rprime = I**2 - 2 * I0 * I
    checks["overlap_recursion_overall_counterterm_positive_reference_square"] = (
        -Rprime.subs(I, I0) - I0**2
    )
    checks["overlap_recursion_equals_product_only_after_overall_subtraction"] = (
        sp.expand(Rprime - Rprime.subs(I, I0) - (I - I0) ** 2)
    )
    hs = sp.symbols("h_s h_t h_u")
    dL, dg, dM = 3 * L**3 * I0**2 / 4, L**2 * g * I0**2 / 2, L * g * I0**2 / 4
    local = -dL + dg * sum(hs) - g * dM * sum(x**2 for x in hs)
    desired = sum(-L * (-L + g * x) ** 2 * I0**2 / 4 for x in hs)
    checks["shared_core_overall_subtractions_are_local_full_model_parameters"] = (
        sp.expand(local - desired)
    )
    sigma = sp.Symbol("loop_order")
    M, s = sp.symbols("heavy_mass_squared channel_invariant")
    tree = -L + g / (M - s)
    literal = dL * sp.diff(tree, L) / 3 + dg * sp.diff(tree, g) + dM * sp.diff(tree, M)
    checks["single_channel_local_counterterm_literal_tree_variation"] = sp.factor(
        literal - (-dL / 3 + dg * h - g * dM * h**2).subs(h, 1 / (M - s))
    )
    mass_ct = sp.Symbol("one_loop_heavy_mass_counterterm")
    expanded = sp.diff(g / (M + sigma * mass_ct - s), sigma, 2).subs(sigma, 0) / 2
    checks["disjoint_mass_counterterm_product_is_allowed_iterated_propagator"] = (
        sp.factor(expanded - g * mass_ct**2 / (M - s) ** 3)
    )
    return {
        "I": I,
        "I0": I0,
        "C": C,
        "same_full_regulated_reference_with_loop_measure": previous.data()[
            "regulated_zero_momentum_light_bubble"
        ]
        / (16 * sp.pi**2),
        "all_actual_forest_terms": tuple(rows),
        "selected_raw_amplitude_per_channel": target,
        "renormalized_selected_amplitude_per_channel": target.subs(I, I - I0),
        "shared_core_overall_local_counterterms": {
            "delta_polynomial_quartic": dL,
            "delta_cubic_squared": dg,
            "delta_heavy_mass_squared": dM,
        },
        "shared_core_overall_counterterm_amplitude": local,
        "scope": "At one common regulator each logarithmic light-bubble core is subtracted at zero external momentum with the inherited entire I0 reference. The shared-vertex two-loop cores use the explicitly extended zero-momentum recursive prescription. This does not include other finite-potential counterterm insertions or complete two-loop normalization.",
        "checks": checks,
    }
