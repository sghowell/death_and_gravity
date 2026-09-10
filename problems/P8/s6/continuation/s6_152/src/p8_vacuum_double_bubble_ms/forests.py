"""Actual shared, disjoint and single-core MS forests, not a guessed product."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_double_bubble import selection


def poles(expr, e):
    if not isinstance(e, s.Symbol):
        raise TypeError("Require a symbolic regulator")
    expr = s.sympify(expr)
    if expr.has(s.Float):
        raise TypeError("Require exact Laurent coefficients")
    kept = []
    for term in s.Add.make_args(s.expand(expr)):
        power = term.as_powers_dict().get(e, s.Integer(0))
        if power.is_Integer is not True or s.cancel(term / e**power).has(e):
            raise ValueError("Require a finite Laurent polynomial in the regulator")
        if power < 0:
            kept.append(term)
    return s.Add(*kept)


@cache
def data():
    e, Q, f0, f1, f2, T = s.symbols(
        "epsilon Q finite_bubble evanescent_first evanescent_second full_finite_triangle"
    )
    P = 1 / (Q * e)
    I = P + f0 + e * f1 + e**2 * f2
    raw_poles = poles(I**2, e)
    nested_poles = poles(P * I, e)
    rows, checks = (
        [],
        {
            "whole_raw_poles_retain_finite_bubble": s.expand(
                raw_poles - P**2 - 2 * P * f0
            ),
            "nested_poles_retain_finite_bubble": s.expand(nested_poles - P**2 - P * f0),
            "proper_subtracted_overlap_poles": s.expand(
                poles(I**2 - 2 * P * I, e) + P**2
            ),
            "overlap_overall_CT_positive_pole_square": s.expand(
                -poles(I**2 - 2 * P * I, e) - P**2
            ),
            "six_actual_overlap_terms_not_forbidden_pair": s.expand(
                I**2 - 2 * P * I - raw_poles + 2 * nested_poles - (I - P) ** 2
            ),
            "finite_complete_square_keeps_no_unpaired_pole": s.expand(
                (I - P) ** 2
            ).coeff(e, 0)
            - f0**2,
        },
    )
    old = selection.data()
    for row in old["selected_rows"]:
        mode = row["mode"]
        terms = []
        whole = next(
            (i for i, uv in enumerate(row["UV_subgraphs"]) if len(uv["edges"]) == 4),
            None,
        )
        for f in row["restricted_forest_index_sets"]:
            if mode == "single":
                term = I * T if not f else -P * T
            elif mode == "disjoint":
                term = (-P) ** len(f) * I ** (2 - len(f))
            elif whole in f:
                term = -raw_poles if len(f) == 1 else nested_poles
            else:
                term = I**2 if not f else -P * I
            terms.append({"forest": f, "term": term})
        expected = (I - P) * T if mode == "single" else (I - P) ** 2
        checks["actual_MS_forest_" + "".join(map(str, row["choices"]))] = s.expand(
            sum(v["term"] for v in terms) - expected
        )
        if mode == "overlap":
            small = tuple(i for i in range(len(row["UV_subgraphs"])) if i != whole)
            checks[
                "forbidden_shared_pair_absent_" + "".join(map(str, row["choices"]))
            ] = s.Integer(small in row["restricted_forest_index_sets"])
        rows.append(
            {
                "choices": row["choices"],
                "mode": mode,
                "terms": terms,
                "renormalized_form": expected,
            }
        )
    checks.update(
        {
            "all_twenty_four_refinements": len(rows) - 24,
            "eight_finite_separable_refinements_excluded": old[
                "already_finite_separable_count"
            ]
            - 8,
        }
    )
    return {
        "symbols": {"epsilon": e, "Q": Q},
        "entire_regulated_bubble": I,
        "MS_inner_pole": P,
        "raw_overall_pole_part": raw_poles,
        "nested_overall_pole_part": nested_poles,
        "actual_MS_forest_rows": rows,
        "same_forest_class_counts": old["forest_class_counts"],
        "checks": checks,
        "scope": "All actual parent restricted forests are evaluated with MS pole projections. The shared pair is absent; raw and nested overall projections carry different finite-bubble simple poles before their cancellation. The finite triangle remains whole in a single-core subtraction.",
    }
