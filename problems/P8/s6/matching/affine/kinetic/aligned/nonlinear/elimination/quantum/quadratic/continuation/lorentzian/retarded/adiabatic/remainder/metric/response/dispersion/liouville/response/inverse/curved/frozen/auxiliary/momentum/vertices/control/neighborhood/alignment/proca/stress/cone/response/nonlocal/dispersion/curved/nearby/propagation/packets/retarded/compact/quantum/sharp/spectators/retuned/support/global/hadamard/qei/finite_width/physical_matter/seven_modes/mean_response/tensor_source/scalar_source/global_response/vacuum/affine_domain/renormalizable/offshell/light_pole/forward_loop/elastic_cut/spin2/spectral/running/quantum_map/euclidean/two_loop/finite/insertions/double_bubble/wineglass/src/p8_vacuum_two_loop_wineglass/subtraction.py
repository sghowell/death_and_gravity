"""Fixed common-regulator forest subtraction and its full-model locality."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import subtraction as previous

from . import selection


@cache
def data():
    inner, I0, W = sp.symbols(
        "regulated_inner_bubble same_regulated_reference outer_weight"
    )
    J, Jzero, outer = sp.symbols(
        "regulated_unsubtracted_wineglass regulated_zero_wineglass regulated_outer_light_bubble"
    )
    J_R, J0 = sp.symbols(
        "inner_subtracted_wineglass fixed_inner_subtracted_zero_reference"
    )
    L, g = sp.symbols("positive_quartic positive_cubic_squared", positive=True)
    hP, hQ, t2, t3 = sp.symbols(
        "external_heavy_inverse internal_heavy_inverse triangle_channel2 triangle_channel3"
    )
    A = (-L, g * hP, t2, t3)
    B = {0: -L, 3: g * hQ}
    total = sum(A[c[0]] * B[c[1]] * B[c[2]] for c in selection.cases())
    CP, CQ, T = -L + g * hP, -L + g * hQ, t2 + t3
    finite_numerator = L**2 * T + (CP + T) * (CQ**2 - L**2)
    checks = {
        "sixteen_full_vertex_terms_sum_exactly": sp.expand(total - (CP + T) * CQ**2),
        "asymptotic_and_heavy_decaying_numerator_split": sp.expand(
            total - CP * L**2 - finite_numerator
        ),
        "inner_subtracted_zero_reference_keeps_both_factors": (Jzero - I0**2).subs(
            Jzero, J0 + I0**2
        )
        - J0,
    }
    rows = []
    for row in selection.data()["selected_rows"]:
        terms = []
        for f in row["restricted_forest_index_sets"]:
            if row["has_overall_core"]:
                table = {(): J, (0,): -I0 * outer, (1,): -Jzero, (0, 1): I0**2}
                term = table[f]
            else:
                term = W * inner if not f else -W * I0
            terms.append({"forest": f, "term": term})
        summed = sp.expand(sum(r["term"] for r in terms))
        if row["has_overall_core"]:
            reduced = sp.expand(summed.subs({J: J_R + I0 * outer, Jzero: J0 + I0**2}))
            expected = J_R - J0
        else:
            reduced, expected = summed, W * (inner - I0)
        checks[
            "actual_inner_and_overall_forest_" + "".join(map(str, row["choices"]))
        ] = sp.expand(reduced - expected)
        rows.append(
            {
                "choices": row["choices"],
                "has_overall_core": row["has_overall_core"],
                "forest_terms": tuple(terms),
                "renormalized_form": expected,
            }
        )
    hs = sp.symbols("h_s h_t h_u")
    dL, dg = -3 * L**3 * J0, -(L**2) * g * J0
    local = -dL + dg * sum(hs)
    checks["overall_subtraction_is_local_full_model_tree_variation"] = sp.expand(
        local + L**2 * J0 * sum(-L + g * h for h in hs)
    )
    M, z = sp.symbols("heavy_mass_squared channel_invariant")
    tree = -L + g / (M - z)
    checks["literal_single_channel_local_parameter_variation"] = sp.factor(
        dL * sp.diff(tree, L) / 3
        + dg * sp.diff(tree, g)
        - (L**3 - L**2 * g / (M - z)) * J0
    )
    return {
        "same_regulated_light_reference_with_loop_measure": previous.data()[
            "regulated_zero_momentum_light_bubble"
        ]
        / (16 * sp.pi**2),
        "all_actual_forest_terms": tuple(rows),
        "complete_selected_vertex_numerator": total,
        "asymptotic_external_factor": CP * L**2,
        "heavy_decaying_remainder_numerator": finite_numerator,
        "fixed_inner_subtracted_zero_reference": J0,
        "overall_local_counterterms": {
            "delta_polynomial_quartic": dL,
            "delta_cubic_squared": dg,
            "delta_heavy_mass_squared": sp.Integer(0),
        },
        "overall_local_counterterm_amplitude": local,
        "scope": "The proper light bubble uses the inherited entire I0 reference. The two logarithmic overall cores use an explicitly fixed recursive zero-external-momentum reference J0 after inner subtraction. J0 is not separately assigned a finite value; its cancellation is performed at the same regulator before the convergent difference is integrated. Other finite-potential insertions and complete two-loop normalization remain outside this family.",
        "checks": checks,
    }
