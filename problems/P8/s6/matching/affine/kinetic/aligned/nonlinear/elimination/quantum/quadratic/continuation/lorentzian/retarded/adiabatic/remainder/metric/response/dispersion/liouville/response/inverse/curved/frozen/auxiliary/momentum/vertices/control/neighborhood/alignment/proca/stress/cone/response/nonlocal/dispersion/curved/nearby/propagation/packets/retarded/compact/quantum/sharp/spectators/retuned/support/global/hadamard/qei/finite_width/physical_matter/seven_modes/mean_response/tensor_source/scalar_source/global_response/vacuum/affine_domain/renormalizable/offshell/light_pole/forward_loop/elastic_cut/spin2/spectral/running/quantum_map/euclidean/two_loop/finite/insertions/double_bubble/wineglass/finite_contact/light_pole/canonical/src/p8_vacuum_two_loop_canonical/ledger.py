"""Full two-loop counterterm ledger in fundamental canonical L,G,M variables."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import subtraction as first
from p8_vacuum_two_loop_double_bubble import forests as double
from p8_vacuum_two_loop_finite_contact import subtraction as finite
from p8_vacuum_two_loop_insertions import subtraction as inserted
from p8_vacuum_two_loop_wineglass import subtraction as wine

from . import contractions


@cache
def data():
    L, G, M = sp.symbols("positive_L positive_G positive_M", positive=True)
    I0 = sp.Symbol("regulated_zero_bubble", real=True)
    Ibar = I0 / (16 * sp.pi**2)
    alpha = sp.Symbol("fixed_alpha", positive=True)
    sigma, sigma2, J0 = sp.symbols(
        "once_fixed_finite_contact fixed_second_order_potential_contact fixed_recursive_wineglass_reference",
        real=True,
    )
    hs = sp.symbols("h_s h_t h_u", real=True)
    common = {
        "positive_quartic": L,
        "quartic_L": L,
        "positive_cubic_squared": G**2,
        "cubic_squared": G**2,
        "regulated_zero_bubble": I0,
        "same_regulated_reference": Ibar,
        "fixed_positive_self_energy_asymptotic_multiplier": alpha,
        "fixed_inner_subtracted_zero_reference": J0,
        "once_fixed_finite_contact": sigma,
        **dict(zip(("h_s", "h_t", "h_u"), hs)),
    }

    def normalize(expr):
        value = sp.sympify(expr)
        unknown = {s.name for s in value.free_symbols if s.name not in common}
        if unknown:
            raise ValueError(
                "Unexpected inherited counterterm symbols: " + str(sorted(unknown))
            )
        return value.xreplace({s: common[s.name] for s in value.free_symbols})

    old, ins, db, wg, fin = (
        first.data(),
        inserted.data(),
        double.data(),
        wine.data(),
        finite.data(),
    )
    L1UV = normalize(old["delta_polynomial_quartic_UV"])
    g1 = normalize(old["delta_cubic_squared_UV"])
    G1 = g1 / (2 * G)
    M1 = normalize(old["delta_heavy_mass_squared_UV"])
    L1 = L1UV + sigma
    families = []
    checks = {}
    for name, d, counterterm_key, amplitude_key in (
        (
            "insertion",
            ins,
            "total_outer_local_counterterms_for_this_group",
            "outer_counterterm_amplitude",
        ),
        (
            "double_bubble",
            db,
            "shared_core_overall_local_counterterms",
            "shared_core_overall_counterterm_amplitude",
        ),
        (
            "wineglass",
            wg,
            "overall_local_counterterms",
            "overall_local_counterterm_amplitude",
        ),
        (
            "finite_contact",
            fin,
            "local_reference_continuation_counterterms",
            "total_canonical_counterterm_amplitude_for_this_group",
        ),
    ):
        ct = d[counterterm_key]
        dl, dg, dm = [
            normalize(ct[k])
            for k in (
                "delta_polynomial_quartic",
                "delta_cubic_squared",
                "delta_heavy_mass_squared",
            )
        ]
        amplitude = normalize(d[amplitude_key])
        checks[name + "_literal_inherited_local_reference"] = sp.expand(
            amplitude + dl - dg * sum(hs) + G**2 * dm * sum(h * h for h in hs)
        )
        families.append(
            {
                "family": name,
                "delta_L2_linear": dl,
                "delta_G2_linear": dg / (2 * G),
                "delta_M2_linear": dm,
                "local_reference_amplitude": amplitude,
            }
        )
    L2 = sum(r["delta_L2_linear"] for r in families) + sigma2
    G2 = sum(r["delta_G2_linear"] for r in families)
    M2 = sum(r["delta_M2_linear"] for r in families)
    linear = -L2 + 2 * G * G2 * sum(hs) - G**2 * M2 * sum(h * h for h in hs)
    products = (
        G1**2 * sum(hs)
        - 2 * G * G1 * M1 * sum(h * h for h in hs)
        + G**2 * M1**2 * sum(h**3 for h in hs)
    )
    literal_products = sp.S.Zero
    channel_sums = {}
    for row in contractions.data()["all_disjoint_tree_cograph_weight_matches"]:
        pL, pG = row["coupling_powers_L_G"]
        endpoints = [
            set(labels)
            for _, labels in row["vertex_types_and_external_labels"]
            if labels
        ]
        assert len(endpoints) == 2 and all(len(v) == 2 for v in endpoints)
        label = (
            0
            if endpoints[0] in ({0, 1}, {2, 3})
            else 1
            if endpoints[0] in ({0, 2}, {1, 3})
            else 2
        )
        edges = row["typed_cograph_edges"]
        assert all(t == "H" for _, _, t in edges)
        term = (
            row["contracted_signed_coefficient"]
            * L**pL
            * G**pG
            * Ibar**2
            * hs[label] ** len(edges)
        )
        literal_products += term
        key = (label, pL, pG, len(edges))
        channel_sums[key] = (
            channel_sums.get(key, 0) + row["contracted_signed_coefficient"]
        )
    checks["actual_all_disjoint_graphs_equal_full_tree_CT_products"] = sp.expand(
        literal_products - products
    )
    for j in range(3):
        checks["channel_" + str(j) + "_two_cubic_counterterms"] = channel_sums[
            (j, 2, 2, 1)
        ] - sp.Rational(1, 4)
        checks["channel_" + str(j) + "_cubic_mass_cross_counterterms"] = channel_sums[
            (j, 1, 4, 2)
        ] + sp.Rational(1, 2)
        checks["channel_" + str(j) + "_two_mass_counterterms"] = channel_sums[
            (j, 0, 6, 3)
        ] - sp.Rational(1, 4)
    tau = sp.Symbol("formal_loop_order")
    l1, l2, c1, c2, m1, m2 = sp.symbols(
        "formal_L1 formal_L2 formal_G1 formal_G2 formal_M1 formal_M2"
    )
    tree = -(L + tau * l1 + tau**2 * l2) + (G + tau * c1 + tau**2 * c2) ** 2 * sum(
        h / (1 + (tau * m1 + tau**2 * m2) * h) for h in hs
    )
    tree_first = sp.diff(tree, tau).subs(tau, 0)
    tree_second = sp.diff(tree, tau, 2).subs(tau, 0) / 2
    target = (
        -l2
        + (2 * G * c2 + c1**2) * sum(hs)
        - (G**2 * m2 + 2 * G * c1 * m1) * sum(h * h for h in hs)
        + G**2 * m1**2 * sum(h**3 for h in hs)
    )
    checks.update(
        {
            "first_order_tree_parameter_variation": sp.expand(
                tree_first
                - (-l1 + 2 * G * c1 * sum(hs) - G**2 * m1 * sum(h * h for h in hs))
            ),
            "full_second_order_fundamental_G_tree_expansion": sp.expand(
                tree_second - target
            ),
            "actual_full_tree_CT_ledger": sp.expand(
                target.subs({l1: L1, l2: L2, c1: G1, c2: G2, m1: M1, m2: M2})
                - linear
                - products
            ),
            "complete_one_loop_UV_vertex_reference": sp.expand(
                -L1UV
                + 2 * G * G1 * sum(hs)
                - G**2 * M1 * sum(h * h for h in hs)
                + Ibar * sum((-L + G**2 * h) ** 2 for h in hs) / 2
            ),
            "first_cubic_vertex_CT_not_squared_parameter_CT": sp.factor(
                2 * G * G1 - g1
            ),
            "second_squared_coupling_includes_disjoint_cubic_product": sp.expand(
                sp.diff((G + tau * G1 + tau**2 * G2) ** 2, tau, 2).subs(tau, 0) / 2
                - (2 * G * G2 + G1**2)
            ),
            "insertion_reference_is_twice_alpha_old_UV_vertex": sp.expand(
                families[0]["local_reference_amplitude"]
                - 2 * alpha * normalize(old["full_UV_counterterm_amplitude"])
            ),
            "finite_contact_reference_is_quartic_variation": sp.expand(
                families[3]["local_reference_amplitude"]
                - sigma * sp.diff(normalize(old["full_UV_counterterm_amplitude"]), L)
            ),
            "linear_reference_sum_keeps_fixed_new_potential_contact": sp.expand(
                linear - sum(r["local_reference_amplitude"] for r in families) + sigma2
            ),
        }
    )
    s = sp.Symbol("forward_invariant", real=True)
    missing = G1**2 * (1 / (M - s) + 1 / M + 1 / (M - 4 + s))
    missing_b2 = sp.factor(sp.diff(missing, s, 2).subs(s, 2) / 2)
    checks["nonlocal_cubic_product_second_forward_coefficient"] = sp.factor(
        missing_b2 - 2 * G1**2 / (M - 2) ** 3
    )
    checks["new_fixed_second_order_contact_has_zero_b2"] = sp.diff(-sigma2, s, 2)
    return {
        "L": L,
        "G": G,
        "M": M,
        "regulated_I0_with_loop_measure": Ibar,
        "first_order_canonical_vertex_counterterms": {"L": L1, "G": G1, "M": M1},
        "second_order_inherited_linear_reference_families": families,
        "fixed_second_order_potential_contact": sigma2,
        "second_order_canonical_vertex_counterterms": {"L": L2, "G": G2, "M": M2},
        "linear_second_order_tree_counterterms": linear,
        "actual_disjoint_first_order_tree_products": literal_products,
        "full_second_order_tree_counterterm_amplitude": linear + products,
        "derived_total_squared_coupling_second_order": 2 * G * G2 + G1**2,
        "diagnostic_omitted_cubic_product_b2": missing_b2,
        "scope": "G1 and G2 are total CANONICAL cubic-vertex counterterms, before squaring the canonical vertex into g. The inherited second-order family delta-g labels describe linear tree variations 2G G2. All disjoint first-order products are retained separately and exactly once. In a derived g parametrization they must instead be included in total delta-g2. I0 and the recursive wineglass reference are common regulated symbols, never separately assigned finite values.",
        "checks": checks,
    }
