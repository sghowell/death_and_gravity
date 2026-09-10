"""Disjoint primitive/forest ownership and regulator-sensitive finite terms."""

from functools import cache

import sympy as sp


@cache
def data():
    h, eps = sp.symbols("loop_order regulator_epsilon")
    a0, a1, a2, pole, finite = sp.symbols(
        "subgraph_0 subgraph_epsilon subgraph_epsilon_squared outer_pole outer_finite"
    )
    product = (a0 + eps * a1 + eps * eps * a2) * (pole / eps + finite)
    double_pole = sp.symbols("full_two_loop_double_pole")
    double_product = (a0 + eps * a1 + eps * eps * a2) * (
        double_pole / eps**2 + pole / eps + finite
    )
    z1, z2 = sp.symbols("first_field_increment second_field_increment")
    A0, A1, A2 = sp.symbols("tree_amplitude first_amplitude second_amplitude")
    normalized = (A0 + h * A1 + h * h * A2) / (1 + h * z1 + h * h * z2) ** 2
    second = sp.diff(normalized, h, 2).subs(h, 0) / 2
    G, G1, G2 = sp.symbols("fundamental_G first_G_reference second_G_reference")
    squared_cubic = (G + h * G1 + h * h * G2) ** 2
    Fraw, CT, Rest = sp.symbols(
        "raw_inner_piece assigned_inner_CT other_CT_piece", commutative=False
    )
    primitive_plus_ct = Fraw + CT + Rest
    grouped = (Fraw + CT) + Rest
    a, b, c, d = sp.symbols(
        "fermion_kernel11 fermion_kernel12 fermion_kernel21 fermion_kernel22"
    )
    K = sp.Matrix([[a, b], [c, d]])
    V = sp.Matrix([[1, 2], [3, 4]])
    det_derivative = -sp.diff(sp.log((K + h * V).det()), h).subs(h, 0)
    checks = {
        "finite_epsilon_squared_times_double_pole_retained": sp.expand(
            double_product
        ).coeff(eps, 0)
        - a0 * finite
        - a1 * pole
        - a2 * double_pole,
        "fundamental_cubic_second_order_square_retained": sp.diff(
            squared_cubic, h, 2
        ).subs(h, 0)
        / 2
        - G1 * G1
        - 2 * G * G2,
        "finite_epsilon_times_pole_retained": sp.expand(product).coeff(eps, 0)
        - a0 * finite
        - a1 * pole,
        "missing_epsilon_term_changes_finite_answer": sp.expand(
            product - (a0) * (pole / eps + finite)
        ).coeff(eps, 0)
        - a1 * pole,
        "second_order_canonical_external_field_formula": sp.factor(
            second - A2 + 2 * z1 * A1 - (3 * z1 * z1 - 2 * z2) * A0
        ),
        "regrouped_inner_counterterm_counted_once": sp.expand(
            primitive_plus_ct - grouped
        ),
        "double_counted_inner_counterterm_diagnostic": sp.expand(
            grouped + CT - primitive_plus_ct - CT
        ),
        "fermion_counterterm_determinant_variation_sign": sp.factor(
            det_derivative + sp.trace(K.inv() * V)
        ),
    }
    return {
        "primitive_mixed_functional": "1/2 Tr[D_B F_raw''], F_raw=-Tr log D_F, with the common regulator retained.",
        "bosonic_one_loop_counterterm_insertion": "1/2 Tr[D_B (S_B^(1))'']",
        "fermionic_one_loop_counterterm_insertion": "-Tr[D_F^-1 D_F^(1)]",
        "two_loop_tree_references": "S_B^(2), including allowed local mass, field, interaction, H one-point and vacuum-energy references in their declared scheme.",
        "second_order_canonical_external_field_coefficient": second,
        "second_order_squared_cubic_increment": G1 * G1 + 2 * G * G2,
        "finite_regulator_product": sp.expand(product).coeff(eps, 0),
        "finite_two_loop_double_pole_product": sp.expand(double_product).coeff(eps, 0),
        "partition_rule": "A proper-subgraph counterterm assigned inside a renormalized F0/F2/F4 or gauge block is removed from the separate insertion ledger for that same occurrence. Renormalized blocks plus a second copy of their counterterms are not the raw master expansion.",
        "four_dimensional_projection_boundary": "S6.135 removes its overall local quartic constant before the finite b2 projection; S6.136 removes the whole outer local affine polynomial. Other sectors, especially heavy-exchange prefactors and finite MS conversions, cannot assume all regulator-sensitive finite terms project to zero.",
        "scope": "Required ownership and order rules. The numerical or exact two-loop counterterm coefficients, complete finite MS conversion and full error budget are not supplied by this ledger.",
        "checks": checks,
    }
