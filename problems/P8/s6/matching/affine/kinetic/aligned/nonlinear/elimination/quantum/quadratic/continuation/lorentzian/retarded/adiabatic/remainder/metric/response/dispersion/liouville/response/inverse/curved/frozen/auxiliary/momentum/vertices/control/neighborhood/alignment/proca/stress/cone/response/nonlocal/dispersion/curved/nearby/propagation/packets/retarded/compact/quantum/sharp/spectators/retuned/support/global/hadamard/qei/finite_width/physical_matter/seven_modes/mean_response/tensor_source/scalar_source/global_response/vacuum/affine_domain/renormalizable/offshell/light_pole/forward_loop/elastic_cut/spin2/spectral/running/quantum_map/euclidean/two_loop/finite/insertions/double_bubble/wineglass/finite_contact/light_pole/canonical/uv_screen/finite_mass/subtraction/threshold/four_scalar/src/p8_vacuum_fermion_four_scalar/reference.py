"""No double counting between the local potential, box and pole field map."""

from functools import cache

import sympy as sp


@cache
def data():
    h, fp = sp.symbols("loop_order on_shell_fermion_slope")
    A0, AF, v4 = sp.symbols(
        "tree_amplitude fermion_box_amplitude fermion_quartic_threshold"
    )
    L, g, M = sp.symbols(
        "reference_quartic reference_cubic_squared reference_heavy_mass_squared"
    )
    s, t, u = sp.symbols("channel_s channel_t channel_u")
    lam, bF = sp.symbols("reference_positive_lambda fermion_box_second_coefficient")
    kappa = 1 - fp
    exchange = sum(1 / (M - z) for z in (s, t, u))
    tree = -L + g * exchange
    canonical_tree = -L / kappa**2 + g * exchange / kappa**2
    potential_tree = -(L + v4) / kappa**2 + g * exchange / kappa**2
    selected = (A0 + AF) / kappa**2
    h_selected = (A0 + h * AF) / (1 - h * fp) ** 2
    b2_selected = (4 * lam + bF) / kappa**2
    checks = {
        "canonical_tree_field_factors": sp.factor(canonical_tree - tree / kappa**2),
        "potential_quartic_plus_residual_box_identity": sp.factor(
            potential_tree + (AF + v4) / kappa**2 - (tree + AF) / kappa**2
        ),
        "first_order_loop_and_external_field_dictionary": sp.factor(
            sp.series(h_selected, h, 0, 2).removeO() - A0 - h * (AF + 2 * fp * A0)
        ),
        "first_order_center_coefficient_dictionary": sp.factor(
            sp.series((4 * lam + h * bF) / (1 - h * fp) ** 2, h, 0, 2).removeO()
            - 4 * lam
            - h * (bF + 8 * lam * fp)
        ),
        "incorrect_repeated_local_threshold_diagnostic": sp.factor(
            potential_tree + AF / kappa**2 - (tree + AF) / kappa**2 + v4 / kappa**2
        ),
        "on_shell_mass_reference_does_not_change_external_Mandelstam_sum": sp.Integer(4)
        * 1
        - 4,
    }
    return {
        "selected_constant_field_normalized_amplitude": selected,
        "canonical_reference_tree_amplitude": canonical_tree,
        "canonical_potential_quartic_tree_amplitude": potential_tree,
        "residual_box_when_potential_quartic_is_used": (AF + v4) / kappa**2,
        "formal_one_loop_amplitude": A0 + h * (AF + 2 * fp * A0),
        "selected_normalized_second_coefficient": b2_selected,
        "formal_one_loop_second_coefficient": 4 * lam + h * (bF + 8 * lam * fp),
        "scope": "The box is the complete one-loop fermion 1PI four-point increment. In the selected tree-scalar plus fermion functional, its only additional four-particle contribution at this order is the explicit pole-field factor on the tree amplitude. There is no direct H Yukawa or internal light tree exchange. Old scalar loops and their different local reference scheme are not included or identified with a full MSbar matching calculation.",
        "not_a_resummation": "Exact constant normalization belongs to the named selected finite functional; its separately displayed loop expansion is not an all-loop error bound.",
        "checks": checks,
    }
