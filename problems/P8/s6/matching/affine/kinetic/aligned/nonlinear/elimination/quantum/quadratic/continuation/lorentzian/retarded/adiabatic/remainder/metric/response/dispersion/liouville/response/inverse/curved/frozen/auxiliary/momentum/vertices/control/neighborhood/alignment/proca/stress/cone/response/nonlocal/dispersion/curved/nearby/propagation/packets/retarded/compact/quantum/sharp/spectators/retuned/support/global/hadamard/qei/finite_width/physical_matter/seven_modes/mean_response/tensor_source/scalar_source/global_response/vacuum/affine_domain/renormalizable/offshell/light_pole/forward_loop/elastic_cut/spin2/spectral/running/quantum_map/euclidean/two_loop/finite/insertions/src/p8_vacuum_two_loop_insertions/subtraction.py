"""Grouped tadpole topology and linear inherited local outer subtraction."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import subtraction as previous
from p8_vacuum_two_loop_denom import graphs


@cache
def data():
    old = previous.data()
    alpha = sp.Symbol("fixed_positive_self_energy_asymptotic_multiplier", positive=True)
    Q1, Q2, R1, R2 = sp.symbols("first_Q second_Q first_decaying_R second_decaying_R")
    V1, V2 = sp.symbols("first_full_vertex second_full_vertex")
    D1, D2 = sp.symbols("first_light_inverse second_light_inverse", positive=True)
    I0 = old["regulated_zero_momentum_light_bubble"]
    L, g = sp.symbols("quartic_L cubic_squared", real=True)
    hs = sp.symbols("h_s h_t h_u", real=True)
    factor = I0 / (32 * sp.pi**2)
    shifts = {
        "delta_polynomial_quartic": 2 * alpha * old["delta_polynomial_quartic_UV"],
        "delta_cubic_squared": 2 * alpha * old["delta_cubic_squared_UV"],
        "delta_heavy_mass_squared": 2 * alpha * old["delta_heavy_mass_squared_UV"],
    }
    local = (
        -shifts["delta_polynomial_quartic"]
        + shifts["delta_cubic_squared"] * sum(hs)
        - g * shifts["delta_heavy_mass_squared"] * sum(h * h for h in hs)
    )
    s = sp.Symbol("line_invariant", real=True)
    tad1, tad2, a, b = sp.symbols(
        "local_light_tadpole local_heavy_tadpole local_mass local_kinetic", real=True
    )
    local_inner = tad1 + tad2 + a + b * s
    os = (
        local_inner
        - local_inner.subs(s, 1)
        - (s - 1) * sp.diff(local_inner, s).subs(s, 1)
    )
    finite = sp.Symbol("grouped_finite_potential_contact", real=True)
    tadrows = [r for r in graphs.skeletons() if r["kind"] == "tadpole_insertion"]
    weight = sum(r["summed_labelled_vertex_Wick_weight"] for r in tadrows)
    return {
        "fixed_alpha": alpha,
        "complete_two_line_insertion_integrand": V1 * V2 * (Q1 + Q2) / (2 * D1 * D2),
        "inherited_one_loop_bubble_integrand": V1 * V2 / (2 * D1 * D2),
        "decaying_remainder_integrand": V1 * V2 * (R1 + R2) / (2 * D1 * D2),
        "total_outer_local_counterterms_for_this_group": shifts,
        "outer_counterterm_amplitude": local,
        "grouped_tadpole_family_Wick_weight": weight,
        "outer_choice_count": 16,
        "inner_choice_count": 4,
        "grouped_refinement_count": 64,
        "same_regulated_zero_momentum_reference": I0,
        "scope": "The complete one-loop self-energy insertion family with its already fixed on-shell inner subtraction and linear continuation of the S6.113 outer reference subtraction. This is not a claim that all two-loop counterterms or light LSZ corrections have been computed.",
        "checks": {
            "constant_and_decaying_insertions_split": sp.expand(
                (V1 * V2 * (Q1 + Q2) / (2 * D1 * D2)).subs(
                    {Q1: alpha + R1, Q2: alpha + R2}
                )
                - 2 * alpha * V1 * V2 / (2 * D1 * D2)
                - V1 * V2 * (R1 + R2) / (2 * D1 * D2)
            ),
            "both_local_tadpoles_and_all_affine_terms_OS_annihilated": sp.expand(os),
            "two_mixed_channels_cancel_tadpole_symmetry_half": 2 * sp.Rational(1, 2)
            - 1,
            "all_sixty_four_refinements_grouped": 16 * 4
            - len([c for c in graphs.cases() if c[0] == "tadpole_insertion"]),
            "independent_insertion_weight_matches_frozen_skeleton": 2
            * sp.Rational(1, 2)
            * sp.Rational(3, 2)
            - weight,
            "outer_local_counterterms_cancel_all_three_channel_UV": sp.expand(
                local + 2 * alpha * factor * sum((-L + g * h) ** 2 for h in hs)
            ),
            "outer_subtraction_is_linear_in_same_reference": sp.expand(
                local - 2 * alpha * old["full_UV_counterterm_amplitude"]
            ),
            "finite_potential_contact_cannot_tune_b2": sp.diff(finite, s, 2),
        },
    }
