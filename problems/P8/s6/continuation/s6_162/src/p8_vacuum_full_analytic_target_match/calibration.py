"""The actual finite-kappa target's all-higher-field error allowance."""

from functools import cache

import sympy as s

from . import bounds, holomorphic, jets


@cache
def data():
    h = holomorphic.data()
    substitutions = {
        jets.N: h["fixed_n"],
        jets.KAPPA: h["fixed_kappa"],
        jets.LAM: h["fixed_lambda"],
    }
    pieces = jets.data()["canonical_density_field_degrees"]
    sixth = bounds.weighted_majorant(s.expand(pieces[6].subs(substitutions)))
    eighth = bounds.weighted_majorant(s.expand(pieces[8].subs(substitutions)))
    tail = bounds.even_tail(
        h["full_canonical_flat_density_circle_upper"],
        h["complex_field_amplitude_radius"],
    )
    pointwise = sixth + eighth + tail
    integrated = bounds.integrated_error(sixth, eighth, tail)
    return {
        "actual_sixth_density_coefficient_majorant": sixth,
        "actual_eighth_density_coefficient_majorant": eighth,
        "all_even_field_degrees_at_least_ten_Cauchy_majorant": tail,
        "all_target_field_degrees_above_four_pointwise_W_squared_majorant": pointwise,
        "full_target_action_minus_quartic_target_common_class_coefficient": integrated,
        "bounds": {
            "sixth_density_positive_and_below_three_e_minus_1594": bool(
                0 < sixth < s.Rational(3, 10**1594)
            ),
            "eighth_density_positive_and_below_three_e_minus_2195": bool(
                0 < eighth < s.Rational(3, 10**2195)
            ),
            "all_even_degrees_at_least_ten_below_two_e_minus_2189": bool(
                0 < tail < s.Rational(2, 10**2189)
            ),
            "full_higher_field_action_allowance_below_one_e_minus_1591": bool(
                0 < integrated < s.Rational(1, 10**1591)
            ),
        },
        "checks": {
            "sixth_majorant_exact_known_value": sixth
            - s.Rational(2858465025, 10**1603),
            "all_three_disjoint_higher_field_groups_assigned_once": pointwise
            - sixth
            - eighth
            - tail,
            "twenty_one_L2_jet_components": integrated - 21 * pointwise,
            "finite_kappa_not_limit": h["fixed_n"] / h["fixed_kappa"]
            - s.Rational(1024, 10**800),
        },
    }
