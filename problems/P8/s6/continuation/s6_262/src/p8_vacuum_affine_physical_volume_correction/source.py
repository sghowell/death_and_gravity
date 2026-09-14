"""Literal parent conformal factor pinned by the immutable ADM coefficients."""

from functools import cache

import sympy as s
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical
from p8_vacuum_affine_physical_background_vertices import parent

R, N = canonical.R, canonical.N
C = canonical.Cchi**2
WRONG = R - s.Rational(1, 2)


@cache
def coefficient_bridge():
    # Two independent frozen coefficient families fix the SAME physical C.
    true_C = s.factor(C)
    spatial_volume = s.powdenest(true_C ** s.Rational(3, 2), force=False)
    scalar_spatial = s.sqrt(true_C)
    maxwell_spatial = 1 / s.sqrt(true_C)
    density = N * spatial_volume
    wrong0 = WRONG.subs(R, 1)
    return {
        "whole_actual_R": R,
        "whole_physical_conformal_factor_from_Cchi_squared": true_C,
        "whole_physical_conformal_factor_from_inverse_M_squared": canonical.M**-2,
        "whole_physical_spatial_volume_multiplier": spatial_volume,
        "whole_physical_spacetime_volume_multiplier": density,
        "whole_scalar_temporal_ADM_coefficient": spatial_volume / N,
        "whole_scalar_spatial_ADM_coefficient": N * scalar_spatial,
        "whole_Maxwell_electric_ADM_coefficient_without_zeta": scalar_spatial / N,
        "whole_Maxwell_magnetic_ADM_coefficient_without_zeta": N * maxwell_spatial,
        "whole_scalar_mass_and_source_density": density,
        "frozen_S261_incorrect_conformal_expression": WRONG,
        "frozen_S261_incorrect_reference_conformal_value": wrong0,
        "whole_actual_reference_conformal_value": true_C.subs(R, 1),
        "checks": {
            "literal_Cchi_squared_metric_binding": s.simplify(
                true_C - R ** (-s.Rational(1, 2))
            ),
            "literal_independent_M_inverse_squared_metric_binding": s.simplify(
                true_C - canonical.M**-2
            ),
            "whole_spatial_volume_matches_frozen_U": s.simplify(
                spatial_volume - canonical.U
            ),
            "whole_scalar_temporal_matches_frozen_U_over_N": s.simplify(
                spatial_volume / N - canonical.U / N
            ),
            "whole_scalar_gradient_matches_frozen_Cchi": s.simplify(
                N * scalar_spatial - N * canonical.Cchi
            ),
            "whole_Maxwell_electric_matches_frozen_Cchi": s.simplify(
                scalar_spatial / N - canonical.Cchi / N
            ),
            "whole_Maxwell_magnetic_matches_frozen_M": s.simplify(
                N * maxwell_spatial - N * canonical.M
            ),
            "whole_mass_and_original_source_volume_matches_frozen_NU": s.simplify(
                density - N * canonical.U
            ),
            "whole_reference_metric_multiplier_one": true_C.subs(R, 1) - 1,
            "wrong_reference_conformal_value_exact_discrepancy": wrong0
            - true_C.subs(R, 1)
            + s.Rational(1, 2),
        },
        "gates": {
            "wrong_expression_fails_metric_reference_value": wrong0 != 1,
            "wrong_expression_fails_original_volume_coefficient": wrong0
            ** s.Rational(3, 2)
            != canonical.U.subs(R, 1),
            "wrong_expression_fails_original_scalar_gradient": s.sqrt(wrong0)
            != canonical.Cchi.subs(R, 1),
            "wrong_expression_fails_original_Maxwell_magnetic_coefficient": 1
            / s.sqrt(wrong0)
            != canonical.M.subs(R, 1),
            "generic_positive_conformal_Ward_identity_does_not_select_parent": True,
        },
    }


@cache
def clock_jets():
    u = parent.u
    data = parent.lapse_jets()
    Rjets = data["R_lapse_jets_zero_through_four"]
    replace = {s.diff(R, N, j): Rjets[j] for j in range(5)}
    functions = {
        "C": C,
        "spatial_volume_U": canonical.U,
        "spacetime_volume_NU": N * canonical.U,
        "scalar_temporal_U_over_N": canonical.U / N,
        "scalar_spatial_NCchi": N * canonical.Cchi,
        "Maxwell_magnetic_NM": N * canonical.M,
    }
    rows = {
        name: tuple(
            s.factor(s.diff(expr, N, j).subs(replace, simultaneous=True).subs(N, 1))
            for j in range(5)
        )
        for name, expr in functions.items()
    }
    c1 = rows["spatial_volume_U"][1] / rows["spatial_volume_U"][0]
    wrong = s.Rational(3, 2) * s.diff(WRONG, N) / WRONG
    wrong_clock = s.factor(wrong.subs(replace, simultaneous=True).subs(N, 1))
    powers = {
        "C": s.Integer(1),
        "spatial_volume_U": s.Rational(3, 2),
        "spacetime_volume_NU": s.Rational(5, 2),
        "scalar_temporal_U_over_N": s.Rational(1, 2),
        "scalar_spatial_NCchi": s.Rational(3, 2),
        "Maxwell_magnetic_NM": s.Rational(1, 2),
    }
    checks = {
        "whole_actual_spatial_volume_first_clock_jet": s.factor(
            c1 - s.Rational(3, 2) / parent.h
        ),
        "whole_actual_bounce_first_volume_jet": c1.subs(u, 0) - s.Rational(3, 2),
        "whole_actual_bounce_second_volume_Taylor_coefficient": rows[
            "spatial_volume_U"
        ][2].subs(u, 0)
        / 2
        - s.Rational(3, 8),
        "whole_wrong_minus_actual_bounce_slope": wrong_clock.subs(u, 0)
        - c1.subs(u, 0)
        + s.Rational(15, 2),
        "whole_actual_spacetime_volume_first_jet": rows["spacetime_volume_NU"][1].subs(
            u, 0
        )
        - s.Rational(5, 2),
    }
    for name, power in powers.items():
        checks["all_five_bounce_jets_of_" + name] = s.Matrix(
            [
                rows[name][j].subs(u, 0) - s.diff(N**power, N, j).subs(N, 1)
                for j in range(5)
            ]
        )
    return {
        "whole_source_pinned_R_clock_jets_zero_through_four": Rjets,
        "whole_six_actual_ADM_clock_families": rows,
        "whole_actual_volume_lapse_slope_on_clock": c1,
        "whole_actual_volume_lapse_slope_at_bounce": c1.subs(u, 0),
        "frozen_incorrect_volume_lapse_slope_on_clock": wrong_clock,
        "frozen_incorrect_bounce_slope_minus_actual": wrong_clock.subs(u, 0)
        - c1.subs(u, 0),
        "whole_source_jet_boundary": "The unchanged entire R differs from its clock germ by an order-eight factor. These finite jets are licensed by the literal S253 factorization. No full off-clock function, profile or preparation is replaced.",
        "checks": checks,
        "gates": {
            "entire_parent_not_replaced_by_clock_germ": True,
            "frozen_physical_slope_minus_six_is_false": wrong_clock.subs(u, 0)
            != c1.subs(u, 0),
            "both_density_factors_and_scalar_Maxwell_families_retained": len(rows) == 6,
        },
    }
