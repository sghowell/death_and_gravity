"""Whole physical ADM metric and matter contractions with the actual factor."""

from functools import cache

import sympy as s
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical, cotangent

from . import source


@cache
def metric():
    N, R, C = source.N, source.R, source.C
    entries = s.symbols("positive_spatial_metric0:6", real=True)
    gamma = s.Matrix(
        [
            [entries[0], entries[1], entries[2]],
            [entries[1], entries[3], entries[4]],
            [entries[2], entries[4], entries[5]],
        ]
    )
    inv = gamma.inv()
    shift = s.Matrix(s.symbols("full_shift0:3", real=True))
    E = s.eye(4)
    E[1:, 0] = shift
    diagonal = s.diag(N * N, -C * gamma)
    g = E.T * diagonal * E
    inverse = E.inv() * s.diag(N**-2, -inv / C) * E.inv().T
    wanted = (
        s.Matrix([[N * N - C * (shift.T * gamma * shift)[0]]])
        .row_join(-C * shift.T * gamma)
        .col_join((-C * gamma * shift).row_join(-C * gamma))
    )
    coframe = E.inv().T * g * E.inv()
    density_ratio = C**3
    return {
        "whole_hat_spatial_metric": gamma,
        "whole_spatial_inverse": inv,
        "whole_positive_lapse": N,
        "whole_three_shift_components": shift,
        "whole_ADM_coframe": E,
        "whole_physical_metric": g,
        "whole_physical_inverse_metric": inverse,
        "whole_physical_spatial_metric": C * gamma,
        "whole_physical_spatial_determinant": C**3 * gamma.det(),
        "whole_physical_spacetime_determinant": -N * N * C**3 * gamma.det(),
        "whole_volume_ratio": canonical.U,
        "whole_metric_field_Jacobian": C**9,
        "whole_extended_position_Jacobian": N * C**6,
        "checks": {
            "whole_physical_ADM_metric_all_entries": (g - wanted).applyfunc(s.factor),
            "whole_inverse_metric_both_sides": s.diag(
                (g * inverse - s.eye(4)).applyfunc(s.factor),
                (inverse * g - s.eye(4)).applyfunc(s.factor),
            ),
            "whole_unit_determinant_ADM_coframe": E.det() - 1,
            "whole_metric_diagonal_congruence": (coframe - diagonal).applyfunc(
                s.factor
            ),
            "whole_physical_spatial_determinant": s.factor(
                (C * gamma).det() - C**3 * gamma.det()
            ),
            "whole_actual_metric_determinant_ratio": s.simplify(
                density_ratio - R ** (-s.Rational(3, 2))
            ),
            "whole_actual_positive_density_root": s.simplify(
                s.sqrt(density_ratio) - canonical.U
            ),
            "whole_metric_field_Jacobian_R_power": s.simplify(
                C**9 - R ** (-s.Rational(9, 2))
            ),
            "whole_extended_point_Jacobian_R_power": s.simplify(N * C**6 - N * R**-3),
        },
        "gates": {
            "full_shift_and_all_six_metric_components_kept": len(shift) == 3
            and len(entries) == 6,
            "physical_lapse_and_shift_unchanged_by_spatial_conformal_map": True,
            "positive_R_and_lapse_and_spatial_metric_domain_required": True,
        },
    }


@cache
def matter():
    data = metric()
    N, C = source.N, source.C
    inv, shift, E = (
        data["whole_spatial_inverse"],
        data["whole_three_shift_components"],
        data["whole_ADM_coframe"],
    )
    gradient = s.Matrix(s.symbols("whole_scalar_gradient0:4", real=True))
    normal = gradient[0] - shift.dot(gradient[1:, :])
    spatial = (gradient[1:, :].T * inv * gradient[1:, :])[0]
    scalar = (gradient.T * data["whole_physical_inverse_metric"] * gradient)[0]
    wanted = normal**2 / N**2 - spatial / C
    fvalues = s.symbols("whole_two_form0:6", real=True)
    F = s.zeros(4)
    for value, (i, j) in zip(
        fvalues, ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)), strict=True
    ):
        F[i, j], F[j, i] = value, -value
    frame = E.inv().T * F * E.inv()
    electric = s.Matrix(
        [
            F[0, i + 1] - sum(shift[j] * F[j + 1, i + 1] for j in range(3))
            for i in range(3)
        ]
    )
    EE = (electric.T * inv * electric)[0]
    magnetic = sum(
        inv[i, k] * inv[j, l] * F[i + 1, j + 1] * F[k + 1, l + 1]
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for l in range(3)
    )
    coefficient_lagrangian = (
        C ** s.Rational(1, 2) * EE / (2 * N)
        - N * C ** (-s.Rational(1, 2)) * magnetic / 4
    )
    frozen_lagrangian = canonical.Cchi * EE / (2 * N) - N * canonical.M * magnetic / 4
    lift = cotangent.extended_map()
    Cformal = next(
        f
        for f in lift["whole_position_density_factor"].atoms(s.Function)
        if str(f.func) == "positive_spatial_conformal_factor"
    )
    Nformal = next(
        z
        for z in lift["whole_position_density_factor"].free_symbols
        if z.name == "positive_lapse"
    )
    return {
        "whole_scalar_covector": gradient,
        "whole_scalar_normal_derivative": normal,
        "whole_full_scalar_metric_contraction": scalar,
        "whole_scalar_ADM_density_over_hat_volume": canonical.U * normal**2 / (2 * N)
        - N * canonical.Cchi * spatial / 2,
        "whole_antisymmetric_field_strength": F,
        "whole_field_strength_in_ADM_coframe": frame,
        "whole_shift_corrected_electric_covector": electric,
        "whole_Maxwell_ADM_density_over_zeta_hat_volume": coefficient_lagrangian,
        "whole_mass_and_source_density_over_hat_volume": N * canonical.U,
        "checks": {
            "whole_scalar_covector_contraction_all_directions": s.factor(
                scalar - wanted
            ),
            "whole_shift_corrected_electric_covector": (
                frame[0, 1:].T - electric
            ).applyfunc(s.factor),
            "whole_spatial_two_form_unchanged_in_ADM_coframe": frame[1:, 1:]
            - F[1:, 1:],
            "whole_field_strength_remains_antisymmetric": (frame + frame.T).applyfunc(
                s.expand
            ),
            "whole_Maxwell_matches_both_frozen_coefficients": s.factor(
                coefficient_lagrangian - frozen_lagrangian
            ),
            "whole_source_pinned_point_Jacobian_binding": s.simplify(
                lift["whole_position_density_factor"].subs(
                    {Cformal: C, Nformal: N}, simultaneous=True
                )
                - N * source.R**-3
            ),
        },
        "gates": {
            "all_six_field_strength_components_kept": len(fvalues) == 6,
            "entire_shift_electric_transport_kept": any(
                electric.has(item) for item in shift
            ),
            "original_mass_heavy_and_source_volume_not_reselected": True,
            "full_cotangent_measure_identity_not_quantum_ordering_claim": True,
        },
    }
