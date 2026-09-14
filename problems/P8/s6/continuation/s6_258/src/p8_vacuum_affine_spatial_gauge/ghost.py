"""Whole off-gauge Fourier operator and retained finite ghost vertices."""

from functools import cache

import sympy as s

from .gauge import WEIGHT


def spectral_block(Q, k, ell):
    k, ell = s.Matrix(k), s.Matrix(ell)
    p = k + ell
    return (
        (p.T * Q * k)[0] * s.eye(3)
        + (Q * k) * p.T
        - (Q * p) * ell.T
        - WEIGHT * (Q * p) * k.T
    )


def assemble(modes, coefficients):
    matrix = s.zeros(3 * len(modes))
    for i, p in enumerate(modes):
        for j, k in enumerate(modes):
            ell = tuple(p[t] - k[t] for t in range(3))
            Q = coefficients.get(ell)
            if Q is not None:
                matrix[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = spectral_block(Q, k, ell)
    return matrix


@cache
def full_fourier():
    a, b, c, d, e, f = s.symbols(
        "density11 density12 density13 density22 density23 density33", real=True
    )
    Q = s.Matrix([[a, b, c], [b, d, e], [c, e, f]])
    k = s.Matrix(s.symbols("incoming_ghost_momentum0:3", real=True))
    ell = s.Matrix(s.symbols("coefficient_momentum0:3", real=True))
    xi = s.Matrix(s.symbols("ghost_vector0:3", real=True))
    lie = s.I * (
        ell.dot(xi) * Q - (Q * k) * xi.T - xi * (Q * k).T + WEIGHT * Q * k.dot(xi)
    )
    direct = s.I * lie * (k + ell)
    whole = spectral_block(Q, k, ell)
    constant = spectral_block(Q, (0, 0, 0), ell)
    zero_out = spectral_block(Q, k, -k)
    positive = spectral_block(Q, k, (0, 0, 0))
    return {
        "whole_symmetric_density_Fourier_coefficient": Q,
        "whole_incoming_ghost_momentum": k,
        "whole_coefficient_momentum": ell,
        "whole_off_gauge_Fourier_operator": whole,
        "whole_constant_gauge_vector_off_slice_action": constant,
        "whole_homogeneous_symbol": positive,
        "Fourier_boundary": "The outgoing momentum is p=k+ell and the entire block is (p^T Q k)I+(Qk)p^T-(Qp)ell^T-(2/3)(Qp)k^T. No on-gauge simplification is made before differentiation or projection. The p=0 row vanishes exactly because the gauge condition is a divergence. The k=0 column need not vanish off the slice; it vanishes on div Q=0. Do not drop an off-gauge constant column by claiming a global translation kernel everywhere.",
        "checks": {
            "whole_independent_fourier_Lie_derivative": (direct - whole * xi).applyfunc(
                s.expand
            ),
            "whole_zero_outgoing_row_is_exact_divergence": zero_out.applyfunc(s.expand),
            "whole_constant_incoming_transport_of_gauge_residual": (
                constant + (Q * ell) * ell.T
            ).applyfunc(s.expand),
            "whole_homogeneous_full_symbol": (
                positive
                - (k.T * Q * k)[0] * s.eye(3)
                - s.Rational(1, 3) * (Q * k) * k.T
            ).applyfunc(s.expand),
        },
        "gates": {
            "whole_off_gauge_constant_column_is_nonzero": constant != s.zeros(3),
            "all_three_ghost_components_in_every_block": whole.shape == (3, 3),
            "whole_fourier_operator_contains_coefficient_momentum": all(
                whole.has(item) for item in ell
            ),
        },
    }


@cache
def finite_vertices():
    modes = tuple((a, b, 0) for a in (-1, 0, 1) for b in (-1, 0, 1) if a != 0 or b != 0)
    S = s.diag(0, 1, -1)
    first = {(1, 0, 0): -S / 2, (-1, 0, 0): -S / 2}
    second = {(0, 0, 0): S * S / 2, (2, 0, 0): S * S / 4, (-2, 0, 0): S * S / 4}
    base = {(0, 0, 0): s.eye(3)}
    A, B, C = assemble(modes, base), assemble(modes, first), assemble(modes, second)
    inverse = A.inv()
    first_trace = s.trace(inverse * B)
    second_trace = s.factor(s.trace(inverse * C - inverse * B * inverse * B))
    no_contact = s.factor(-s.trace(inverse * B * inverse * B))
    derivative_inverse = -inverse * B * inverse
    return {
        "explicit_finite_nonzero_Fourier_modes": modes,
        "whole_three_ghost_reference_matrix": A,
        "whole_first_shape_ghost_matrix": B,
        "whole_second_shape_ghost_matrix": C,
        "whole_reference_ghost_inverse": inverse,
        "whole_first_inverse_variation": derivative_inverse,
        "whole_normalized_logdet_first_variation": first_trace,
        "whole_normalized_logdet_second_variation": second_trace,
        "incorrect_second_variation_without_shape_contact": no_contact,
        "entire_reference_determinant": s.factor(A.det()),
        "exact_shape_family": "Use gamma=a_hat^2 exp(2v) exp(t S cos x1), S=diag(0,1,-1). Its entire conformal density is Q=diag(1,exp(-t cos x1),exp(t cos x1)), with det Q=1 and div Q=0 exactly. The displayed finite matrices are the full reference, first and second derivatives of the projected entire operator, not an action truncated to a quadratic shape polynomial. For |t|<=1/100, ||Q-I||op<=|t|/(1-|t|)=1/99<1/8.",
        "finite_determinant_boundary": "This projection has eight nonzero spatial Fourier modes and all three spatial ghosts, hence a 24-by-24 matrix. The relative determinant divides by the explicitly retained reference determinant only. At finite dimension d log det M=Tr(M^-1 dM) and d1 d2 log det M=Tr(M^-1 d1d2 M-M^-1 d1M M^-1 d2M). The second variation here is 3 and is not zero. The omitted shape contact gives the different displayed value. No continuum determinant, time regulator, physical counterterm or full gauge-invariant quantum measure is inferred.",
        "checks": {
            "whole_finite_reference_inverse": A * inverse - s.eye(24),
            "whole_finite_inverse_first_variation": A * derivative_inverse
            + B * inverse,
            "whole_finite_first_logdet_variation": first_trace,
            "whole_finite_second_logdet_variation": second_trace - 3,
            "whole_retained_second_shape_contact_difference": s.factor(
                second_trace - no_contact - s.trace(inverse * C)
            ),
        },
        "gates": {
            "all_24_finite_ghost_coordinates_retained": A.shape
            == B.shape
            == C.shape
            == (24, 24),
            "entire_reference_matrix_regular": A.det() != 0,
            "whole_second_shape_ghost_vertex_nonzero": C != s.zeros(24),
            "omitting_second_shape_contact_changes_logdet_vertex": no_contact
            != second_trace,
            "exact_shape_smallness_inside_coercive_box": s.Rational(1, 99)
            < s.Rational(1, 8),
        },
    }


@cache
def mixed_off_gauge_contact():
    # Two individually TT plane waves whose nonlinear shape contact is not TT.
    ell1 = s.Matrix([1, 0, 0])
    ell2 = s.Matrix([0, 1, 0])
    A = s.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    B = s.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    second = (A * B + B * A) / 8
    ell = ell1 + ell2
    k = s.Matrix([1, -1, 1])
    full = spectral_block(second, k, ell)
    # The gauge-surface-only differential operator wrongly omits all chi terms.
    truncated = (k.T * second * k)[0] * s.eye(3) + s.Rational(1, 3) * (second * k) * k.T
    ell_squared = ell.dot(ell)
    inverse = (s.eye(3) - ell * ell.T / (4 * ell_squared)) / ell_squared
    correction = -s.I * inverse * second * ell
    shape_correction = s.I * (
        ell * correction.T
        + correction * ell.T
        - s.Rational(2, 3) * s.eye(3) * ell.dot(correction)
    )
    restored = second - shape_correction
    volume_correction = s.I * ell.dot(correction) / 3
    restored_block = spectral_block(restored, k, ell)
    restored_surface = (k.T * restored * k)[0] * s.eye(3) + s.Rational(1, 3) * (
        restored * k
    ) * k.T
    return {
        "two_distinct_TT_shape_polarizations": (A, B),
        "two_distinct_shape_wavevectors": (ell1, ell2),
        "whole_mixed_exponential_shape_Fourier_contact": second,
        "whole_mixed_shape_gauge_residual_coefficient": s.I * second * ell,
        "whole_second_order_off_gauge_operator_block": full,
        "incorrect_gauge_surface_only_second_vertex": truncated,
        "whole_second_order_gauge_restoring_vector": correction,
        "whole_second_order_metric_shape_correction": shape_correction,
        "whole_gauge_restored_mixed_density_coefficient": restored,
        "whole_second_order_scalar_volume_contact": volume_correction,
        "whole_gauge_restored_second_ghost_vertex": restored_block,
        "nonlinear_gauge_boundary": "Each first shape wave is TT in its own direction, but their mixed exp(-t) contact has a nonzero divergence. The displayed second-order coordinate correction restores the gauge and has a nonzero scalar-volume contact 1/32 for this coefficient mode. Its restored full ghost block agrees with the surface formula only after that correction. Thus neither the raw mixed family nor its gauge-surface-only ghost vertex is silently used as the nonlinear reference chart. Keep this coordinate change in all physical source and canonical boundary maps; an off-shell reference can retain one-point contacts, as S257 showed.",
        "checks": {
            "first_shape_wave_is_transverse": A * ell1,
            "second_shape_wave_is_transverse": B * ell2,
            "both_shape_waves_are_tracefree": s.Matrix([s.trace(A), s.trace(B)]),
            "whole_mixed_exponential_shape_contact": second - (A * B + B * A) / 8,
            "whole_second_order_gauge_restoration": (restored * ell).applyfunc(
                s.factor
            ),
            "whole_restored_full_and_surface_ghost_vertices": (
                restored_block - restored_surface
            ).applyfunc(s.factor),
            "whole_nonzero_second_order_scalar_volume_contact": volume_correction
            - s.Rational(1, 32),
        },
        "gates": {
            "mixed_shape_contact_is_not_transverse": second * ell != s.zeros(3, 1),
            "off_gauge_second_vertex_differs_from_surface_only_vertex": full
            != truncated,
            "nonlinear_gauge_correction_not_assumed_zero": True,
            "restoring_coordinate_change_has_scalar_volume_contact": volume_correction
            != 0,
        },
    }
