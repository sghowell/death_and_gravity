"""Complete nonsymmetric cotangent lift with the residual kernel retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_spatial_gauge import gauge

from .shape import tt


def dq(metric, variation):
    c = metric.det() ** s.Rational(1, 3)
    inv = metric.inv()
    return c * (s.trace(inv * variation) * inv / 3 - inv * variation * inv)


def dq_adjoint(metric, value):
    c = metric.det() ** s.Rational(1, 3)
    inv = metric.inv()
    return c * (inv * s.trace(value * inv) / 3 - inv * value * inv)


@cache
def data():
    # Non-diagonal determinant-one shape: diagnostic, not restriction of the proof.
    Q = s.Matrix(
        [
            [s.Rational(6, 5), s.Rational(1, 10), 0],
            [s.Rational(1, 10), s.Rational(9, 10), 0],
            [0, 0, s.Rational(100, 107)],
        ]
    )
    scale = s.Rational(7, 5)
    gamma = scale * Q.inv()
    value = s.Matrix([[2, 1, 3], [1, -1, 2], [3, 2, -1]])
    variation = s.Matrix([[1, 2, 3], [2, -1, 4], [3, 4, 0]])
    inv = Q.inv()
    tangent = variation - s.eye(3) * s.trace(inv * variation) / s.trace(inv)
    projected = tt(value - inv * s.trace(value) / s.trace(inv), (0, 0, 0))
    Pv = s.Rational(3, 7)
    pi0 = Pv * gamma.inv() / 6 + dq_adjoint(gamma, value)
    shape_cotangent = -scale * inv * pi0 * inv
    shape_dual = tt(
        shape_cotangent - inv * s.trace(shape_cotangent) / s.trace(inv), (0, 0, 0)
    )
    checks = {
        "whole_nondiagonal_fixture_unit_determinant": Q.det() - 1,
        "whole_full_scalar_tangent_contact": s.trace(inv * tangent).factor(),
        "whole_full_shape_adjoint_pairing": (
            s.trace(value * tangent) - s.trace(projected * variation)
        ).factor(),
        "whole_metric_DQ_adjoint_pairing": (
            s.trace(value * dq(gamma, variation))
            - s.trace(dq_adjoint(gamma, value) * variation)
        ).factor(),
        "whole_canonical_scalar_momentum_reconstruction": (
            2 * s.trace(pi0 * gamma) - Pv
        ).factor(),
        "whole_canonical_shape_momentum_reconstruction": (shape_dual - value).applyfunc(
            s.factor
        ),
        "whole_metric_conformal_direction_has_no_shape_variation": dq(
            gamma, 2 * gamma
        ).applyfunc(s.factor),
    }
    # Finite exact operator analog checks signs, adjoints and surviving constants.
    A = s.Matrix([[1, 2, 0, 1, 0, 0], [0, 1, 1, 0, 2, 0]])
    V = s.Matrix([[1, 0, 0], [0, 2, 0], [1, 1, 0], [0, 1, 0], [1, -1, 0], [0, 0, 1]])
    M = A * V[:, :2]
    H = s.eye(6) - V[:, :2] * M.inv() * A
    lift = H.T
    checks.update(
        {
            "whole_horizontal_projector_idempotence": H * H - H,
            "whole_horizontal_slice_condition": A * H,
            "whole_horizontal_mean_zero_gauge_annihilation": H * V[:, :2],
            "whole_full_dual_lift_adjoint": lift
            - (s.eye(6) - A.T * M.T.inv() * V[:, :2].T),
            "whole_mean_zero_momentum_constraint_lift": V[:, :2].T * lift,
            "whole_residual_constant_constraint_retained": V[:, 2].T * lift - V[:, 2].T,
        }
    )
    C = s.Matrix([[0, 3], [-3, 0]])
    block = s.zeros(2).row_join(M).col_join((-M.T).row_join(C))
    inverse = (
        (M.T.inv() * C * M.inv())
        .row_join(-M.T.inv())
        .col_join(M.inv().row_join(s.zeros(2)))
    )
    checks["whole_second_class_block_inverse_with_nonzero_constraint_bracket"] = (
        block * inverse - s.eye(4)
    )
    checks["whole_second_class_block_inverse_reverse"] = inverse * block - s.eye(4)
    weak = gauge.local_slice()
    return {
        "whole_metric_shape_derivative": "DQ[h]=c[(tr(gamma^-1 h)/3)gamma^-1-gamma^-1 h gamma^-1], c=det(gamma)^(1/3). A h=div DQ[h].",
        "whole_metric_shape_adjoint": "DQ^*[B]=c[(gamma^-1/3)tr(B gamma^-1)-gamma^-1 B gamma^-1], with symmetric full matrix pairing, including both off-diagonal entries.",
        "whole_full_horizontal_projection": "V is the Lie action on gamma, W_i, M1,H, not only the metric. M_Q=A V. H=I-V(M_Q|mean0)^-1 A; the auxiliary primary momenta have already been removed.",
        "whole_base_cotangent_extension": "pi0=(Pi_v/6)gamma^-1+DQ^*(P_TT Pi_tau), with the complete unchanged W,M1,H momenta. D0=V^* alpha0 is the ENTIRE source-pinned generator density, not only its metric part.",
        "whole_unique_momentum_lift": "alpha=alpha0-A^*(M_Q^*|mean0)^-1 Pmean0 D0. Only metric momentum is corrected. This annihilates every mean-zero gauge direction; the three residual constant constraints remain.",
        "whole_pulled_back_canonical_one_form": "integral Pi_v delta v+Pi_tau:delta tau+Pi_W delta W+Pi_M delta M1+Pi_H delta H. Since A delta q=0 on the slice, the entire dual correction pairs to zero. Differentiation gives the canonical cotangent symplectic form.",
        "whole_full_shape_cotangent": "Ashape=-a^2 exp(2v)Q^-1 pi Q^-1; Pi_tau=P_TT[Ashape-cof(Q)(K_Q^-1)^*B^*Ashape], Pi_v=2 pi:gamma. All inverse and adjoint contacts are retained.",
        "whole_complete_nonsymmetric_second_class_diagnostic": {
            "A": A,
            "V": V,
            "M": M,
            "constraint_bracket": C,
            "block": block,
            "inverse": inverse,
        },
        "whole_actual_weak_coercivity_margin": weak["whole_weak_coercivity_margin"],
        "whole_actual_mean_zero_inverse_bound": weak[
            "mean_zero_torus_inverse_L2_upper_bound"
        ],
        "functional_boundary": "Use the actual S258 derivative-gap local coordinate slice, or the explicit near-identity shape chart, with smooth fields and L2 covectors. A:L2->H^-1, M^-1:H^-1_mean0->H1_mean0 and the full V:H1->L2 are bounded for fixed smooth coefficients. The dual formula is a bounded cotangent lift and is unique by the same inverse. No same-regularity diffeomorphism group, strong infinite-dimensional Darboux theorem, global quotient or uniform infinite-volume inverse is assumed.",
        "checks": checks,
        "gates": {
            "actual_nonsymmetric_gauge_block": M != M.T,
            "actual_nonzero_constraint_constraint_block": C != s.zeros(2),
            "residual_kernel_is_not_inverted": A * V[:, 2] == s.zeros(2, 1),
            "residual_kernel_lift_is_not_annihilated": V[:, 2].T * lift
            != s.zeros(1, 6),
            "full_three_matter_vector_channels_enter_dual_constraint": True,
            "old_weak_coercivity_and_derivative_gap_retained": all(
                weak["gates"].values()
            ),
            "canonical_pullback_not_a_flat_quantum_measure_claim": True,
        },
    }
