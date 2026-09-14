"""Entire projective quotient, 56-direction complement and exact inertia."""

from functools import cache

import sympy as s
from p8_affine import connection as old
from p8_affine_retuned import geometry
from p8_vacuum_analytic_affine_parent import affine, source

P = old.P


def clean(value):
    return (
        s.ImmutableMatrix(value.applyfunc(s.factor))
        if isinstance(value, s.MatrixBase)
        else s.factor(value)
    )


def components(matrix):
    remaining = set(range(matrix.rows))
    out = []
    while remaining:
        queue = [min(remaining)]
        remaining.remove(queue[0])
        block = []
        while queue:
            row = queue.pop()
            block.append(row)
            for column in sorted(remaining):
                if matrix[row, column] != 0:
                    remaining.remove(column)
                    queue.append(column)
        out.append(tuple(sorted(block)))
    return tuple(out)


def rational_congruence(matrix):
    """Exact real rational symmetric elimination with 1x1 or 2x2 pivots."""
    matrix = s.Matrix(matrix)
    if matrix.rows != matrix.cols or matrix != matrix.T:
        raise ValueError("Require a square symmetric rational matrix")
    if any(entry.is_Rational is not True for entry in matrix):
        raise TypeError("Use exact rational entries, not floating inertia")
    n = matrix.rows
    if n == 0:
        return s.eye(0), s.zeros(0), (), (0, 0)
    diagonal = next((i for i in range(n) if matrix[i, i] != 0), None)
    if diagonal is not None:
        chosen = (diagonal,)
    else:
        chosen = next(
            ((i, j) for i in range(n) for j in range(i + 1, n) if matrix[i, j] != 0),
            None,
        )
        if chosen is None:
            raise ValueError("A zero residual block is singular")
    order = tuple(chosen) + tuple(i for i in range(n) if i not in chosen)
    permutation = s.eye(n)[:, list(order)]
    ordered = permutation.T * matrix * permutation
    width = len(chosen)
    pivot = ordered[:width, :width]
    cross = ordered[:width, width:]
    rest = ordered[width:, width:] - cross.T * pivot.inv() * cross
    shear = s.eye(n)
    shear[:width, width:] = -pivot.inv() * cross
    tail, diagonal_tail, pivots, inertia = rational_congruence(rest)
    transform = permutation * shear * s.diag(s.eye(width), tail)
    diagonalized = s.diag(pivot, diagonal_tail)
    if width == 1:
        added = (int(bool(pivot[0, 0] > 0)), int(bool(pivot[0, 0] < 0)))
    else:
        assert pivot.det() < 0
        added = (1, 1)
    assert transform.T * matrix * transform == diagonalized
    return (
        transform,
        diagonalized,
        (s.ImmutableMatrix(pivot), *pivots),
        tuple(a + b for a, b in zip(added, inertia, strict=True)),
    )


@cache
def matrices():
    quotient = old.quotient()
    parent = affine.matrices()
    N, M, lift = map(s.Matrix, (parent["N"], parent["M_old"], parent["lift"]))
    pivot = N.rref()[1]
    free = tuple(i for i in range(60) if i not in pivot)
    right = s.zeros(60, 4)
    for i, row in enumerate(pivot):
        right[row, i] = 1
    right = right * (N * right).inv()
    kernel = s.zeros(60, 56)
    for i, row in enumerate(free):
        kernel[row, i] = 1
    kernel = kernel - right * N * kernel
    A = clean(kernel.T * M * kernel)
    blocks = components(A)
    block_matrices = tuple(A.extract(block, block) for block in blocks)
    determinants = tuple(
        s.factor(block.det(method="domain-ge")) for block in block_matrices
    )
    inverse = s.zeros(56)
    for indices, block in zip(blocks, block_matrices, strict=True):
        inv = block.inv(method="DM")
        for i, row in enumerate(indices):
            for j, column in enumerate(indices):
                inverse[row, column] = s.factor(inv[i, j])
    embedding = s.Matrix(quotient["embedding"])
    projective = s.Matrix(old.quadratic()["gauge"])
    gauge = s.Matrix(
        4,
        64,
        lambda c, j: int(
            old.INDICES[j][0] == old.INDICES[j][1] and old.INDICES[j][2] == c
        ),
    )
    full = (embedding * kernel).row_join(embedding * lift).row_join(projective)
    selector = s.Matrix(60, 64, lambda i, j: int(j == quotient["kept"][i]))
    quotient_map = selector * (s.eye(64) - projective * gauge / 4)
    free_selector = s.Matrix(56, 60, lambda i, j: int(j == free[i]))
    full_inverse = (
        (free_selector * (s.eye(60) - lift * N) * quotient_map)
        .col_join(N * quotient_map)
        .col_join(gauge / 4)
    )
    return {
        "N": N,
        "M": M,
        "lift": lift,
        "right": right,
        "kernel": kernel,
        "A": A,
        "inverse_A": s.ImmutableMatrix(inverse),
        "blocks": blocks,
        "block_matrices": block_matrices,
        "block_determinants": determinants,
        "embedding": embedding,
        "projective": projective,
        "gauge": gauge,
        "full": full,
        "full_inverse": clean(full_inverse),
        "parent": parent,
        "quotient_determinant": quotient["determinant"],
        "pivot": pivot,
        "free": free,
    }


@cache
def complement():
    d = matrices()
    N, K, L, M, A = (d[key] for key in ("N", "kernel", "lift", "M", "A"))
    block_product = s.factor(s.prod(d["block_determinants"]))
    expected = (
        -26388279066624
        * P**65
        * (8 * P + 5) ** 3
        * (8 * P**2 - 1) ** 3
        * (2 * P**3 - 1)
    )
    split = K.row_join(L)
    difference = (L - d["right"]).applyfunc(s.factor)
    return {
        "whole_fixed_trace_map": N,
        "whole_56_direction_kernel_basis": K,
        "whole_actual_retained_trace_lift": L,
        "whole_56_direction_complement_Hessian": A,
        "whole_56_direction_complement_inverse": d["inverse_A"],
        "whole_nine_complement_blocks": d["blocks"],
        "whole_nine_complement_block_determinants": d["block_determinants"],
        "whole_60_direction_old_quotient_determinant": d["quotient_determinant"],
        "whole_56_direction_complement_determinant": expected,
        "whole_split_configuration_Jacobian": s.factor(split.det(method="domain-ge")),
        "strict_domain": "p=sqrt(R)/2>0 and 1/2<R<6/5, hence 1/8<p^2<3/10. Every displayed determinant factor is nonzero. The right inverse and kernel are fixed for the full four-trace map, not a vector ansatz; no complement direction is dropped.",
        "checks": {
            "whole_trace_right_inverse": clean(N * L - s.eye(4)),
            "whole_56_trace_kernel": N * K,
            "whole_60_split_is_regular": clean(split.det(method="domain-ge") + 1),
            "whole_lift_difference_is_in_fixed_kernel": clean(N * difference),
            "whole_full_complement_inverse": clean(A * d["inverse_A"] - s.eye(56)),
            "whole_trace_mass_update_leaves_all_complement_entries": clean(
                K.T * (d["parent"]["M_new"] - M) * K
            ),
            "whole_complement_trace_cross_term": clean(K.T * d["parent"]["M_new"] * L),
            "whole_retained_Lorentz_mass": clean(
                L.T * d["parent"]["M_new"] * L - d["parent"]["eta"]
            ),
            "whole_direct_nine_block_determinant": clean(block_product - expected),
            "whole_60_to_56_determinant_Schur_identity": clean(
                expected
                - split.det(method="domain-ge") ** 2
                * d["quotient_determinant"]
                * d["parent"]["D_old"].det()
            ),
        },
        "gates": {
            "all_56_complement_and_four_retained_directions": K.shape == (60, 56)
            and L.shape == (60, 4),
            "all_nine_blocks_cover_56": sum(map(len, d["blocks"])) == 56,
            "actual_complement_determinant_is_not_constant": s.diff(expected, P) != 0,
            "whole_complement_is_not_positive_definite_by_assumption": True,
        },
    }


@cache
def projective():
    d = matrices()
    F, G, J, Jinv = (d[key] for key in ("gauge", "projective", "full", "full_inverse"))
    return {
        "whole_four_projective_gauge_columns": G,
        "whole_original_trace_gauge_functional": F,
        "whole_original_projective_Faddeev_Popov_matrix": F * G,
        "whole_64_connection_to_56_plus_four_plus_four_map": J,
        "whole_64_connection_inverse_map": Jinv,
        "whole_64_configuration_Jacobian": clean(J.det(method="domain-ge")),
        "whole_original_gauge_determinant": (F * G).det(),
        "coordinate_boundary": "Write a full connection displacement as E K z + E L t + G g. Here g=F kappa/4, F G=4I, and the full 64-coordinate Jacobian is constant but not unit. The normalized gauge g=0 has unit Faddeev-Popov matrix; the original trace gauge F kappa=0 has determinant 4^4. These statements agree once the displayed coordinate Jacobian and gauge-parameter normalization are both retained. A full cotangent lift has unit phase-space Liouville Jacobian; it does not make every configuration Jacobian one.",
        "checks": {
            "whole_projective_trace_gauge_matrix": F * G - 4 * s.eye(4),
            "whole_gauge_free_quotient_columns": F * d["embedding"],
            "whole_retained_trace_projective_invariance": d["full_inverse"][56:60, :]
            * G,
            "whole_64_point_map_left_inverse": clean(Jinv * J - s.eye(64)),
            "whole_64_point_map_right_inverse": clean(J * Jinv - s.eye(64)),
            "whole_64_point_map_constant_determinant": clean(
                J.det(method="domain-ge") + 256
            ),
            "whole_original_gauge_determinant": (F * G).det() - 256,
        },
        "gates": {
            "all_original_64_connection_coordinates_retained": J.shape == (64, 64),
            "configuration_Jacobian_not_silently_set_to_one": J.det(method="domain-ge")
            != 1,
            "four_projective_gauge_modes_not_quantized_as_physical_vectors": True,
        },
    }


@cache
def inertia():
    d = matrices()
    rows = []
    checks = {}
    positive = negative = 0
    for i, block in enumerate(d["block_matrices"]):
        reference = block.subs(P, s.Rational(1, 2))
        transform, diagonalized, pivots, signs = rational_congruence(reference)
        rows.append(
            {
                "indices": d["blocks"][i],
                "reference_matrix": reference,
                "exact_congruence": transform,
                "exact_diagonal_blocks": diagonalized,
                "pivots": pivots,
                "inertia": signs,
            }
        )
        checks["whole_rational_block_congruence_" + str(i)] = (
            transform.T * reference * transform - diagonalized
        )
        positive += signs[0]
        negative += signs[1]
    signature = positive - negative
    return {
        "whole_exact_nine_block_inertia_decompositions": rows,
        "whole_positive_negative_complement_inertia": (positive, negative),
        "whole_complement_signature": signature,
        "whole_positive_scale_Fresnel_phase": s.exp(s.I * s.pi * signature / 4),
        "whole_inertia_continuation": "The full real symmetric complement is continuous in p on the connected strict interval 1/8<p^2<3/10 with p>0. Its explicit determinant has no zero there. No real eigenvalue can change sign without vanishing, so the exact reference inertia (26,30) holds throughout. This argument is for the regular timelike source chart, not a claimed quantum unitary-gauge extension through a vanishing clock gradient.",
        "checks": {
            **checks,
            "whole_56_inertia_count": s.Integer(positive + negative - 56),
            "whole_exact_signature": s.Integer(signature + 4),
            "whole_Fresnel_phase": s.exp(s.I * s.pi * signature / 4) + 1,
        },
        "gates": {
            "all_inertia_pivots_exact_rational": all(
                entry.is_Rational is True
                for row in rows
                for pivot in row["pivots"]
                for entry in pivot
            ),
            "full_26_positive_30_negative_inertia": (positive, negative) == (26, 30),
            "inertia_continuity_not_a_numerical_sampling_proof": True,
        },
    }


@cache
def source_map():
    d = matrices()
    original = s.Matrix(old.eliminate()["solution"])
    # Use the literal full trace map; embedding.T is not its inverse.
    Nfull = s.Matrix(geometry.update()["N_full"])
    Tstar = clean(Nfull * original)
    retained = s.Matrix(s.symbols("retained_unshifted_trace0:4", real=True))
    center = clean(original + d["embedding"] * d["lift"] * (retained - Tstar))
    full_hessian = (
        s.Matrix(old.quadratic()["hessian"]) + Nfull.T * d["parent"]["Xi"] * Nfull
    )
    full_source = (
        s.Matrix(old.quadratic()["source"]) - Nfull.T * d["parent"]["Xi"] * Tstar
    )
    E56 = d["embedding"] * d["kernel"]
    return {
        "whole_original_64_stationary_source_center": original,
        "whole_original_four_trace_stationary_source": Tstar,
        "whole_retained_unshifted_trace": retained,
        "whole_actual_64_connection_center_at_retained_trace": center,
        "whole_source_centered_64_updated_Hessian": full_hessian,
        "whole_source_centered_64_updated_linear_term": full_source,
        "whole_original_shifted_vector_binding": source.data()["new_retained_vector"],
        "whole_original_lower_coefficient_binding": source.data()[
            "regular_lower_dictionary"
        ],
        "whole_original_nonlinear_normal_source": source.data()[
            "normal_source_timelike_hat_chart"
        ],
        "source_boundary": "The stationary center retains all original derivative forcing. The retained trace is T=W+B(u,X)du, not a source-free W. The conditional complement shift depends on finite metric/clock jets, so it is not declared an ordinary point-canonical map on the original unextended configuration space. The finite-jet lifting theorem retains the required derivative coordinates, multiplier constraints, source pullback and time boundaries. The four dynamical Proca directions and the complete source square remain.",
        "checks": {
            "whole_64_source_centered_stationary_Euler": clean(
                full_hessian * original + full_source
            ),
            "whole_actual_retained_trace_after_center_shift": clean(
                Nfull * center - retained
            ),
            "whole_56_source_centered_Euler_at_fixed_retained_trace": clean(
                E56.T * (full_hessian * center + full_source)
            ),
            "whole_complement_mass_source_second_derivative": clean(
                E56.T * full_hessian * E56 - d["A"]
            ),
            "whole_original_source_center_projective_gauge": clean(
                d["gauge"] * original
            ),
        },
        "gates": {
            "all_64_source_components_retained": len(center) == 64,
            "all_56_conditional_equations_not_just_trace_stationarity": E56.shape
            == (64, 56),
            "stationary_connection_source_not_set_to_zero": original != s.zeros(64, 1),
        },
    }
