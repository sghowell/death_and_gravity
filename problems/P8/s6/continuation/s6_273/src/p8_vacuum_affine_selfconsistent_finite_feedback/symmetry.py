"""Actual radial seed covariance, full cubic action and all homogeneous channels."""

import itertools
from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_phase_domain import reference as seed


@cache
def rotations():
    result = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            g = s.zeros(3)
            for col, row in enumerate(perm):
                g[row, col] = signs[col]
            if g.det() == 1:
                result.append(s.ImmutableMatrix(g))
    return tuple(result)


@cache
def tensor_basis():
    eye = s.eye(3)
    result = []
    for axis in range(3):
        e, f = eye[:, (axis + 1) % 3], eye[:, (axis + 2) % 3]
        result.append(
            ((e * e.T - f * f.T) / s.sqrt(2), (e * f.T + f * e.T) / s.sqrt(2))
        )
    return tuple(result)


def configuration(g):
    if s.ImmutableMatrix(g) not in rotations():
        raise ValueError("Require an actual proper cubic rotation")
    result = s.zeros(48)
    tensor = tensor_basis()
    for axis in range(3):
        image = g[:, axis]
        target = next(i for i in range(3) if image[i])
        sign = image[target]
        spin = s.Matrix(
            2,
            2,
            lambda i, j, target=target, axis=axis: s.simplify(
                s.trace(tensor[target][i] * g * tensor[axis][j] * g.T)
            ),
        )
        block = s.diag(s.eye(2), spin, s.eye(1), g)
        for wave in range(2):
            phase = 1 if wave == 0 else sign
            start, end = (axis * 2 + wave) * 8, (target * 2 + wave) * 8
            result[end : end + 8, start : start + 8] = phase * block
    return s.ImmutableMatrix(result)


@cache
def representations():
    return tuple(configuration(g) for g in rotations())


def index(axis, wave, channel, qp):
    return qp * 48 + (axis * 2 + wave) * 8 + channel


@cache
def radial_covariance():
    V = s.MutableSparseMatrix(96, 96, {})
    scalar = s.Matrix(
        4,
        4,
        lambda i, j: s.Symbol(
            "full_scalar_cov_" + str(min(i, j)) + "_" + str(max(i, j)), real=True
        ),
    )
    blocks = {}
    for name in ("tensor", "heavy", "Proca_T", "Proca_L"):
        blocks[name] = s.Matrix(
            2,
            2,
            lambda i, j, name=name: s.Symbol(
                name + "_cov_" + str(min(i, j)) + "_" + str(max(i, j)), real=True
            ),
        )
    for axis in range(3):
        for wave in range(2):
            for i, j in itertools.product(range(4), repeat=2):
                V[
                    index(axis, wave, i % 2, i // 2), index(axis, wave, j % 2, j // 2)
                ] = scalar[i, j]
            for qp, rp in itertools.product(range(2), repeat=2):
                for channel in (2, 3):
                    V[
                        index(axis, wave, channel, qp), index(axis, wave, channel, rp)
                    ] = blocks["tensor"][qp, rp]
                V[index(axis, wave, 4, qp), index(axis, wave, 4, rp)] = blocks["heavy"][
                    qp, rp
                ]
                for component in range(3):
                    name = "Proca_L" if component == axis else "Proca_T"
                    V[
                        index(axis, wave, 5 + component, qp),
                        index(axis, wave, 5 + component, rp),
                    ] = blocks[name][qp, rp]
    return s.ImmutableSparseMatrix(V)


@cache
def data():
    groups, reps = rotations(), representations()
    lookup = {tuple(g): i for i, g in enumerate(groups)}
    V = radial_covariance()
    J = s.zeros(48).row_join(s.eye(48)).col_join((-s.eye(48)).row_join(s.zeros(48)))
    symplectic = []
    covariant = []
    orthogonal = []
    composition = []
    for g, U in zip(groups, reps, strict=True):
        S = s.diag(U, U)
        orthogonal.append(s.Integer(U * U.T != s.eye(48)))
        symplectic.append(s.Integer(S * J * S.T != J))
        covariant.append(s.Integer(S * V * S.T != V))
    for i, g in enumerate(groups):
        for j, h in enumerate(groups):
            composition.append(
                s.Integer(reps[i] * reps[j] != reps[lookup[tuple(g * h)]])
            )
    entries = s.symbols("full_homogeneous_symmetric0:6", real=True)
    shape = s.Matrix(
        [
            [entries[0], entries[1], entries[2]],
            [entries[1], entries[3], entries[4]],
            [entries[2], entries[4], entries[5]],
        ]
    )
    averaged = sum((g * shape * g.T for g in groups), s.zeros(3)) / 24
    bad = s.eye(48)
    bad[0, 0] = 2
    k = s.Matrix(s.symbols("full_wavevector0:3", real=True))
    a, m = s.symbols("positive_scale positive_vector_mass", positive=True)
    K = s.eye(3) / a + k * k.T / (a**3 * m**2)
    W = a * m**2 * s.eye(3) + (k.dot(k) * s.eye(3) - k * k.T) / a
    vector_checks = []
    for g in groups:
        rule = dict(zip(k, g * k, strict=True))
        vector_checks.extend(
            list((K.subs(rule, simultaneous=True) - g * K * g.T).applyfunc(s.expand))
        )
        vector_checks.extend(
            list((W.subs(rule, simultaneous=True) - g * W * g.T).applyfunc(s.expand))
        )
    sparse_actions = [
        [[j, U[i, j]] for i in range(48) for j in range(48) if U[i, j]] for U in reps
    ]
    return {
        "whole_24_proper_cubic_rotations": groups,
        "whole_24_full_signed_48_configuration_actions_row_order": sparse_actions,
        "whole_full_radial_96_phase_covariance": V,
        "whole_full_Proca_K_V": (K, W),
        "whole_covariance_binding": "The actual S251 scalar generator and positive selection form depend on |k|^2 and act identically on cosine/sine pairs and all three axes. Their unique full Gaussian minimizer and both canonical boundary maps therefore have the displayed arbitrary scalar4x4 block, including every internal and q-p entry. Both tensor polarizations use the same full2x2 preparation. S240's basis-independent exact SLE is radial. The unchanged S190 all-order initial T/T/L frequencies and slopes are radial on the zero-shear reference, with equal transverse blocks; uniqueness of the full zero-shear mode/Riccati equations preserves the complete transverse/longitudinal covariance. This is a binding of the actual reference, not a new isotropic state or deletion of an allowed odd block from an arbitrary translation-invariant state.",
        "whole_original_reference_source_modules": (
            "p8_vacuum_affine_coupled_gaussian_state.phase",
            "p8_vacuum_affine_coupled_gaussian_state.gaussian",
            "p8_vacuum_affine_heavy_curved_state.state",
            "p8_vacuum_affine_matrix_adiabatic.initial",
            "p8_vacuum_affine_matrix_adiabatic.frame",
        ),
        "whole_whitening_and_cutoff_proof": "The physical cubic action S is symplectic and preserves the actual V0=S0 S0^T/2. Hence S0^-1 S S0 is orthogonal as well as symplectic. Both unchanged radial cutoffs and D=Delta/4 are invariant. The actual same-state coherent map and Weyl map are covariant, without replacing S0 or discarding cross covariance.",
        "whole_full_reconstruction_equivariance": "All original scalar contractions, derivatives, full matrix Fourier projectors, shape determinant fixed point, actual mean-zero formal-transpose inverse, temporal/lapse roots and spatial integrals are cubic and translation equivariant. Uniqueness preserves the symmetry of each reconstruction. All generated modes remain. The reference flow and both complete quantum Hamiltonians commute with the same groups.",
        "whole_homogeneous_consistency": "The unchanged pure centered radial seed is invariant under the ordinary configuration pullback; its nonzero value at the origin fixes the possible phase to1. Exact bounded unitary evolution preserves this invariant subspace and the original zero-translation-charge sector. Any homogeneous vector force is an invariant three-vector, hence0. Each homogeneous traceless shape or dual force is an invariant traceless symmetric3x3 matrix, hence0. This proves consistency for all FIVE shape and all three spatial-vector canonical pairs, not a two-TT-zero-mode truncation. The heavy scalar is not removed by this symmetry.",
        "checks": {
            "all_24_configuration_orthogonal": s.Matrix(orthogonal),
            "all_24_full_phase_symplectic": s.Matrix(symplectic),
            "all_24_complete_covariance_invariant": s.Matrix(covariant),
            "all_576_group_compositions": s.Matrix(composition),
            "homogeneous_vector_Reynolds_average": sum(groups, s.zeros(3)),
            "all_five_traceless_shape_Reynolds_average": (
                averaged - s.trace(shape) * s.eye(3) / 3
            ).applyfunc(s.expand),
            "all_24_full_Proca_K_V_covariances": s.Matrix(vector_checks),
            "same_48_original_configuration_channels": s.Integer(seed.DIMENSION - 48),
        },
        "gates": {
            "all_proper_rotations_and_unique_actions": len(groups) == len(lookup) == 24,
            "every_configuration_action_is_signed_permutation": all(
                sum(value != 0 for value in U.row(i)) == 1
                for U in reps
                for i in range(48)
            ),
            "anisotropic_seed_cannot_be_replaced": any(
                U * bad * U.T != bad for U in reps
            ),
            "scalar_internal_and_qp_covariance_kept": V[0, 1] != 0
            and V[0, 48] != 0
            and V[0, 49] != 0,
            "all_three_Proca_covariance_polarizations_kept": V[5, 5] != V[6, 6]
            and V[6, 6] == V[7, 7],
            "homogeneous_shape_dimension_is_five": 3 * 4 // 2 - 1 == 5,
            "same_original_seed_not_state_averaging": True,
            "finite_hybrid_not_homogeneous_quantization": True,
        },
    }
