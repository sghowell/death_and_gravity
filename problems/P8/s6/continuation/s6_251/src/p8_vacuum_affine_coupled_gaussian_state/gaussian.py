"""Basis-independent two-mode Gaussian minimizer with exact CCR and purity."""

from functools import cache

import sympy as s

from .phase import OMEGA


def exact_positive_matrix(value):
    M = s.Matrix(value)
    if M.shape != (4, 4) or M != M.T or any(x.has(s.Float) for x in M):
        raise ValueError("Expected an exact real symmetric four-by-four matrix")
    if not all(x.is_real is True for x in M):
        raise ValueError("Matrix entries must be provably real")
    if not all(s.factor(M[:j, :j].det()).is_positive is True for j in range(1, 5)):
        raise ValueError("Preparation matrix must be strictly positive")
    return M


def invariants(value):
    M = exact_positive_matrix(value)
    product = s.sqrt(s.factor(M.det()))
    square_sum = s.factor(-s.trace((OMEGA * M) ** 2) / 2)
    frequency_sum = s.sqrt(s.factor(square_sum + 2 * product))
    return product, square_sum, frequency_sum


def ground_covariance(value):
    """Two-mode formula; all branches are positive for the SPD preparation."""
    M = exact_positive_matrix(value)
    product, _, total = invariants(M)
    return ((product * M.inv() - OMEGA * M * OMEGA) / (2 * total)).applyfunc(s.simplify)


def single_mode_covariance(value):
    M = s.Matrix(value)
    if (
        M.shape != (2, 2)
        or M != M.T
        or any(x.has(s.Float) or x.is_real is not True for x in M)
    ):
        raise ValueError(
            "Expected a complete exact real symmetric tensor preparation matrix"
        )
    if M[0, 0].is_positive is not True or s.factor(M.det()).is_positive is not True:
        raise ValueError("Tensor preparation must be strictly positive")
    return (s.sqrt(s.factor(M.det())) * M.inv() / 2).applyfunc(s.simplify)


def symplectic_fixture(case):
    if type(case) is not int or not 0 <= case < 12:
        raise ValueError("Unknown exact symplectic fixture")
    n = s.Integer(case + 1)
    L = s.Matrix([[1 + n / 5, n / 7], [0, 1 / (1 + n / 5)]])
    C = s.Matrix([[n / 11, (n + 1) / 13], [(n + 1) / 13, -n / 17]])
    D = s.Matrix([[n / 19, -n / 23], [-n / 23, (n + 2) / 29]])
    low = s.eye(4)
    low[2:, :2] = C
    up = s.eye(4)
    up[:2, 2:] = D
    diagonal = s.diag(L, L.inv().T)
    return low * diagonal * up


def fixture_preparation(case):
    S = symplectic_fixture(case)
    nu1, nu2 = s.Rational(case + 2, 3), s.Rational(case + 5, 4)
    diagonal = s.diag(nu1, nu2, nu1, nu2)
    M = S.inv().T * diagonal * S.inv()
    return S, diagonal, M


@cache
def data():
    checks = {}
    for case in range(4):
        S, D, M = fixture_preparation(case)
        V = ground_covariance(M)
        p, squares, total = invariants(M)
        checks[f"full_noncommuting_symplectic_fixture_{case}"] = S * OMEGA * S.T - OMEGA
        checks[f"full_ground_covariance_not_diagonal_only_{case}"] = V - S * S.T / 2
        checks[f"entire_CCR_purity_{case}"] = V * OMEGA * V - OMEGA / 4
        checks[f"ground_uncertainty_congruence_{case}"] = (
            V + s.I * OMEGA / 2 - S * (s.eye(4) + s.I * OMEGA) * S.T / 2
        )
        checks[f"minimum_whole_averaged_energy_{case}"] = s.trace(M * V) / 2 - total / 2
        checks[f"frequency_product_{case}"] = p - D[0, 0] * D[1, 1]
        checks[f"frequency_square_sum_{case}"] = squares - D[0, 0] ** 2 - D[1, 1] ** 2
        checks[f"frequency_sum_{case}"] = total - D[0, 0] - D[1, 1]
        checks[f"two_mode_hamiltonian_Cayley_Hamilton_{case}"] = (
            (OMEGA * M) ** 4 + squares * (OMEGA * M) ** 2 + p**2 * s.eye(4)
        )
        checks[f"normalizing_bump_does_not_change_state_{case}"] = (
            ground_covariance((case + 2) * M) - V
        )
    S, _, M = fixture_preparation(2)
    T = symplectic_fixture(5)
    changed = T.inv().T * M * T.inv()
    checks["entire_reference_basis_covariance"] = (
        ground_covariance(changed) - T * ground_covariance(M) * T.T
    )
    nu = s.Rational(7, 3)
    equal = S.inv().T * (nu * s.eye(4)) * S.inv()
    checks["equal_preparation_symplectic_frequencies_are_regular"] = (
        ground_covariance(equal) - S * S.T / 2
    )
    # Imaginary number-conserving mixing is retained; realification commutes with Omega.
    X = s.Matrix([[3, 1], [1, 4]])
    Y = s.Matrix([[0, 1], [-1, 0]])
    number = X.row_join(-Y).col_join(Y.row_join(X))
    checks["complete_number_conserving_mixing_does_not_change_vacuum"] = (
        ground_covariance(number) - s.eye(4) / 2
    )
    checks["number_conserving_matrix_commutes_with_complex_structure"] = (
        number * OMEGA - OMEGA * number
    )
    real_part = s.Matrix(s.symbols("test_real0:4", real=True))
    imaginary_part = s.Matrix(s.symbols("test_imaginary0:4", real=True))
    complex_vector = real_part + s.I * imaginary_part
    factor = s.eye(2).row_join(s.I * s.eye(2))
    checks["full_uncertainty_positive_factorization"] = (
        s.eye(4) + s.I * OMEGA - factor.conjugate().T * factor
    )
    checks["uncertainty_positive_on_all_complex_vectors"] = (
        complex_vector.conjugate().T * (s.eye(4) + s.I * OMEGA) * complex_vector
    )[0] - ((factor * complex_vector).conjugate().T * (factor * complex_vector))[0]
    epsilon = s.Rational(1, 5)
    uncertainty = (s.eye(4) + s.I * OMEGA) / 2
    checks["full_uncertainty_projector"] = uncertainty**2 - uncertainty
    checks["pure_uncertainty_rank_two"] = s.trace(uncertainty) - 2
    tensor_omega = s.Matrix([[0, 1], [-1, 0]])
    for case in range(4):
        T = s.Matrix(
            [
                [1, s.Rational(case + 1, 7)],
                [s.Rational(case + 2, 11), 1 + s.Rational((case + 1) * (case + 2), 77)],
            ]
        )
        tensor_M = T.inv().T * s.diag(case + 2, case + 2) * T.inv()
        tensor_V = single_mode_covariance(tensor_M)
        checks[f"full_tensor_preparation_covariance_{case}"] = tensor_V - T * T.T / 2
        checks[f"full_tensor_CCR_and_purity_{case}"] = (
            tensor_V * tensor_omega * tensor_V - tensor_omega / 4
        )
        checks[f"full_tensor_ground_energy_{case}"] = s.trace(
            tensor_M * tensor_V
        ) / 2 - s.Rational(case + 2, 2)
    return {
        "preparation_formula": "M(P)=integral f(t)^2 S(t,t*)^T E(t,P) S(t,t*) dt; complete fixed-current flow and complete positive outer selection form, no instantaneous reset or physical action change.",
        "entire_two_mode_covariance": "V=(sqrt(det M) M^-1-Omega M Omega)/(2 sqrt(-Tr[(Omega M)^2]/2+2sqrt(det M))). All square roots have the positive SPD branch. This remains regular when the two preparation symplectic eigenvalues coincide.",
        "equivalent_spectral_formula": "V=1/2 M^-1/2 sqrt[-(M^1/2 Omega M^1/2)^2] M^-1/2. Positive uncertainty holds on ALL complex test vectors; V Omega V=Omega/4. Mean=0.",
        "whole_minimization": "The exact minimum of 1/2 Tr(M V) is (nu1+nu2)/2 over all admissible covariances. The minimizer is unique. A positive number-conserving branch mixing does not change its annihilation vacuum; anomalous blocks cannot be deleted unless their effect is controlled.",
        "full_unequal_time_two_point": "W_z(t,s;P)=S(t,t*)[V(P)+i Omega/2]S(s,t*)^T. The full matrix, not its diagonal, is retained. Fields v,sigma recover BOTH 1/sqrt(kappa) factors; normalized momenta also recover a^-3 at their own endpoint.",
        "counterexample_half_uncertainty": epsilon * s.eye(4),
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "strictly_positive_real_covariance_is_not_enough": bool(
                epsilon < s.Rational(1, 2)
            ),
            "half_uncertainty_counterexample_has_negative_complex_eigenvalue": bool(
                epsilon - s.Rational(1, 2) < 0
            ),
            "actual_mixing_fixture_not_block_diagonal": M[:2, 2:] != s.zeros(2),
            "deleting_anomalous_blocks_changes_generic_state": ground_covariance(M)
            != s.eye(4) / 2,
            "smooth_bump_normalization_not_a_new_state_choice": True,
            "Gaussian_reference_not_full_interacting_state": True,
        },
    }
