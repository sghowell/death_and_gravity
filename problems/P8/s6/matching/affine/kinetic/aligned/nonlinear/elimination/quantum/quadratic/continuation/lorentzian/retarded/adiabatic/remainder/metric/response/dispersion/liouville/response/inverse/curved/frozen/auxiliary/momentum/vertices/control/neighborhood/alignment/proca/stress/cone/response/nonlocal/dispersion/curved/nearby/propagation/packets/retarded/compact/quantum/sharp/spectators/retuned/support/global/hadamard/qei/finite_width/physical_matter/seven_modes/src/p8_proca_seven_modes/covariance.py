"""Exact Minkowski physical covariances and the Proca constraint normalization."""

from functools import cache

import sympy as sp

from . import model

omega = sp.Symbol("positive_massive_frequency", positive=True)


def on_shell(value):
    def scalar(v):
        numerator, denominator = sp.cancel(v).as_numer_denom()
        numerator = sp.rem(numerator, omega**2 - model.m**2 - model.q, omega)
        return sp.factor(numerator / denominator)

    return (
        value.applyfunc(scalar) if isinstance(value, sp.MatrixBase) else scalar(value)
    )


@cache
def data():
    m, k, q = model.m, model.kvec, model.q
    eye = sp.eye(3)
    Q = (m * m * eye + k * k.T) / omega
    P = omega * (eye - k * k.T / omega**2) / (m * m)
    C = sp.zeros(6)
    C[:3, :3] = Q / 2
    C[3:, 3:] = P / 2
    C[:3, 3:] = sp.I * eye / 2
    C[3:, :3] = -sp.I * eye / 2
    # A positive square root with no angular polarization frame or q divisor.
    root = (m * eye + k * k.T / (omega + m)) / sp.sqrt(omega)
    inverse = sp.sqrt(omega) * (eye - k * k.T / (omega * (omega + m))) / m
    factor = root.col_join(-sp.I * inverse) / sp.sqrt(2)
    M = model.data()["Proca_cartesian_generator"].subs(model.a, 1)
    L = (
        model.data()["Proca_temporal_reconstruction"]
        .subs(model.a, 1)
        .col_join(sp.eye(3).row_join(sp.zeros(3)))
    )
    physical = on_shell(L * C * L.conjugate().T)
    eta = sp.diag(-1, 1, 1, 1)
    xi = sp.Matrix([-omega, *list(k)])
    expected = m * m * (eta + xi * xi.T / (m * m)) / (2 * omega)
    eps = (-k.T / m).col_join(eye + k * k.T / (m * (omega + m)))
    polcov = m * m * eps * eps.T / (2 * omega)
    # Two tensor polarizations are T/sqrt(2) ordinary minimal scalars.
    kt = sp.Symbol("positive_tensor_frequency", positive=True)
    CT = sp.Matrix([[1 / kt, sp.I / 2], [-sp.I / 2, kt / 4]])
    GT = sp.Matrix([1 / sp.sqrt(kt), -sp.I * sp.sqrt(kt) / 2])
    MT = model.data()["tensor_polarization_generator"].subs(
        {model.a: 1, model.q: kt * kt}, simultaneous=True
    )
    # A concrete null covector exhibits the nonelliptic but nonzero symbol.
    null = sp.Matrix([1, 0, 0, 1])
    principal = null * (null.T * eta) / (m * m)
    return {
        "mass_shell_relation": omega**2 - m * m - q,
        "Minkowski_Proca_density_covariance_without_hbar_over_kappa": C,
        "positive_Proca_Gram_factor_without_sqrt_hbar_over_kappa": factor,
        "physical_four_component_Proca_covariance": physical,
        "physical_Proca_polarization_factor": eps,
        "Minkowski_tensor_mode_covariance_without_hbar_over_kappa": CT,
        "nonelliptic_nonzero_null_Proca_projection_symbol": principal,
        "checks": {
            "Proca_initial_covariance_is_Hermitian": on_shell(C - C.conjugate().T),
            "Proca_initial_covariance_has_actual_density_CCR": on_shell(
                C - C.T - sp.I * model.data()["Proca_symplectic_form"]
            ),
            "Proca_covariance_is_positive_Gram_product": on_shell(
                C - factor * factor.conjugate().T
            ),
            "Proca_positive_square_root_is_exact": on_shell(root * root - Q),
            "Proca_positive_square_root_inverse_is_exact": on_shell(
                root * inverse - eye
            ),
            "Proca_covariance_is_stationary_in_actual_flat_Hamiltonian": on_shell(
                M * C + C * M.T
            ),
            "constrained_four_component_covariance_has_correct_metric_mass_polarization": on_shell(
                physical - expected
            ),
            "constrained_four_component_covariance_has_positive_three_column_factor": on_shell(
                physical - polcov
            ),
            "three_massive_polarizations_are_mass_shell_transverse": on_shell(
                xi.T * eta * eps
            ),
            "three_massive_polarizations_have_positive_unit_norm": on_shell(
                eps.T * eta * eps - sp.eye(3)
            ),
            "tensor_mode_covariance_has_actual_CCR": sp.simplify(
                CT - CT.T - sp.I * model.J2
            ),
            "tensor_mode_covariance_is_positive_Gram_product": sp.simplify(
                CT - GT * GT.conjugate().T
            ),
            "tensor_flat_covariance_is_stationary": sp.simplify(MT * CT + CT * MT.T),
            "Proca_projection_principal_symbol_is_not_elliptic": sp.factor(
                principal.det()
            ),
            "Proca_projection_symbol_is_nonzero_rank_one_on_a_null_covector": sp.Integer(
                principal.rank() - 1
            ),
            "Proca_null_projection_symbol_is_nilpotent": sp.simplify(
                principal * principal
            ),
        },
    }
