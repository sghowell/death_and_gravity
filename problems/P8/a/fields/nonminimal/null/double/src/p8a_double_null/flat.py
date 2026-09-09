"""Exact physical stress identity, timelike-plane spectrum and Fourier factors."""

from functools import cache

import sympy as sp


@cache
def data():
    f, df, phi, dphi, ddphi, xi = sp.symbols("f Df phi Dphi DDphi xi", real=True)
    physical = f * f * ((1 - 2 * xi) * dphi * dphi - 2 * xi * phi * ddphi)
    positive = (
        (1 - 2 * xi) * f * f * dphi * dphi
        + 2 * xi * (df * phi + f * dphi) ** 2
        - 2 * xi * df * df * phi * phi
    )
    boundary = (
        -2 * xi * (2 * f * df * phi * dphi + f * f * dphi * dphi + f * f * phi * ddphi)
    )
    T = sp.Symbol("Fourier_time", positive=True)
    Z, kz = sp.symbols("Fourier_space kz", real=True)
    omega = sp.Symbol("omega", nonnegative=True)
    q = omega - kz
    inner = sp.integrate(
        sp.integrate(
            (1 - 2 * xi) * q * q + 2 * xi * (T + Z - q) ** 2, (kz, -omega, omega)
        ),
        (omega, 0, T),
    )
    polynomial = (
        sp.Rational(2, 3) * (1 - xi) * T**4
        + sp.Rational(4, 3) * xi * T**3 * Z
        + 2 * xi * T * T * Z * Z
    )
    form = sp.Matrix(
        [
            [sp.Rational(2, 3) * (1 - xi), sp.Rational(2, 3) * xi],
            [sp.Rational(2, 3) * xi, 2 * xi],
        ]
    )
    # d²p=2pi p dp; p dp=omega d omega in the actual massless measure.
    spectral_factor = 2 * sp.pi / (2 * (2 * sp.pi) ** 3)
    fourier_half_factor = 2 / (2 * sp.pi) ** 2
    half_parseval_factor = (2 * sp.pi) ** 2 / 2
    L, M = sp.symbols("L M", real=True)
    derivative_form = (sp.Matrix([L, M]).T * form * sp.Matrix([L, M]))[0]
    return {
        "xi": xi,
        "reference_inner_polynomial": sp.expand(inner),
        "proper_derivative_quadratic_form": form,
        "timelike_plane_vacuum_spectral_factor_without_hbar": spectral_factor,
        "quantum_derivative_form_factor_without_hbar": spectral_factor
        * fourier_half_factor
        * half_parseval_factor,
        "checks": {
            "full_physical_nonminimal_null_sampler_identity": sp.expand(
                physical - positive - boundary
            ),
            "complete_timelike_plane_mass_shell_integral": sp.expand(
                inner - polynomial
            ),
            "plane_vacuum_spectral_normalization": sp.factor(
                spectral_factor - 1 / (8 * sp.pi**2)
            ),
            "two_dimensional_half_Fourier_and_Parseval_normalization": sp.factor(
                spectral_factor * fourier_half_factor * half_parseval_factor
                - 1 / (8 * sp.pi**2)
            ),
            "minimal_null_contraction_known_answer": sp.factor(
                form[0, 0].subs(xi, 0) / (8 * sp.pi**2) - 1 / (12 * sp.pi**2)
            ),
            "nonminimal_derivative_form_determinant": sp.factor(
                form.det() - sp.Rational(4, 9) * xi * (3 - 4 * xi)
            ),
            "full_derivative_form_positive_square_decomposition": sp.expand(
                derivative_form
                - 2 * xi * (M + L / 3) ** 2
                - sp.Rational(2, 9) * (3 - 4 * xi) * L * L
            ),
            "conformal_derivative_matrix": (
                form.subs(xi, sp.Rational(1, 6))
                - sp.Matrix(
                    [
                        [sp.Rational(5, 9), sp.Rational(1, 9)],
                        [sp.Rational(1, 9), sp.Rational(1, 3)],
                    ]
                )
            ),
            "reference_integrand_retains_full_null_derivative_of_sampler": sp.expand(
                (1 - 2 * xi) * q * q
                + 2 * xi * (T + Z - q) ** 2
                - q * q
                - 2 * xi * ((T + Z) ** 2 - 2 * (T + Z) * q)
            ),
        },
    }
