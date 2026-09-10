"""Literal full Hessian and the same fixed quartic Taylor subtraction."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_polynomial_vacuum import potential as previous


@cache
def data():
    M, L, g, G = sp.symbols(
        "positive_heavy_mass_squared positive_quartic positive_cubic_squared positive_cubic",
        positive=True,
    )
    phi, H = sp.symbols("constant_light constant_heavy", real=True)
    y = sp.Symbol("nonnegative_radial", nonnegative=True)
    V = phi**2 / 2 + M * H**2 / 2 + G * H * phi**2 / 2 + L * phi**4 / 24
    Hstar = -G * phi**2 / (2 * M)
    matrix = sp.hessian(V, (phi, H)) + sp.diag(y, y)
    Schur = (
        (matrix[0, 0] - matrix[0, 1] * matrix[1, 0] / matrix[1, 1])
        .subs(H, Hstar)
        .subs(G**2, g)
    )
    A = (L - g / M) / 2
    F = A - g / (y + M)
    C0 = -L + g / M
    V0 = C0 + 2 * g / (y + M)
    raw = sp.log(1 + F * phi**2 / (y + 1)) / 2
    fourth = sp.diff(raw, phi, 4).subs(phi, 0)
    old, par = previous.data(), model.data()
    old_y = sp.Symbol("Euclidean_momentum_squared", nonnegative=True)
    actual_F = F.subs(
        {
            M: par["heavy_mass_squared"],
            L: par["polynomial_quartic"],
            g: par["cubic_coupling_squared"],
            y: old_y,
        }
    )
    density = 6 * y * (F**2 - A**2) / ((y + 1) ** 2 * 16 * sp.pi**2)
    return {
        "M": M,
        "L": L,
        "g": g,
        "y": y,
        "literal_full_two_field_Hessian": matrix,
        "classical_stationary_heavy_field": Hstar,
        "stationary_light_Schur_kernel": Schur,
        "constant_background_kernel_F": F,
        "high_momentum_kernel_A": A,
        "complete_zero_external_vertex": V0,
        "fixed_finite_quartic_contact_radial_density": density,
        "fixed_contact_sign": "Strictly negative from 0<F<A at all finite radii",
        "same_potential_condition": old["one_loop_integral_convention"],
        "scope": "The constant-field quartic Taylor condition was already fixed in S6.110 and retained in S6.113. Only the finite contact remaining after the latter's entire I0 subtraction is evaluated here. The subsequent light on-shell quadratic conditions do not introduce a new quartic Taylor choice.",
        "checks": {
            "literal_stationary_full_Hessian_Schur": sp.factor(
                Schur - y - 1 - F * phi**2
            ),
            "same_frozen_actual_potential_kernel": sp.factor(
                actual_F - old["vertex_kernel_F"]
            ),
            "constant_external_full_vertex_equals_minus_twice_F": sp.factor(V0 + 2 * F),
            "full_log_fourth_derivative": sp.factor(fourth + 6 * F**2 / (y + 1) ** 2),
            "counterterm_fourth_derivative_cancels_potential_Taylor_term": sp.factor(
                fourth + sp.diff(6 * F**2 * phi**4 / (24 * (y + 1) ** 2), phi, 4)
            ),
            "remaining_finite_quartic_integrand_decomposition": sp.factor(
                F**2 - A**2 - C0 * g / (y + M) - g**2 / (y + M) ** 2
            ),
            "strict_negative_contact_integrand_factorization": sp.factor(
                F**2 - A**2 + g * (2 * A - g / (y + M)) / (y + M)
            ),
            "kernel_positive_lower_margin": sp.factor(
                F - (L - 3 * g / M) / 2 - g * y / (M * (y + M))
            ),
            "on_shell_quadratic_adjustment_has_no_quartic_Taylor_term": sp.diff(
                sp.Symbol("fixed_OS_quadratic") * phi**2, phi, 4
            ),
            "same_high_momentum_potential_quartic_subtraction": sp.factor(
                6 * A**2 - sp.Rational(3, 2) * C0**2
            ),
        },
    }
