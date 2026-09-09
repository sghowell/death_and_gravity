"""Tensor action variation before restriction to the physical clock."""

from functools import cache

import sympy as sp
from p8_affine_nonlinear import adm
from p8_proca_seven_modes import model

u = adm.u
N = sp.Symbol("positive_unchanged_physical_lapse", positive=True)
a = sp.Symbol("positive_physical_scale", positive=True)
A = sp.Symbol("positive_hat_scale", positive=True)
h = sp.Symbol("positive_clock_h", positive=True)
q = sp.Symbol("positive_comoving_momentum_squared", positive=True)
T, P, V = sp.symbols(
    "tensor_coordinate tensor_density_momentum tensor_velocity", real=True
)


@cache
def data():
    ratio = (h - 1 + N**-2) / h
    e = ratio ** (-sp.Rational(1, 4))
    L = ratio * a**3 * V * V / (4 * N) - N * ratio * a * q * T * T / 4
    velocity = 2 * N * P / (ratio * a**3)
    Hamiltonian = N * P * P / (ratio * a**3) + N * ratio * a * q * T * T / 4
    hat = N * e * P * P / A**3 + N * A * q * T * T / (4 * e**3)
    r = P * P / a**6 + q * T * T / (4 * a * a)
    s = P * P / a**6 - q * T * T / (12 * a * a)
    rho = sp.factor(sp.diff(Hamiltonian, N).subs(N, 1) / a**3)
    pressure = sp.factor(-sp.diff(Hamiltonian, a).subs(N, 1) / (3 * a * a))
    force = sp.factor(sp.diff(hat, N).subs(N, 1) / A**3).subs(A, a)
    actual_GT = -2 * adm.coefficients()["B"]
    matrix = sp.Matrix([[0, 2 / a**3], [-a * q / 2, 0]])
    x, y = P * P / a**6, q * T * T / (4 * a * a)
    return {
        "actual_physical_tensor_Lagrangian": L,
        "actual_physical_tensor_Hamiltonian": Hamiltonian,
        "actual_spatial_hat_tensor_Hamiltonian": hat,
        "positive_on_clock_tensor_energy": r,
        "on_clock_spatial_pressure": s,
        "fixed_physical_metric_lapse_source": rho,
        "actual_hat_lapse_source": force,
        "clock_coupled_kinetic_scalar": r - 3 * s,
        "checks": {
            "literal_actual_DHOST_tensor_coefficient_before_clock_restriction": sp.factor(
                actual_GT - (1 + (adm.X - 1) / (1 + u * u) ** 3)
            ),
            "actual_tensor_Legendre_transform_before_clock_restriction": sp.factor(
                P * velocity - L.subs(V, velocity) - Hamiltonian
            ),
            "actual_tensor_spatial_frame_Hamiltonian": sp.powsimp(
                Hamiltonian.subs(a, e * A) - hat, force=True
            ),
            "actual_tensor_lapse_source_is_not_minimal_matter_density": sp.factor(
                rho - (1 - 1 / h) * r - 3 * s / h
            ),
            "actual_tensor_spatial_pressure": sp.factor(pressure - s),
            "actual_tensor_hat_lapse_source_keeps_clock_coupling": sp.factor(
                force - (1 - 1 / h) * r - 3 * s / (2 * h)
            ),
            "actual_tensor_frame_lapse_chain_retains_pressure": sp.factor(
                force - rho + 3 * s / (2 * h)
            ),
            "actual_tensor_lapse_force_upper_uses_h_at_least_one": sp.factor(
                sp.Rational(3, 2) * r
                - force
                - (1 - 1 / h) * x / 2
                - (1 + 3 / h) * y / 2
            ),
            "actual_tensor_lapse_force_lower_uses_h_at_least_one": sp.factor(
                sp.Rational(3, 2) * r
                + force
                - (5 + 1 / h) * x / 2
                - (5 - 3 / h) * y / 2
            ),
            "on_clock_tensor_generator_matches_existing_state_normalization": sp.Matrix(
                [
                    [0, sp.diff(Hamiltonian, P, 2).subs(N, 1)],
                    [-sp.diff(Hamiltonian, T, 2).subs(N, 1), 0],
                ]
            )
            - matrix,
            "on_clock_tensor_generator_is_literal_S6_103_generator": (
                matrix
                - model.data()["tensor_polarization_generator"].subs(
                    {model.a: a, model.q: q}, simultaneous=True
                )
            ).applyfunc(sp.factor),
            "on_clock_tensor_generator_has_actual_CCR": matrix * model.J2
            + model.J2 * matrix.T,
        },
    }


@cache
def clock():
    r, s, rd, sd, H, hp = sp.symbols(
        "proxy_density proxy_pressure proxy_density_derivative proxy_pressure_derivative H h_derivative",
        real=True,
    )
    K = r - 3 * s
    Kd = rd - 3 * sd
    source = Kd / h + (3 * H / h - hp / h**2) * K
    physical = (1 - 1 / h) * r + 3 * s / h
    physical_d = (1 - 1 / h) * rd + 3 * sd / h + (r - 3 * s) * hp / h**2
    residual = (physical_d + 3 * H * (physical + s) + source).subs(rd, -3 * H * (r + s))
    return {
        "proxy_density": r,
        "proxy_pressure": s,
        "clock_source": source,
        "physical_fixed_metric_lapse_source": physical,
        "clock_source_is_not_generically_zero": True,
        "checks": {
            "actual_homogeneous_tensor_source_Ward_identity_keeps_clock_equation": sp.factor(
                residual
            ),
            "tensor_source_difference_from_minimal_Proca_is_nonzero_clock_kinetic_term": sp.factor(
                (1 - 1 / h) * r
                + 3 * s / (2 * h)
                - (r - 3 * s / (2 * h))
                + (r - 3 * s) / h
            ),
        },
    }
