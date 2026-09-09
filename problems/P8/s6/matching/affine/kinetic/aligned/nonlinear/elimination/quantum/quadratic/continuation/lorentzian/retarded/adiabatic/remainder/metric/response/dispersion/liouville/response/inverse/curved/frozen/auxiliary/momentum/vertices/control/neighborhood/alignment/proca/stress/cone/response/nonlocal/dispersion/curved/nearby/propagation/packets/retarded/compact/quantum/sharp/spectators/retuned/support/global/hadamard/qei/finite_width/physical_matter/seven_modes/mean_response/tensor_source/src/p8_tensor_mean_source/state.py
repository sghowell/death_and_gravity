"""Positive two-polarization band stress and the actual clock exchange."""

from functools import cache

import sympy as sp

from . import operator

a, q, T, P = operator.a, operator.q, operator.T, operator.P
H = sp.Symbol("physical_Hubble", real=True)


@cache
def data():
    M = sp.Matrix([[0, 2 / a**3], [-a * q / 2, 0]])
    R = sp.diag(q / (2 * a * a), 2 / a**6)
    S = sp.diag(-q / (6 * a * a), 2 / a**6)
    x = P * P / a**6
    y = q * T * T / (4 * a * a)
    r = x + y
    s = x - y / 3
    cross = q * T * P / a**5
    K = r - 3 * s
    Dt = lambda f: sp.factor(
        a * H * sp.diff(f, a)
        + 2 * P * sp.diff(f, T) / a**3
        - a * q * T * sp.diff(f, P) / 2
    )
    omega = sp.sqrt(q) / a
    return {
        "tensor_generator": M,
        "proxy_density_Hessian": R,
        "pressure_Hessian": S,
        "clock_kinetic_Hessian": R - 3 * S,
        "two_polarization_anchor_density": sp.trace(R.subs(a, 1)),
        "two_polarization_anchor_pressure": sp.trace(S.subs(a, 1)),
        "positive_band_state": "Add eta*b(k)*I2 to each of the two existing TT polarization covariances, and evolve on actual a on both legs",
        "Cartesian_TT_anchor_blocks": "diag(2*eta*b*Pi_TT,eta*b*Pi_TT/2), with vanishing mixed block",
        "two_polarization_initial_density_upper_per_eta": sp.Integer(4),
        "clock_kinetic_absolute_upper_per_proxy_density": sp.Integer(2),
        "clock_kinetic_derivative_bound": "abs(K_dot)<=(4*omega+12*abs(H))*r",
        "checks": {
            "two_polarization_proxy_stress_obeys_exact_evolution_conservation": (
                a * H * sp.diff(R, a) + M.T * R + R * M + 3 * H * (R + S)
            ).applyfunc(sp.factor),
            "positive_tensor_pressure_lower_keeps_kinetic_piece": R
            + 3 * S
            - sp.diag(0, 8 / a**6),
            "positive_tensor_pressure_upper_keeps_gradient_piece": R
            - S
            - sp.diag(2 * q / (3 * a * a), 0),
            "clock_kinetic_scalar_is_density_minus_three_pressure": sp.factor(
                K + 2 * x - 2 * y
            ),
            "clock_kinetic_time_derivative_keeps_transfer_cross_term": sp.factor(
                Dt(K) - 4 * cross - 12 * H * x + 4 * H * y
            ),
            "tensor_cross_term_bound_is_exact_square": sp.factor(
                r - cross / omega - (P / a**3 - sp.sqrt(q) * T / (2 * a)) ** 2
            ),
            "tensor_opposite_cross_term_bound_is_exact_square": sp.factor(
                r + cross / omega - (P / a**3 + sp.sqrt(q) * T / (2 * a)) ** 2
            ),
            "tensor_relative_covariance_reversal_is_anti_Hamiltonian": sp.diag(1, -1)
            * M
            * sp.diag(1, -1)
            + M,
            "tensor_actual_generator_is_even_in_clock": (
                M.subs(a, (1 + operator.u**2) ** 2).subs(operator.u, -operator.u)
                - M.subs(a, (1 + operator.u**2) ** 2)
            ).applyfunc(sp.factor),
            "tensor_proxy_density_equals_actual_on_clock_positive_Hamiltonian": sp.factor(
                r - operator.data()["positive_on_clock_tensor_energy"]
            ),
        },
    }


@cache
def variational_clock():
    vphi = sp.Symbol("homogeneous_clock_velocity", real=True)
    N, h = operator.N, operator.h
    V = operator.V
    ratio = 1 + (vphi * vphi / N**2 - 1) / h
    L = ratio * (a**3 * V * V / (4 * N) - N * a * q * T * T / 4)
    K = -V * V / 2 + q * T * T / (2 * a * a)
    return {
        "homogeneous_tensor_clock_momentum_on_background": -(a**3) * K / h,
        "checks": {
            "literal_clock_velocity_variation_is_actual_kinetic_scalar": sp.factor(
                sp.diff(L, vphi).subs({vphi: 1, N: 1}) + a**3 * K / h
            ),
            "explicit_clock_coefficient_derivative_vanishes_only_after_variation": sp.factor(
                sp.diff(L, h).subs({vphi: 1, N: 1})
            ),
        },
    }
