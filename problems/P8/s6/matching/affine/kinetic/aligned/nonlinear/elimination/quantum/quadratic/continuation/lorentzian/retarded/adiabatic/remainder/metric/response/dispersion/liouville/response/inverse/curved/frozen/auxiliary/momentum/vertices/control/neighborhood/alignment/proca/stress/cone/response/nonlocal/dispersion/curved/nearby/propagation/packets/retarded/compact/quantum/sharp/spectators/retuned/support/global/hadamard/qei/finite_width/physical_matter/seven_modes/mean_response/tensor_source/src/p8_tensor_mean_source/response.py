"""Actual tensor-forced homogeneous mean, not an ordinary Proca transfer."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import mean as prior

from . import operator, state


@cache
def data():
    d = prior.data()
    u = prior.u
    r, s = prior.rho, prior.pressure
    xi, dp = prior.xi, prior.dp
    F = (1 - 1 / d["h"]) * r + 3 * s / (2 * d["h"])
    n = sp.factor((d["alpha"] * dp - 3 * d["ell"] * d["beta"] * xi + F) / (2 * d["J"]))
    flow = sp.Matrix(
        [
            -dp / 2 + d["alpha"] * n / 3,
            -3 * d["H"] * dp + d["ell"] * d["beta"] * n - 3 * d["ell"] ** 2 * xi + s,
        ]
    )
    variables = sp.Matrix([xi, dp])
    A = flow.jacobian(variables).applyfunc(sp.factor)
    b = flow.subs({xi: 0, dp: 0}).applyfunc(sp.factor)
    B = xi + n / (2 * d["h"])
    reflection = {xi: xi, dp: -dp, r: r, s: s}
    return {
        "actual_tensor_lapse_forcing": F,
        "actual_tensor_induced_lapse": n,
        "actual_tensor_mean_flow": flow,
        "same_mean_generator": A,
        "new_tensor_mean_forcing": b,
        "physical_linear_log_scale_response": B,
        "physical_matter_density_and_pressure_response": -3 * d["ell"] ** 2 * B,
        "induced_matter_field_derivative": d["beta"] * n - 3 * d["ell"] * xi,
        "checks": {
            "actual_tensor_source_does_not_change_classical_mean_operator": (
                A - d["two_component_mean_generator"]
            ).applyfunc(sp.factor),
            "actual_tensor_mean_flow_retains_complete_new_forcing": (
                flow - A * variables - b
            ).applyfunc(sp.factor),
            "actual_tensor_lapse_constraint_retains_clock_source": sp.factor(
                -2 * d["J"] * n + d["alpha"] * dp - 3 * d["ell"] * d["beta"] * xi + F
            ),
            "tensor_forcing_difference_is_actual_clock_kinetic_coupling": sp.factor(
                F - d["lapse_force"] + (r - 3 * s) / d["h"]
            ),
            "tensor_mean_scale_equation_has_reflection_parity": sp.factor(
                flow[0].subs(u, -u).subs(reflection, simultaneous=True) + flow[0]
            ),
            "tensor_mean_trace_equation_has_reflection_parity": sp.factor(
                flow[1].subs(u, -u).subs(reflection, simultaneous=True) - flow[1]
            ),
            "tensor_mean_lapse_has_reflection_parity": sp.factor(
                n.subs(u, -u).subs(reflection, simultaneous=True) - n
            ),
        },
    }


@cache
def ward():
    d = prior.data()
    u = prior.u
    r, s = sp.Function("tensor_proxy_density")(u), sp.Function("tensor_pressure")(u)
    B = sp.Function("linear_physical_log_scale")(u)
    K = r - 3 * s
    clock = sp.diff(K / d["h"], u) + 3 * d["H"] * K / d["h"]
    physical = r - K / d["h"]
    matter = -3 * d["ell"] ** 2 * B
    connection = 3 * d["ell"] ** 2 * sp.diff(B, u)
    total = (
        sp.diff(physical + matter, u)
        + 3 * d["H"] * (physical + s + 2 * matter)
        + connection
        + clock
    )
    scale = (1 + u * u) ** 2
    return {
        "actual_tensor_clock_source": clock,
        "actual_tensor_fixed_metric_density": physical,
        "actual_induced_matter_density_and_pressure": matter,
        "retained_background_matter_connection_variation": connection,
        "volume_weighted_clock_boundary_function": scale**3 * K / d["h"],
        "checks": {
            "induced_matter_Ward_identity_retains_connection_variation": sp.factor(
                sp.diff(matter, u) + 6 * d["H"] * matter + connection
            ),
            "joint_tensor_and_induced_matter_Ward_identity_retains_actual_clock_source": sp.factor(
                total.subs(sp.diff(r, u), -3 * d["H"] * (r + s))
            ),
            "actual_volume_weighted_clock_exchange_is_exact_boundary_derivative": sp.factor(
                scale**3 * clock - sp.diff(scale**3 * K / d["h"], u)
            ),
        },
    }


@cache
def center():
    d = prior.data()
    u = prior.u
    a, q = operator.a, operator.q
    H = state.H
    Hd = sp.Symbol("physical_Hubble_derivative", real=True)
    st = state.data()
    M = st["tensor_generator"]

    def transport(G):
        return (a * H * sp.diff(G, a) + Hd * sp.diff(G, H) + M.T * G + G * M).applyfunc(
            sp.factor
        )

    at = {a: 1, H: 0, Hd: 4}
    R, S = st["proxy_density_Hessian"], st["pressure_Hessian"]
    r, s = (sp.factor(sp.trace(G.subs(at))) for G in (R, S))
    r2, s2 = (sp.factor(sp.trace(transport(transport(G)).subs(at))) for G in (R, S))
    F = sp.Rational(3, 2) * s
    F2 = 6 * r + sp.Rational(3, 2) * s2 - 9 * s
    J = d["J"].subs(u, 0)
    J2 = sp.diff(d["J"], u, 2).subs(u, 0)
    ap = sp.diff(d["alpha"], u).subs(u, 0)
    ell, beta = d["ell"].subs(u, 0), d["beta"].subs(u, 0)
    n = F / (2 * J)
    pdot = s + ell * beta * n
    xiddot = -pdot / 2 + ap * n / 3
    nddot = sp.factor(
        (2 * ap * pdot - 3 * ell * beta * xiddot + F2) / (2 * J) - n * J2 / J
    )
    acceleration = sp.factor(xiddot + nddot / 2 - 11 * n)
    clock_derivative = sp.factor(r2 - 3 * s2 + 6 * (r - 3 * s))
    return {
        "anchor_proxy_density": r,
        "anchor_pressure": s,
        "anchor_proxy_density_second_derivative": r2,
        "anchor_pressure_second_derivative": s2,
        "anchor_lapse_forcing": F,
        "anchor_lapse_forcing_second_derivative": F2,
        "anchor_lapse_response": sp.factor(n),
        "anchor_hat_scale_second_derivative": sp.factor(xiddot),
        "anchor_lapse_response_second_derivative": nddot,
        "anchor_proper_Hubble_derivative_response": acceleration,
        "anchor_clock_source_first_derivative": clock_derivative,
        "checks": {
            "new_tensor_proxy_stress_center_conservation": sp.factor(r2 + 12 * (r + s)),
            "new_tensor_source_center_keeps_both_clock_coefficients": sp.factor(
                F2 - (2 * q * q - sp.Rational(3, 2) * q - 78)
            ),
            "actual_retuned_mean_pivot_second_derivative": J2 + sp.Rational(633, 200),
            "actual_selected_tensor_band_has_nonzero_clock_source_derivative": sp.factor(
                clock_derivative - (72 + 14 * q - 4 * q * q)
            ),
            "actual_clock_source_derivative_lower_bound_factorization": sp.factor(
                clock_derivative - 64 - 2 * (2 * q + 1) * (4 - q)
            ),
            "actual_clock_source_derivative_upper_bound_square": sp.factor(
                sp.Rational(337, 4)
                - clock_derivative
                - 4 * (q - sp.Rational(7, 4)) ** 2
            ),
        },
    }
