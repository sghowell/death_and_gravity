"""Independent center nonlinear spatial-chart correction to actual matter density."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import mean
from p8_proca_physical_matter import center


@cache
def data():
    d = center.data()
    fields = d["fields"]
    T, PT, S, PS = fields[4:8]
    k = next(iter(d["center_packet_map"].free_symbols))
    p = PT * PT + PS * PS - k * k * (T * T + S * S) / 12
    Ft = sp.Rational(3, 2) * p
    expected = -sp.Rational(2, 405) * Ft
    keep = {field: 0 for field in fields if field not in (T, PT, S, PS)}
    additive = sp.factor(d["rho"].subs(keep, simultaneous=True))
    v2 = (T * T + S * S) / 6
    dc = -sp.Rational(3, 10) * v2
    n1 = (-dc / 20) / (sp.Rational(243, 80))
    linear = dc / 10 - sp.Rational(3, 200) * n1
    return {
        "additive_tensor_density_second_variation": additive,
        "required_exponential_to_additive_scalar_mean_shift": v2,
        "induced_linear_homogeneous_scalar_density": linear,
        "actual_exponential_chart_density_response": expected,
        "checks": {
            "complete_nonlinear_spatial_chart_correction_recovers_tensor_mean_source": sp.factor(
                additive + linear - expected
            ),
            "omitted_scalar_mean_shift_is_the_exact_nonzero_constant_mismatch": sp.factor(
                additive - expected - sp.Rational(203, 40500) * (T * T + S * S)
            ),
        },
    }


@cache
def metric():
    xi, T, S, epsilon = sp.symbols(
        "homogeneous_log_scale tensor_plus tensor_cross perturbation_order", real=True
    )
    gamma = sp.Matrix([[T, S, 0], [S, -T, 0], [0, 0, 0]])
    trace = sp.trace(gamma * gamma)
    v = epsilon * xi + epsilon**2 * (xi**2 + trace / 12)
    tt = epsilon * gamma + epsilon**2 * (
        2 * xi * gamma + (gamma * gamma - sp.eye(3) * trace / 3) / 2
    )
    target = (
        sp.eye(3)
        + epsilon * (2 * xi * sp.eye(3) + gamma)
        + epsilon**2 * (2 * xi**2 * sp.eye(3) + 2 * xi * gamma + gamma * gamma / 2)
    )
    additive = (1 + 2 * v) * sp.eye(3) + tt
    determinant = sp.expand(additive.det())
    return {
        "nonlinear_additive_scalar_chart": v,
        "nonlinear_additive_TT_chart": tt,
        "exponential_metric_through_second_order": target,
        "chart_bridge_scope": "The second-order spatial metric and center mean-density bridge only; no all-time canonical momentum bridge is inferred",
        "checks": {
            "literal_second_order_exponential_spatial_metric_chart": (
                additive - target
            ).applyfunc(sp.expand),
            "second_order_chart_tensor_remainder_is_traceless": sp.expand(sp.trace(tt)),
            "second_order_determinant_has_no_pure_tensor_volume_term": sp.expand(
                determinant.coeff(epsilon, 2) - 18 * xi**2
            ),
            "first_order_determinant_is_homogeneous_log_scale_only": sp.expand(
                determinant.coeff(epsilon, 1) - 6 * xi
            ),
            "center_mean_density_does_not_depend_on_trace_momentum_chart_shift": mean.data()[
                "alpha"
            ].subs(mean.u, 0),
        },
    }
