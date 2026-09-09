"""Complete leading scalar mean forcing, without an unproved stress conservation transfer."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import mean as previous

from . import model, source

FN, FX, FP, FC = sp.symbols(
    "scalar_lapse_source scalar_scale_source scalar_trace_source scalar_matter_field_source",
    real=True,
)


@cache
def data():
    d = previous.data()
    xi, dp = previous.xi, previous.dp
    n = (d["alpha"] * dp - 3 * d["ell"] * d["beta"] * xi + FN) / (2 * d["J"])
    flow = sp.Matrix(
        [
            -dp / 2 + d["alpha"] * n / 3 + FX,
            -3 * d["H"] * dp + d["ell"] * d["beta"] * n - 3 * d["ell"] ** 2 * xi + FP,
        ]
    )
    field = d["beta"] * n - 3 * d["ell"] * xi + FC
    matrix = flow.jacobian((xi, dp))
    return {
        "mean_lapse": n,
        "mean_flow": flow,
        "induced_matter_field_derivative": field,
        "linear_physical_log_scale_response": xi + n / (2 * d["h"]),
        "checks": {
            "same_classical_mean_matrix_with_four_actual_scalar_sources": (
                matrix - d["two_component_mean_generator"]
            ).applyfunc(sp.factor),
            "actual_scalar_mean_lapse_constraint": sp.factor(
                -2 * d["J"] * n + d["alpha"] * dp - 3 * d["ell"] * d["beta"] * xi + FN
            ),
            "nonzero_scalar_scale_source_not_discarded": sp.factor(
                flow[0] + dp / 2 - d["alpha"] * n / 3 - FX
            ),
            "nonzero_scalar_matter_field_source_not_discarded": sp.factor(
                field - d["beta"] * n + 3 * d["ell"] * xi - FC
            ),
        },
    }


@cache
def center():
    d = source.data()
    M = d["natural_phase_generator"]
    u = model.u
    kernels = d["source_Hessians"]
    transport = lambda G: (sp.diff(G, u) + M.T * G + G * M).applyfunc(sp.factor)
    trace0 = lambda G: sp.factor(sp.trace(G.subs(u, 0)) / 2)
    f = {key: trace0(value) for key, value in kernels.items()}
    F2 = trace0(transport(transport(kernels["actual_mean_lapse_source"])))
    X1 = trace0(transport(kernels["mean_hat_scale_direct_source"]))
    bg = previous.data()
    J = bg["J"].subs(u, 0)
    J2 = sp.diff(bg["J"], u, 2).subs(u, 0)
    ap = sp.diff(bg["alpha"], u).subs(u, 0)
    ell, beta = bg["ell"].subs(u, 0), bg["beta"].subs(u, 0)
    n = f["actual_mean_lapse_source"] / (2 * J)
    pdot = ell * beta * n + f["mean_trace_direct_source"]
    xi2 = -pdot / 2 + ap * n / 3 + X1
    n2 = (2 * ap * pdot - 3 * ell * beta * xi2 + F2) / (2 * J) - n * J2 / J
    accel = sp.factor(xi2 + n2 / 2 - 11 * n)
    return {
        "initial_source_traces": f,
        "lapse_source_second_derivative": F2,
        "scale_source_first_derivative": X1,
        "anchor_induced_lapse": sp.factor(n),
        "anchor_trace_first_derivative": sp.factor(pdot),
        "anchor_hat_scale_second_derivative": sp.factor(xi2),
        "anchor_lapse_second_derivative": sp.factor(n2),
        "anchor_proper_Hubble_derivative_response": accel,
        "checks": {
            "zero_anchor_direct_scale_source_for_chosen_I4_addition": f[
                "mean_hat_scale_direct_source"
            ],
            "zero_anchor_lapse_source_first_derivative": trace0(
                transport(kernels["actual_mean_lapse_source"])
            ),
        },
    }


@cache
def parity():
    d = source.data()
    u = model.u
    R = sp.diag(1, -1, -1, 1)
    M = d["natural_phase_generator"]
    checks = {
        "actual_scalar_state_reversal_keeps_background_matter_sign": (
            M.subs(u, -u) * R + R * M
        ).applyfunc(sp.factor)
    }
    for key, G in d["source_Hessians"].items():
        sign = -1 if key == "mean_hat_scale_direct_source" else 1
        checks["actual_scalar_source_reflection_" + key] = (
            R * G.subs(u, -u) * R - sign * G
        ).applyfunc(sp.factor)
    return {"natural_phase_reversal": R, "checks": checks}
