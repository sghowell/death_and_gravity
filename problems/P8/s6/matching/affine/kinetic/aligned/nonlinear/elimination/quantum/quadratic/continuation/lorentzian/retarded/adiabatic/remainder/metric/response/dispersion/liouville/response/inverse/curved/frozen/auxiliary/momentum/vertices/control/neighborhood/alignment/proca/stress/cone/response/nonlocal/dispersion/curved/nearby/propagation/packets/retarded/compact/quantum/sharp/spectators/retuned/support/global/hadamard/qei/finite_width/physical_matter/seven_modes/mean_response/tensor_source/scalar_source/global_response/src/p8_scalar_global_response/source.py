"""Exact compensated scalar mean kernels and global source integrals."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import mean
from p8_scalar_mean_source import geometry, model
from p8_scalar_mean_source import source as previous

from . import phase


@cache
def data():
    u = model.u
    t = 1 + u * u
    bg = mean.data()
    p0 = -2 * bg["H"]
    c = sp.Rational(9, 2) * p0
    B = geometry.canonical()["old_to_natural_phase"]
    raw = previous.data()["source_Hessians"]
    old = {key: (B.T * value * B).applyfunc(sp.factor) for key, value in raw.items()}
    V = sp.diag(2, 0, 0, 0)
    M = phase.data()["old_generator"]
    Vdot = (M.T * V + V * M).applyfunc(sp.factor)
    comp = {
        "lapse": (old["actual_mean_lapse_source"] + bg["alpha"] * c * V).applyfunc(
            sp.factor
        ),
        "scale": (old["mean_hat_scale_direct_source"] - c * V / 2).applyfunc(sp.factor),
        "trace": (
            old["mean_trace_direct_source"]
            - (3 * bg["H"] * c + sp.diff(c, u)) * V
            - c * Vdot
        ).applyfunc(sp.factor),
        "matter": old["mean_matter_field_direct_source"],
    }
    Ti = phase.data()["old_to_weighted_phase"].inv()
    kernels = {
        "mean_scale_forcing": bg["alpha"] * comp["lapse"] / (6 * bg["J"])
        + comp["scale"],
        "weighted_mean_trace_forcing": 100
        * t**3
        * (bg["ell"] * bg["beta"] * comp["lapse"] / (2 * bg["J"]) + comp["trace"]),
        "lapse_state_term": comp["lapse"] / (2 * bg["J"]),
        "matter_field_state_term": comp["matter"],
    }
    weighted = {
        key: (Ti.T * value * Ti).applyfunc(sp.factor) for key, value in kernels.items()
    }
    return {
        "compensating_trace_coefficient": c,
        "curvature_variance_Hessian": V,
        "curvature_variance_derivative_Hessian": Vdot,
        "compensated_old_phase_sources": comp,
        "actual_weighted_response_kernels": weighted,
    }


@cache
def envelopes():
    result = {}
    for key, matrix in data()["actual_weighted_response_kernels"].items():
        rows = []
        for i in range(4):
            for j in range(4):
                extra = int(i >= 2) + int(j >= 2)
                try:
                    row = phase.envelope(matrix[i, j], extra)
                except ValueError as exc:
                    raise ValueError(
                        key + " entry " + str((i, j)) + ": " + str(exc)
                    ) from exc
                rows.append(row)
        result[key] = {
            "entries": rows,
            "half_entry_L1_upper": sum(row["integral_upper"] for row in rows) / 2,
        }
    return result
