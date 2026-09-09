"""Actual natural-phase scalar source kernels and independent center observable bridge."""

from functools import cache

import sympy as sp
from p8_proca_physical_matter import center as physical_center

from . import model


@cache
def data():
    d = model.on_clock()
    a = (1 + model.u**2) ** 2
    fields = model.FIELDS
    kernels = {}
    for key in (
        "actual_mean_lapse_source",
        "mean_hat_scale_direct_source",
        "mean_trace_direct_source",
        "mean_matter_field_direct_source",
    ):
        kernels[key] = sp.hessian(d[key] / a**3, fields).applyfunc(sp.factor)
    H = sp.hessian(d["raw_density_quadratic_Hamiltonian"], fields).applyfunc(sp.factor)
    Omega = sp.zeros(4)
    Omega[:2, 2:] = sp.eye(2)
    Omega[2:, :2] = -sp.eye(2)
    return {
        "source_Hessians": kernels,
        "natural_phase_Hamiltonian_Hessian": H,
        "natural_phase_generator": Omega * H,
        "Omega": Omega,
        "checks": {
            "actual_natural_scalar_generator_preserves_CCR": Omega * H * Omega
            + Omega * (Omega * H).T
        },
    }


@cache
def center():
    d = data()
    values = {
        key: matrix.subs(model.u, 0).applyfunc(sp.factor)
        for key, matrix in d["source_Hessians"].items()
    }
    traces = {key: sp.factor(sp.trace(matrix) / 2) for key, matrix in values.items()}
    return {"source_Hessians": values, "source_traces_for_positive_I4_addition": traces}


@cache
def density_center():
    d = physical_center.data()
    fields = d["fields"]
    v, p, chi, pc = fields[:4]
    keep = {field: 0 for field in fields[4:]}
    additive = sp.factor(d["rho"].subs(keep, simultaneous=True))
    # v_add=zeta+zeta^2; p_add=P_zeta-2*zeta*P_zeta at H=0.
    # The latter has no linear density coefficient (alpha(0)=0).
    scalar_shift = -sp.Rational(3, 10) * v * v
    lapse_shift = (-scalar_shift / 20) / sp.Rational(243, 80)
    correction = scalar_shift / 10 - sp.Rational(3, 200) * lapse_shift
    natural = sp.factor(additive + correction)
    natural_matrix = sp.hessian(natural, (v, chi, p, pc)).applyfunc(sp.factor)
    trace = sp.factor(sp.trace(natural_matrix) / 2)
    natural_lapse = sp.factor(d["n2"].subs(keep, simultaneous=True) + lapse_shift)
    new_lapse = model.on_clock()["actual_mean_lapse_source"].subs(
        model.u, 0
    ) * sp.Rational(80, 243)
    k = next(iter(d["center_packet_map"].free_symbols))
    new_lapse = new_lapse.subs(
        {model.v: v, model.chi: chi, model.Pv: p, model.Pc: pc, model.k: k},
        simultaneous=True,
    )
    return {
        "additive_physical_matter_density_second_variation": additive,
        "required_scalar_chart_correction": correction,
        "natural_exponential_physical_density_jet": natural,
        "natural_density_trace_for_positive_I4_addition": trace,
        "independent_natural_second_order_lapse": natural_lapse,
        "checks": {
            "all_center_phase_pairs_of_scalar_mean_lapse_match_independent_constraint": sp.factor(
                natural_lapse - new_lapse
            ),
        },
    }
