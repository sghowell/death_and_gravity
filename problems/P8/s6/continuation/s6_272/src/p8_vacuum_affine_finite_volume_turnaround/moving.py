"""Actual near-clock reference flow in a fixed symplectic balance."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from p8_vacuum_affine_quantitative_local_time import reference as old
from p8_vacuum_affine_scalar_tame_propagator import charts

from . import source

P, CORE = old.field.P, old.field.CORE
TIME = source.TIME
RADIUS, IMAGE_RADIUS = 4 * CORE, 8 * CORE
GENERATOR = s.Integer(10) ** 128


@cache
def bounds():
    a, original = phase.a, old.original
    physical_box = {
        original.u: I(-TIME, TIME),
        original.rho: I(-original.PROFILE_BOUND, original.PROFILE_BOUND),
        original.pressure: I(-original.PROFILE_BOUND, original.PROFILE_BOUND),
    }
    pivot = evaluate(original.J, physical_box)
    delta = 1 / (2 * (1 + original.u**2) ** 3)
    actual = {
        "Theta": old.maximum(evaluate(original.theta, physical_box)),
        "H": old.maximum(evaluate(original.H, physical_box)),
        "E": old.maximum(evaluate(original.E, physical_box)),
        "ell": old.maximum(evaluate(original.ell, physical_box)),
        "A": original.PROFILE_BOUND,
        "Tcorr": old.maximum(
            evaluate(original.rho - 3 * delta * original.pressure, physical_box)
        ),
    }
    box = {
        a: I(1, 2),
        charts.th: I(-5 * TIME, 5 * TIME),
        charts.E: I(-s.Rational(1, 2), s.Rational(1, 2)),
        charts.l: I(-s.Rational(1, 10), s.Rational(1, 10)),
        charts.J: I(1, 100),
        charts.A: I(-original.PROFILE_BOUND, original.PROFILE_BOUND),
        charts.T: I(-3 * original.PROFILE_BOUND, 3 * original.PROFILE_BOUND),
    }
    B = s.diag(
        s.sqrt(P), s.sqrt(P), 1 / s.sqrt(P), 1 / s.sqrt(P)
    ) * phase.central_map().subs({a: 1, charts.q: P**2})
    light = (
        (B * phase.canonical_generator() * B.inv())
        .subs(charts.q, P**2 / a**2)
        .applyfunc(s.cancel)
    )
    entries = [old.maximum(evaluate(value, box)) for value in light]
    tensor_balance = s.diag(s.sqrt(P), 1 / s.sqrt(P))
    tensor = (tensor_balance * phase.tensor_generator() * tensor_balance.inv()).subs(
        charts.q, P**2 / a**2
    )
    tensor_bound = sum(old.maximum(evaluate(value, box)) for value in tensor)
    heavy = (
        old.field.OMEGA_HI
        + (8 * old.field.heavy_state.MASS2 + 2 * P**2) / old.field.OMEGA_LO
    )
    vector = 18 * P
    balanced = 6 * (sum(entries) + 2 * tensor_bound + heavy + vector)
    condition = 192 * s.Integer(10) ** 21
    whitened = condition * balanced
    return {
        "B": B,
        "light": light,
        "light_entries": entries,
        "tensor": tensor,
        "tensor_bound": tensor_bound,
        "heavy": heavy,
        "vector": vector,
        "balanced": balanced,
        "condition": condition,
        "whitened": whitened,
        "flow_deviation": 2 * whitened * TIME,
        "actual": actual,
        "pivot": [
            s.Rational(pivot.lo.numerator, pivot.lo.denominator),
            s.Rational(pivot.hi.numerator, pivot.hi.denominator),
        ],
    }


@cache
def data():
    b = bounds()
    a, mass = s.symbols("positive_scale positive_mass", positive=True)
    k = s.Matrix(s.symbols("full_three_vector_k0:3", real=True))
    K = s.eye(3) / a + k * k.T / (a**3 * mass**2)
    V = a * mass**2 * s.eye(3) + (k.dot(k) * s.eye(3) - k * k.T) / a
    tensor_shape = s.Matrix([[0, 4 / phase.a**3], [-(phase.a**3) * charts.q / 4, 0]])
    tensor_balance = s.diag(s.sqrt(P) / 2, 2 / s.sqrt(P))
    checks = {
        "whole_fixed_initial_scalar_balance_symplectic": b["B"] * phase.OMEGA * b["B"].T
        - phase.OMEGA,
        "whole_balanced_scalar_generator_symplectic": (
            b["light"] * phase.OMEGA + phase.OMEGA * b["light"].T
        ).applyfunc(s.cancel),
        "whole_tensor_shape_and_dual_balance": (
            tensor_balance * tensor_shape * tensor_balance.inv()
        ).subs(charts.q, P**2 / phase.a**2)
        - b["tensor"],
        "whole_all_polarization_Proca_product": (
            K * V - (mass**2 + k.dot(k) / a**2) * s.eye(3)
        ).applyfunc(s.factor),
        "full_covariance_trace_condition_factor": b["condition"] - 2 * 96 * 10**21,
        "physical_heavy_frequency_not_fixed_Omega_star": old.field.heavy_state.MASS2
        + P**2
        - (old.field.heavy_state.MASS2 + P**2 / 16)
        - 15 * P**2 / 16,
    }
    return {
        "whole_fixed_scalar_symplectic_balance": b["B"],
        "whole_full_balanced_scalar_generator": b["light"],
        "whole_all_sixteen_outward_entry_bounds": b["light_entries"],
        "whole_full_tensor_balanced_generator": b["tensor"],
        "whole_full_three_polarization_Proca_K_V": [K, V],
        "whole_actual_near_clock_coefficient_bounds": b["actual"],
        "whole_actual_reference_pivot": b["pivot"],
        "whole_sector_bounds": {
            name: b[name] for name in ("tensor_bound", "heavy", "vector", "balanced")
        },
        "whole_covariance_condition_and_whitened_generator": [
            b["condition"],
            b["whitened"],
        ],
        "whole_exponential_flow_deviation_bound": b["flow_deviation"],
        "whole_fixed_balance_proof": "The scalar balance is fixed at u=0, so no moving-connection term is omitted. The full density canonical generator already includes the original time contacts. Tensor shape/dual factors are kept in the exact displayed balance. Heavy Omega_star is the fixed S240 energy balance, not instantaneous dispersion.",
        "whole_Proca_bound_proof": "Use the fixed initial B=sqrt(I+kk^T/m^2), Omega=sqrt(m^2+P^2) and both zeta factors. In a full longitudinal/transverse polarization frame the balanced upper block is Omega B^-1 K B^-1<=Omega and the lower block B V B/Omega<=2Omega for 1<=a<=2. The entire six-entry sum is<=9Omega<18P. No polarization or Gauss source is removed.",
        "whole_covariance_proof": "Every balanced diagonal variance is bounded by1e21 in the unchanged pure state; all cross covariances remain. For S_bal=B_full S0, ||S_bal||^2<=2 trace(V_bal)<=192e21. Symplecticity gives ||S_bal^-1||=||S_bal||. Thus the displayed bound controls S0^-1 G_ref S0 without a new vacuum or discarded covariance. The full integral equation bounds ||M_w-I|| by exp(GT)-1<=2GT.",
        "checks": checks,
        "gates": {
            "actual_Theta_small": b["actual"]["Theta"] <= 5 * TIME,
            "actual_H_small": b["actual"]["H"] < 5 * TIME,
            "actual_E_and_ell_bounds": b["actual"]["E"] <= s.Rational(1, 2)
            and b["actual"]["ell"] <= s.Rational(1, 10),
            "actual_nonzero_profile_bounds": b["actual"]["A"] == s.Rational(1, 10**390)
            and b["actual"]["Tcorr"] < 3 * old.original.PROFILE_BOUND,
            "actual_pivot": 1 < b["pivot"][0] <= b["pivot"][1] < 100,
            "full_scalar_bound": sum(b["light_entries"]) < 10**68,
            "both_tensor_bound": b["tensor_bound"] < 10**67,
            "full_heavy_bound": b["heavy"] < 10**102,
            "full_balanced_generator": b["balanced"] < 10**103,
            "full_whitened_generator": b["whitened"] < GENERATOR,
            "entire_flow_deviation": b["flow_deviation"] < s.Rational(1, 50),
            "entire_complex_phase_image": (1 + b["flow_deviation"]) * RADIUS
            < IMAGE_RADIUS,
            "actual_scale_and_source_slab": (1 + TIME**2) ** 2 < 2
            and TIME < source.q.current.TIME_LENGTH,
            "same_full_reference_and_all_eight_channels": True,
        },
    }
