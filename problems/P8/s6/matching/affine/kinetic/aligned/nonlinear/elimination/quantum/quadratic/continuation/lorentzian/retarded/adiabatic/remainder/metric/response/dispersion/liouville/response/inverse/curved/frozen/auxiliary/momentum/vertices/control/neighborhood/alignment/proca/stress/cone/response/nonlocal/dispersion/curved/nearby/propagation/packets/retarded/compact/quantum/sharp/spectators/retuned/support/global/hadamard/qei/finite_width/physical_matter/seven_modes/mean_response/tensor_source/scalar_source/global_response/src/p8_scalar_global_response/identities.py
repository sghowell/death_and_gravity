"""Independent variable-change, weighted-CCR and asymptotic controls."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import mean
from p8_scalar_mean_source import geometry, model
from p8_scalar_mean_source import source as original

from . import phase, source


@cache
def data():
    bg = mean.data()
    u = model.u
    t = 1 + u * u
    xi, d, V, Vp, FN, FX, FP, FC = sp.symbols(
        "xi compensated_trace variance variance_derivative FN FX FP FC", real=True
    )
    c = -9 * bg["H"]
    dp = d + c * V
    n_old = (bg["alpha"] * dp - 3 * bg["ell"] * bg["beta"] * xi + FN) / (2 * bg["J"])
    Fn = FN + bg["alpha"] * c * V
    Fx = FX - c * V / 2
    Fp = FP - (3 * bg["H"] * c + sp.diff(c, u)) * V - c * Vp
    n_new = (bg["alpha"] * d - 3 * bg["ell"] * bg["beta"] * xi + Fn) / (2 * bg["J"])
    x_old = -dp / 2 + bg["alpha"] * n_old / 3 + FX
    x_new = -d / 2 + bg["alpha"] * n_new / 3 + Fx
    d_old = (
        -3 * bg["H"] * dp
        + bg["ell"] * bg["beta"] * n_old
        - 3 * bg["ell"] ** 2 * xi
        + FP
        - sp.diff(c, u) * V
        - c * Vp
    )
    d_new = (
        -3 * bg["H"] * d + bg["ell"] * bg["beta"] * n_new - 3 * bg["ell"] ** 2 * xi + Fp
    )
    psi_old = bg["beta"] * n_old - 3 * bg["ell"] * xi + FC
    psi_new = bg["beta"] * n_new - 3 * bg["ell"] * xi + FC
    ph = phase.data()
    T = ph["old_to_weighted_phase"]
    K = ph["weighted_generator"]
    Omega = model.phase.data()["constant_symplectic_form"]
    transformed = T * Omega * T.T
    primitive = u + sp.Rational(2, 3) * u**3 + u**5 / 5
    checks = {
        "exact_compensated_lapse_is_same_physical_mean_lapse": sp.factor(n_old - n_new),
        "exact_compensated_scale_equation_is_original": sp.factor(x_old - x_new),
        "exact_compensated_trace_equation_keeps_variance_derivative": sp.factor(
            d_old - d_new
        ),
        "exact_compensated_matter_field_equation_is_original": sp.factor(
            psi_old - psi_new
        ),
        "compensated_trace_has_same_zero_anchor_data": c.subs(u, 0),
        "weighted_CCR_is_not_reset_to_canonical_identity": (
            transformed - phase.MOMENTUM_WEIGHT * Omega / t**3
        ).applyfunc(sp.factor),
        "weighted_generator_evolves_its_actual_time_dependent_CCR": (
            K * transformed + transformed * K.T - sp.diff(transformed, u)
        ).applyfunc(sp.factor),
        "exact_momentum_damping_integrating_factor": sp.factor(
            sp.diff(t**3, u) / t**3 - 6 * u / t
        ),
        "exact_polynomial_momentum_source_primitive": sp.expand(
            sp.diff(primitive, u) - t * t
        ),
        "momentum_source_primitive_positive_remainder": sp.expand(
            u * t * t - primitive - sp.Rational(4, 3) * u**3 - sp.Rational(4, 5) * u**5
        ),
        "improved_momentum_square_root_comparison": t - u * u - 1,
        "half_integer_integral_base": sp.integrate(
            (1 + u * u) ** (-sp.Rational(3, 2)), (u, 0, sp.oo)
        )
        - 1,
    }
    for n in range(2, 21):
        s = sp.Rational(2 * n + 1, 2)
        checks["half_integer_integral_recurrence_" + str(n)] = sp.factor(
            phase.half_line_integral(s)
            - (2 * s - 3) * phase.half_line_integral(s - 1) / (2 * s - 2)
        )
    reversal = sp.diag(1, -1, -1, 1)
    for name, kernel in source.data()["actual_weighted_response_kernels"].items():
        checks["compensated_kernel_is_symmetric_" + name] = kernel - kernel.T
        sign = -1 if name == "mean_scale_forcing" else 1
        checks["compensated_kernel_has_actual_reflection_" + name] = (
            kernel.subs(u, -u) - sign * reversal * kernel * reversal
        ).applyfunc(sp.factor)
    return {
        "same_mean_variable_change": "d=delta_p+(9H)*<v^2>; delta_p=d-9H*<v^2>",
        "checks": checks,
    }


def finite_infinite_limit(value):
    numerator, denominator = sp.fraction(sp.cancel(value))
    n, d = sp.Poly(numerator, model.u), sp.Poly(denominator, model.u)
    if n.is_zero or n.degree() < d.degree():
        return sp.Integer(0)
    if n.degree() != d.degree():
        raise ValueError("The requested rational limit is not finite")
    return sp.factor(n.LC() / d.LC())


@cache
def asymptotic_controls():
    u = model.u
    t = 1 + u * u
    bg = mean.data()
    B = geometry.canonical()["old_to_natural_phase"]
    Ti = phase.data()["old_to_weighted_phase"].inv()
    raw = original.data()["source_Hessians"]
    fN = B.T * raw["actual_mean_lapse_source"] * B
    fX = B.T * raw["mean_hat_scale_direct_source"] * B
    fP = B.T * raw["mean_trace_direct_source"] * B
    oldx = (Ti.T * (bg["alpha"] * fN / (6 * bg["J"]) + fX) * Ti)[0, 0]
    oldp = (
        100 * t**3 * Ti.T * (bg["ell"] * bg["beta"] * fN / (2 * bg["J"]) + fP) * Ti
    )[0, 0]
    new = source.data()["actual_weighted_response_kernels"]
    limits = {
        "uncompensated_scale_diagonal_times_u": finite_infinite_limit(u * oldx),
        "uncompensated_weighted_trace_diagonal_over_u_four": finite_infinite_limit(
            oldp / u**4
        ),
        "compensated_scale_diagonal_times_u": finite_infinite_limit(
            u * new["mean_scale_forcing"][0, 0]
        ),
        "compensated_weighted_trace_diagonal_over_u_four": finite_infinite_limit(
            new["weighted_mean_trace_forcing"][0, 0] / u**4
        ),
    }
    return {
        "actual_limits": limits,
        "no_physical_divergence_inferred_from_uncompensated_forcing": True,
        "checks": {
            "actual_uncompensated_scale_has_nonzero_nonintegrable_leading_term": limits[
                "uncompensated_scale_diagonal_times_u"
            ]
            - 396,
            "actual_uncompensated_weighted_trace_keeps_canonical_volume_growth": limits[
                "uncompensated_weighted_trace_diagonal_over_u_four"
            ]
            + 79200,
            "compensation_cancels_scale_diagonal_leading_term": limits[
                "compensated_scale_diagonal_times_u"
            ],
            "compensation_cancels_trace_diagonal_leading_term": limits[
                "compensated_weighted_trace_diagonal_over_u_four"
            ],
        },
    }
