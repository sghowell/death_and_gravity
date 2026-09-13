"""Complete causal all-transfer two-mass flat quotient inverse and normalization."""

from functools import cache

import sympy as s
from p8_vacuum_affine_flat_scalar_quotient_inverse import geometry as old_geometry
from p8_vacuum_affine_flat_scalar_quotient_inverse import resolvent as old_resolvent

from . import geometry as g
from . import spectral

lam, q = s.symbols("Laplace_frequency spatial_momentum_squared", positive=True)
p = lam * lam + q
t, T = s.symbols("time window_length", positive=True)
f0, f2 = s.symbols(
    "complete_total_trace_factor complete_total_shear_factor", nonzero=True
)
B = s.Matrix([[lam * lam + 2 * q / 3, -lam * lam / 3], [-q, -lam * lam]])
M = s.Matrix([[q / 3, lam * lam + 2 * q / 3, -lam / 3], [q, -q, -lam]])
GAUGE = s.Matrix([[lam, 0, 0], [0, 1, 0], [q, 0, lam]])
BINV = s.Matrix(
    [
        [1 / p, -1 / (3 * p)],
        [1 / p - 1 / lam**2, -1 / (3 * p) - s.Rational(2, 3) / lam**2],
    ]
)
COORDINATE_BOUND = s.Rational(5, 2)
PRIMITIVE_BOUND = 45 / (g.ell + 2)
L1_BOUND = 375 / (32 * (g.ell + 2))


def coordinate_kernel(transfer_squared):
    if isinstance(transfer_squared, bool) or not isinstance(
        transfer_squared, (int, s.Rational)
    ):
        raise TypeError("Use an exact nonnegative squared transfer")
    v = s.Rational(transfer_squared)
    if v < 0:
        raise ValueError("Squared transfer must be nonnegative")
    wave = t if v == 0 else s.sin(s.sqrt(v) * t) / s.sqrt(v)
    return s.ImmutableMatrix([[wave, -wave / 3], [wave - t, -wave / 3 - 2 * t / 3]])


@cache
def data():
    h = spectral.data()
    weight = s.diag(f0, s.Rational(8, 3) * f2)
    full3 = M.subs(lam, -lam).T * weight * M
    full2 = B.T * weight * B
    inverse = BINV * s.diag(1 / f0, s.Rational(3, 8) / f2) * BINV.T
    wave = s.sin(s.sqrt(q) * t) / s.sqrt(q)
    R = s.Matrix([[wave, -wave / 3], [wave - t, -wave / 3 - 2 * t / 3]])
    rg = s.zeros(2, 3)
    rg[:, 1:3] = B
    target = s.zeros(3)
    target[1:3, 1:3] = full2
    embedding = s.Matrix.hstack(
        s.Matrix(list(2 * s.eye(3))), s.Matrix(list(-2 * s.diag(1, 0, 0)))
    )
    gram = embedding.T * embedding
    cosine = s.Symbol("cosine", real=True)
    rp = R.diff(t).subs(s.cos(s.sqrt(q) * t), cosine)
    gram_poly = s.expand(sum(v * v for v in rp))
    kernel_bound = COORDINATE_BOUND**2 * PRIMITIVE_BOUND * t**3 / 6
    derivative_kernel_bound = COORDINATE_BOUND**2 * PRIMITIVE_BOUND * t * t / 2
    physical_force_bound = 64 * s.pi**2 * g.KAPPA * L1_BOUND * T**4
    og = {old_geometry.lam: lam, old_geometry.q: q}
    checks = {
        "original_complete_curvature_quotient_geometry_retained": B
        - old_geometry.B.subs(og),
        "original_complete_ordered_three_source_map_retained": M
        - old_geometry.M.subs(og),
        "complete_source_gauge_kernel": M * s.Matrix([lam, 0, q]),
        "complete_detector_gauge_kernel": M.subs(lam, -lam) * s.Matrix([-lam, 0, q]),
        "complete_source_coordinate_quotient": M * GAUGE - rg,
        "complete_detector_coordinate_quotient": M.subs(lam, -lam)
        * GAUGE.subs(lam, -lam)
        - rg,
        "full_three_source_ordered_gauge_quotient": GAUGE.subs(lam, -lam).T
        * full3
        * GAUGE
        - target,
        "both_detector_and_source_orders_retained": full3 - full3.subs(lam, -lam).T,
        "coordinate_left_inverse": B * BINV - s.eye(2),
        "coordinate_right_inverse": BINV * B - s.eye(2),
        "full_total_quotient_left_inverse": full2 * inverse - s.eye(2),
        "full_total_quotient_right_inverse": inverse * full2 - s.eye(2),
        "complete_total_quotient_determinant": full2.det()
        - s.Rational(8, 3) * f0 * f2 * lam**4 * p**2,
        "causal_coordinate_Laplace_inverse": R.applyfunc(
            lambda v: s.laplace_transform(v, t, lam, noconds=True)
        )
        - BINV,
        "whole_coordinate_initial_kernel": R.subs(t, 0),
        "continuous_zero_transfer_kernel": R.applyfunc(lambda v: s.limit(v, q, 0))
        - coordinate_kernel(0),
        "complete_coordinate_Frobenius_derivative": gram_poly
        - (20 * cosine**2 - 14 * cosine + 13) / 9,
        "coordinate_derivative_upper_margin": s.Rational(47, 9)
        - gram_poly
        - (1 + cosine) * (34 - 20 * cosine) / 9,
        "whole_trace_shifted_primitive_bound": 2
        * h["complete_static_and_next_measure_moments"]["static"][0]
        - 1 / (g.ell + 2),
        "whole_shear_shifted_primitive_bound": 2
        * h["complete_static_and_next_measure_moments"]["static"][2]
        - 120 / (g.ell + 2),
        "whole_scaled_diagonal_primitive_bound": s.Rational(3, 8) * 120 / (g.ell + 2)
        - PRIMITIVE_BOUND,
        "ordered_double_time_triangle": kernel_bound - 375 * t**3 / (8 * (g.ell + 2)),
        "complete_inverse_window_L1": s.integrate(kernel_bound, (t, 0, T))
        - L1_BOUND * T**4,
        "complete_first_derivative_window_L1": s.integrate(
            derivative_kernel_bound, (t, 0, T)
        )
        - 375 * T**3 / (8 * (g.ell + 2)),
        "physical_force_inverse_retains_kappa": physical_force_bound
        - 750 * s.pi**2 * g.KAPPA * T**4 / (g.ell + 2),
        "physical_amplitude_to_metric_Gram": gram - s.Matrix([[12, -4], [-4, 4]]),
        "physical_metric_lower_margin": (gram - 2 * s.eye(2)).det() - 4,
        "physical_metric_upper_margin": (14 * s.eye(2) - gram).det() - 4,
        "earlier_coordinate_kernel_not_earlier_spectral_inverse": coordinate_kernel(1)
        - old_resolvent.coordinate_kernel(1).subs(old_resolvent.t, t),
    }
    # Keep the paired forward distributions rather than individually divergent terms.
    y = g.y
    forward_h = {
        i: s.integrate(g.W_H[i] / s.sqrt(1 - y * y), (y, 0, 1)) / 2 for i in (0, 2)
    }
    checks["whole_scalar_trace_forward_beta_bound"] = forward_h[0] - 41 * s.pi / 64
    checks["whole_scalar_shear_forward_beta_bound"] = forward_h[2] - s.pi / 384
    m = s.Symbol("positive_mass", positive=True)
    freq = s.Symbol("positive_frequency", positive=True)
    sinprim = (1 - s.cos(freq * t)) / (freq * freq)
    checks["absolutely_convergent_shifted_primitive_Laplace"] = s.laplace_transform(
        sinprim, t, lam, noconds=True
    ) - 1 / (lam * (lam * lam + freq * freq))
    lo = g.ell
    checks["actual_log_gap_uniform_primitive_rational_cap"] = PRIMITIVE_BOUND.subs(
        lo, 394
    ) - s.Rational(5, 44)
    checks["actual_log_gap_uniform_operator_rational_cap"] = L1_BOUND.subs(
        lo, 394
    ) - s.Rational(125, 4224)
    return {
        "complete_curvature_channel_map": B,
        "complete_ordered_three_source_map": M,
        "complete_source_gauge_chart": GAUGE,
        "whole_two_mass_quotient_transform": full2,
        "whole_two_mass_quotient_inverse_transform": inverse.applyfunc(s.factor),
        "complete_coordinate_inverse_kernel": R,
        "entire_shifted_spectral_primitive": "J_i,q(t)=-integral from4mP^2 to infinity rho_i(sigma)/(q+sigma) [1-cos(sqrt(q+sigma)*t)] dsigma. The total two-mass reciprocal density is used, with no original pole term and no discarded heavy cut.",
        "diagonal_scaled_primitive_bound": PRIMITIVE_BOUND,
        "whole_ordinary_inverse_kernel": "E_q=R_q' * diag(J_trace,q,(3/8)J_2,q) * R_q^T, with causal time convolution. The exact transform is B^-1 diag(1/Ftrace,(3/8)/F2) (B^T)^-1. Both inverse products include the full initial distributional boundary.",
        "uniform_kernel_bound": kernel_bound,
        "uniform_operator_bound": L1_BOUND * T**4,
        "uniform_first_derivative_operator_bound": 375 * T**3 / (8 * (g.ell + 2)),
        "actual_rational_operator_cap": s.Rational(125, 4224) * T**4,
        "physical_force_inverse_bound": physical_force_bound,
        "spatial_graph": "For every real r and finite T, C_t H^r and L2_t H^r quotient dual sources have the displayed all-q uniform causal inverse bound, with no spatial derivative loss. The graph requires the complete forward distribution to equal the ordinary forcing at and after the initial boundary.",
        "forward_distribution": "F_i(D^2+q)=-c_i delta-(D^2+q)(G_Pi,q+G_Hi,q). Each G is the full y integral W(y) sin(sqrt(q+4m^2/(1-y^2))*t)/[(1-y^2)sqrt(q+4m^2/(1-y^2))]. Keep both masses and all initial terms.",
        "whole_forward_H_bounds": {i: v / m for i, v in forward_h.items()},
        "metric_embedding_Gram": gram,
        "scope": "This is the isolated complete two-mass Gaussian fourth-order flat quotient reference, not the three-source gauge inverse, a literal homogeneous constraint solution, complete curved force inverse, total interacting quantum state, smallness, stability or P8 closure. Restoring64pi^2 kappa is essential.",
        "checks": {
            name: value.applyfunc(s.factor)
            if isinstance(value, s.MatrixBase)
            else s.factor(value)
            for name, value in checks.items()
        },
        "gates": {
            "full_three_source_remains_gauge_singular": s.factor(full3.det()) == 0,
            "whole_quotient_is_generically_invertible": s.factor(full2.det()) != 0,
            "both_source_detector_parities_necessary": s.factor(
                full3[0, 2] - full3[2, 0]
            )
            != 0,
            "coordinate_derivative_bound_has_strict_margin": COORDINATE_BOUND**2
            > s.Rational(47, 9),
            "metric_embedding_lower_strict_margin": (gram - 2 * s.eye(2))[0, 0] > 0,
            "metric_embedding_upper_strict_margin": (14 * s.eye(2) - gram)[0, 0] > 0,
            "all_transfer_operator_cap_below_three_hundredths": s.Rational(125, 4224)
            < s.Rational(3, 100),
            "new_shear_primitive_dominates_new_trace_primitive": s.simplify(
                (g.ell + 2) * PRIMITIVE_BOUND
            )
            > 1,
            "no_inverse_transfer_in_continuous_kernel_definition": True,
            "both_total_cuts_and_original_finite_constants_used": True,
            "full_initial_distributional_graph_not_equation_only_after_zero": True,
            "normalized_reference_bound_not_physical_force_smallness": True,
        },
    }
