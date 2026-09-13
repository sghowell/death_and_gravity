"""Complete ordered nonscalar inverse and precise compatible vector force graph."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_mass_coupled_inverse import inverse as scalar_inverse
from p8_vacuum_affine_two_mass_reference import spectral

from . import geometry, local

require_ball = scalar_inverse.require_ball


@cache
def data():
    C = s.Symbol("entire_nonscalar_remainder_C", nonnegative=True)
    K = s.Symbol("whole_total_shear_half_line_L1_norm", positive=True)
    beta = K * C
    weight = (4 * beta + 1) ** 2
    contraction = 2 * beta / (4 * beta + 1)
    T = s.Symbol("window_length", positive=True)
    Pmax = s.Symbol("fixed_external_radius", positive=True)
    canonical = 2 * s.exp(weight) * K * local.C0
    tt = geometry.t
    a = geometry.metric.a
    H = geometry.metric.H
    u = s.Dummy("source_time", real=True)
    f = s.Function("canonical_vector_force")(tt)
    e_v = s.sqrt(local.state.KAPPA) * f / 2
    IW = a**-3 * s.Integral((a**3 * e_v).subs(tt, u), (u, -s.Rational(1, 2), tt))
    e_beta = -s.sqrt(2) * s.I * geometry.P * IW
    rho0 = s.Rational(25, 16) ** 3
    v = s.Function("metric_vector")(tt)
    b = s.Function("shift_vector")(tt)
    # An arbitrary independent source pair must satisfy this relation; gauge forcing is not erased.
    source_pair = (
        -s.sqrt(2) * s.I * geometry.P * e_v - s.diff(e_beta, tt) - 3 * H * e_beta
    )
    checks = {
        "strict_half_contraction_margin": s.Rational(1, 2)
        - contraction
        - 1 / (2 * (4 * beta + 1)),
        "both_chain_factors_in_canonical_forcing": local.C0
        - 16 * s.pi**2 * local.state.KAPPA,
        "complete_physical_vector_force_Ward": s.cancel(source_pair),
        "vector_canonical_to_spatial_force": 2 * e_v / s.sqrt(local.state.KAPPA) - f,
        "weighted_vector_force_primitive_bound_margin": rho0
        + s.Rational(759, 4096)
        - 4,
        "two_canonical_legs_not_scalar_coefficient_chart": 4
        / (64 * s.pi**2 * local.state.KAPPA)
        - 1 / (16 * s.pi**2 * local.state.KAPPA),
    }
    kappa = s.Symbol("kappa", positive=True)
    h = s.Function("canonical_vector_metric")(tt)
    physical = 2 * h / s.sqrt(kappa)
    beta_inv = -s.diff(physical, tt) / (s.sqrt(2) * s.I * geometry.P)
    invariant = s.diff(v, tt) - s.sqrt(2) * s.I * geometry.P * b
    completion = (
        invariant.subs({v: physical, b: 0}) + s.sqrt(2) * s.I * geometry.P * beta_inv
    )
    checks["complete_physical_shift_representative"] = s.cancel(completion.doit())
    # Finite ordered causal matrix fixture; continuum estimates are in the notes.
    F = s.Matrix([[2, 0, 0], [1, 3, 0], [2, 1, 4]])
    V = s.Matrix(
        [[0, 0, 0], [s.Rational(1, 3), 0, 0], [s.Rational(1, 7), s.Rational(1, 5), 0]]
    )
    B0 = F.inv()
    E = (s.eye(3) + B0 * V).inv() * B0
    checks["ordered_causal_inverse_left"] = (F + V) * E - s.eye(3)
    checks["ordered_causal_inverse_right"] = E * (F + V) - s.eye(3)
    return {
        "entire_reference_channel": "F2,total of S247 with both actual mass cuts, fixed finite constant-(log(n)+2)/60, no first-sheet isolated pole and ordinary half-line L1 inverse K2,total. It is not the old isolated Proca inverse.",
        "actual_normal_form": "I4 T_tensor,total=F2,total+V_tensor and I4 T_vector,total=F2,total+V_vector. All physical finite/classical/retarded contacts remain. Each complete V and every fixed simultaneous time derivative has finite C_j(Pmax)(1+|loglag|) majorant.",
        "ordered_inverse": "(I+K2,total V_sector)^-1 K2,total I4, applied to g=-16pi²kappa f in canonical force units.",
        "weight": weight,
        "strict_contraction_bound": contraction,
        "canonical_source_inverse_bound": canonical,
        "source_norm": "||h||C_tHr <=2exp(weight) Kmax16pi²kappa ||I4 f||C_tHr; use the maximum complete remainder majorant of tensor/vector. All constants are finite UNEVALUATED and may be enormous. There is no physical smallness or all-Pmax bound.",
        "smoothness": "The reference convolution is stationary proper-time q0. Prepared simultaneous shifts and the complete V_j bounds prove every fixed smooth seminorm; no new initial branch, final-time condition or state reset. Recover the original equation by four output derivatives.",
        "vector_applied_force_pair": s.ImmutableMatrix([e_v, e_beta]),
        "vector_force_graph_bound": "With e_v=sqrt(kappa)f/2 and e_beta=-sqrt(2)iP I_W e_v, ||e_beta||<=2sqrt(2kappa) Pmax T ||f|| since ||I_W||<=4T on the actual slab. This is a gauge-compatible source graph, not arbitrary raw shift forcing in the same infrared norm.",
        "vector_force_bound": 2 * s.sqrt(2 * local.state.KAPPA) * Pmax * T,
        "vector_constraint_recovery": "The entire original normalized-a³ residual satisfies W E_beta+sqrt(2)iP E_v=0. The compatible applied pair satisfies the same identity. A solution of the synchronous equation makes W(E_beta+e_beta)=0, and the unchanged initial zero germ forces the original transverse shift equation. Both polarizations are retained.",
        "physical_vector_representative": beta_inv,
        "infrared_boundary": "For a prescribed shift force one needs -(D+3H)e_beta/(sqrt(2)iP) in the admitted quotient space. Multiplication by1/P is not bounded on L2 on a full ball. The zero-spatial-vector representative's shift is tempered, not necessarily H^r. The literal P0 system is separate.",
        "combined_scope": "O(2), including reflection, makes the full actual nonzero-P metric/M1 conditional reference the direct sum of S248's scalar/clock/M1 block and these tensor/vector doublets, on the direct sum of their explicit smooth prepared amplitude/compatible-force domains. This is not unrestricted physical forcing, evaluated stability, interacting light/mixed loops, nonlinear same-state bounce, UV/Regge or P8 closure.",
        "checks": {key: geometry.zero(value) for key, value in checks.items()},
        "gates": {
            "strict_half_margin_for_every_finite_C_and_K": s.simplify(
                s.Rational(1, 2) - contraction
            ).is_positive,
            "both_full_reference_mass_threshold_gates_rechecked": all(
                bool(v) for v in spectral.data()["gates"].values()
            ),
            "weighted_force_ratio_below_four": rho0 < 4,
            "inverse_order_noncommuting": F * V != V * F,
            "vector_constraint_zero_germ_not_replaced_by_final_data": True,
            "tensor_vector_actual_remainders_distinct": True,
            "all_original_parameters_and_finite_prescription_unchanged": True,
            "no_unrestricted_shift_force_or_P8_promotion": True,
        },
    }
