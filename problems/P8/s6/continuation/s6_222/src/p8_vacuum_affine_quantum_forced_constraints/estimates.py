"""Exact force-interface constants and the unresolved compatible-space gap."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_reduced_scalar_hamiltonian import estimates as previous
from p8_vacuum_affine_scalar_tame_propagator import energy

from . import forces

PHASE_FORCE = s.Integer(400)
DIRECT = s.Integer(100)
RESPONSE = s.Integer(10) ** 117
METRIC_TRANSFER = s.Integer(2) * 10**5


@cache
def data():
    d = forces.system()
    p, scale = s.symbols("physical_transfer a", positive=True)
    phase_metric = d["C"].subs(forces.q, p * p / scale**2)
    F = d["F"].subs(forces.q, p * p / scale**2)
    f0 = F.subs(p, 0)
    f2 = F.applyfunc(lambda value: s.diff(value, p, 2) / 2)
    source = forces.old
    actual = source.actual_coefficients()
    mapping = {
        forces.th: actual[source.Theta],
        forces.E: actual[source.E],
        forces.l: source.ell,
        forces.w: actual[source.w],
        forces.J: actual[source.J0] + actual[source.dJ],
        forces.A: actual[source.A],
        forces.T: actual[source.Tc],
        forces.delta: source.delta,
        forces.H: source.H,
    }
    C0, C2 = previous.phase_matrices()
    actual_C = phase_metric.subs(mapping, simultaneous=True).subs(scale, source.a)
    # This data packet is a bridge to the established map constants; it does
    # not invent numerical high stress jets.
    density = RESPONSE / modes.KAPPA
    det_lower = s.Rational(1, 100)
    delta_upper = s.Rational(1, 2)
    direct_bound = (1 + delta_upper**2) / (2 * det_lower)
    r = s.Symbol("Sobolev_index", real=True)
    checks = {
        "actual_fixed_profile_metric_map_matches_S220": (
            actual_C - C0 - p * p * C2
        ).applyfunc(s.factor),
        "same_physical_transfer_phase_metric_force_norm": (
            F + forces.Jcan * phase_metric.T
        ).applyfunc(s.factor),
        "phase_force_degree_two_in_transfer": (F - f0 - p * p * f2).applyfunc(s.factor),
        "exact_output_density_response_constant": density - s.Integer(10) ** -683,
        "direct_rank_two_eigenvalue_majorant": direct_bound - s.Rational(125, 2),
        "force_to_phase_derivative_loss": (r - 2) - energy.SPATIAL_LOSS - (r - 14),
        "force_to_metric_derivative_loss": (r - 14) - 2 - (r - 16),
        "scalar_response_to_phase_target": (r - 3) - 14 - (r - 17),
        "scalar_response_to_metric_target": (r - 3) - 16 - (r - 19),
        "phase_input_to_first_response_spatial_loss": (r + 10) - (r - 17) - 27,
        "scalar_input_to_first_phase_spatial_loss": (r + 8) - (r - 17) - 25,
        "scalar_input_to_next_scalar_spatial_loss": (r + 8) - (r - 19) - 27,
    }
    return {
        "phase_force_matrix_F0": f0,
        "phase_force_matrix_F2": f2,
        "phase_force_amplitude_bound": PHASE_FORCE,
        "direct_auxiliary_amplitude_bound": DIRECT,
        "normalized_scalar_response_bound": density,
        "classical_propagator_log_constant": energy.LOG_PROPAGATOR,
        "full_force_to_metric_constant_before_exponential": METRIC_TRANSFER,
        "known_interface_estimates": "||F(t,P)||<400(1+|P|^2) follows from F=-Jcan C^T and the exact S220 C bound. ||D||<=125/2<100 follows from its two nonzero eigenvalues. Thus G0F:L2H^s->C H^(s-14) has bound400 exp(1e29), and T=D+C G0F:L2H^s->L2H^(s-16) has bound2e5 exp(1e29).",
        "Riesz_and_spatial_shift": "The complete S220 detector bound has no time derivatives, so Riesz gives Rhat:X_0->L2H^-3. Spatial translation invariance makes it commute with every Bessel weight, giving Qbar:X_r->Y_r with bound1e-683 for each real r. Output multiplication by a^-3 costs at most1 and no spatial derivatives.",
        "first_response_statement": "For a prescribed S0 inX_r, g1=Qbar S0 is well defined. Z1=G0F g1 lies inC H^(r-17) with bound exp(1e29)1e-680 ||S0||X_r, and S1=Tg1 lies inL2H^(r-19) with bound exp(1e29)1e-677 ||S0||X_r. These are constructed first responses, not solutions of the feedback equation.",
        "comparison_phase_input": "For S0=C Zref, the unchanged explicitly defined C13 gives ||S0||X_r<=C13 ||Zref||H13_t H^(r+10). The first phase response therefore loses27 spatial derivatives in these displayed norms. No numerical thirteen-jet stress estimate is assumed.",
        "precise_estimate_gap": "The available scalar feedback output is only L2_t H^(r-19), whereas its required source input is H13_t H^(r+8). The stated bounds do not give a bounded self-map or justify Neumann iteration. This is a limitation of the current estimates, not proof that the actual response lacks better regularity or that the candidate is unstable.",
        "normalization_scope": "The explicit1/kappa here is derived for the stated scalar density force interface. It is not transferred to an unspecified canonical noise or physical-source norm. Even this derived factor gives no smallness claim after exp(1e29) and the unevaluated C13.",
        "checks": checks,
        "gates": {
            "direct_auxiliary_uniform_bound": direct_bound < DIRECT,
            "metric_transfer_constant": DIRECT + PHASE_FORCE**2 < METRIC_TRANSFER,
            "first_phase_constant": PHASE_FORCE < s.Integer(10) ** 3,
            "first_metric_constant": METRIC_TRANSFER < s.Integer(10) ** 6,
            "original_fixed_kappa": modes.KAPPA == s.Integer(10) ** 800,
            "original_twelve_derivative_classical_comparison": energy.SPATIAL_LOSS
            == 12,
            "original_classical_log_constant": energy.LOG_PROPAGATOR
            == s.Integer(10) ** 29,
            "S220_C13_definition_not_replaced": callable(previous.phase_matrices),
            "no_time13_self_map_from_L2_estimate": True,
            "no_finite_coupling_remainder_or_quantum_stability": True,
        },
    }
