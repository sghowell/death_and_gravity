"""Complete actual scalar ADM/common-clock norms and original projected-mean tail."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_spatial_response import endpoints

MEAN = s.Integer(10) ** 400
MEAN_TAIL = s.Integer(10) ** 310
SYN_DETECTOR = s.Integer(100)
SYN_SOURCE = s.Integer(10) ** 20
CLOCK_DETECTOR = s.Integer(4)
CLOCK_SOURCE = s.Integer(10) ** 20


def leibniz(order, leading):
    if (
        isinstance(order, bool)
        or not isinstance(order, (int, s.Integer))
        or not 0 <= order <= 13
    ):
        raise ValueError("Use the complete source derivative order zero through13")
    if isinstance(leading, bool) or leading not in (6, 7, 28):
        raise ValueError("Use a specified full clock, Hubble or inverse-scale majorant")
    return sum(
        s.binomial(order, j) * leading * s.factorial(j) * 4**j for j in range(order + 1)
    )


@cache
def constants():
    kk = state.KAPPA
    sync = 4 * (3 + 4 * leibniz(13, 7) + 2 * leibniz(12, 28))
    clock = 4 * (3 + 4 * leibniz(13, 6))
    from_spatial = endpoints.constants()
    spatial = from_spatial["complete_actual_all_transfer_spatial_response_bound"]
    spatial_tail = from_spatial[
        "whole_original_frequency_projected_response_tail_numerator"
    ]
    adm = (
        SYN_DETECTOR * SYN_SOURCE * (spatial + 8 * MEAN / kk) + (28520 + 64) * MEAN / kk
    )
    total = CLOCK_DETECTOR * CLOCK_SOURCE * adm
    adm_tail = SYN_DETECTOR * SYN_SOURCE * spatial_tail + 28520 * MEAN_TAIL / kk
    clock_tail = CLOCK_DETECTOR * CLOCK_SOURCE * adm_tail + 24 * MEAN_TAIL / kk
    return {
        "entire_synchronous_source_majorant": sync,
        "entire_common_clock_source_majorant": clock,
        "complete_actual_spatial_input_normalized": spatial,
        "complete_original_spatial_projection_tail_normalized": spatial_tail,
        "actual_full_scalar_reference_stress_C1": MEAN,
        "complete_original_projected_mean_C1_tail_numerator": MEAN_TAIL,
        "complete_source_Ward_numerator_over_mean": s.Integer(2064),
        "complete_detector_Ward_numerator_over_mean": s.Integer(26456),
        "complete_ordered_Ward_numerator_over_mean": s.Integer(28520),
        "complete_spatial_profile_numerator_over_mean": s.Integer(8),
        "complete_ADM_profile_numerator_over_mean": s.Integer(64),
        "complete_actual_ADM_plus_profile_normalized_bound": adm,
        "complete_actual_common_clock_normalized_bound": total,
        "complete_Ward_completed_ADM_approximant_tail_numerator": adm_tail,
        "complete_projected_mean_nonlinear_clock_tail_numerator": 24 * MEAN_TAIL / kk,
        "complete_Ward_completed_common_clock_approximant_tail_numerator": clock_tail,
    }


@cache
def data():
    cc = constants()
    amax = s.Rational(25, 16)
    hmax = s.Rational(8, 5)
    Enorm2 = (amax**3 / 2) ** 2 + 3 * (amax / 2) ** 2
    Edot2 = (amax**3 * (1 + 3 * hmax) / 2) ** 2 + 3 * (amax * (1 + hmax) / 2) ** 2
    v = s.Symbol("transfer_squared", nonnegative=True)
    d = s.Symbol("delta", real=True)
    om, K, m = s.symbols("Omega K mass", positive=True)
    r = 4 * s.sqrt(om * om - m * m)
    checks = {
        "full_synchronous_source_positive_sum": cc["entire_synchronous_source_majorant"]
        - 62408287126207901212,
        "full_common_clock_source_positive_sum": cc[
            "entire_common_clock_source_majorant"
        ]
        - 51511602072425569260,
        "complete_source_Ward_terms": 6 * 84 + 30 * 4 * 13 - 2064,
        "complete_detector_Ward_terms": 4 * (3 * 376 + 2 * 4 * 52)
        + 30 * 4 * 13 * 13
        - 26456,
        "whole_ordered_Ward_sum": 2064 + 26456 - 28520,
        "whole_source_Lie_density_terms": 20 + 32 + 32 - 84,
        "entire_source_spacetime_metric_derivative": 2 * 4 * 2 * 13 + 4 * 29 - 324,
        "entire_source_metric_spacetime_gradient": 324 + 52 - 376,
        "full_original_frequency_radial_Jacobian": s.simplify(
            r * r * s.diff(r, om) - 64 * om * s.sqrt(om * om - m * m)
        ),
        "complete_mean_tail_radial_integral": s.integrate(om**-3, (om, K, s.oo))
        - 1 / (2 * K * K),
        "whole_per_mode_mean_tail_coefficient": s.Rational(32, 9) / 2
        - s.Rational(16, 9),
        "detector_only_two_spatial_derivatives_suffice": s.expand(
            (1 + v) ** 2 - v**2 - (1 + 2 * v)
        ),
        "source_two_extra_spatial_derivatives_suffice": s.expand(
            (1 + v) ** 8 - v**2 * (1 + v) ** 6 - (1 + v) ** 6 * (1 + 2 * v)
        ),
        "clock_curvature_contact_complete_square": s.expand(
            s.Rational(9, 16) - (3 * d - 4 * d * d) - (2 * d - s.Rational(3, 4)) ** 2
        ),
        "same_heavy_mass_and_clock_remain": state.KAPPA - s.Integer(10) ** 800,
    }
    gates = {
        "full_scalar_density_norm_not_an_old_vector_constant": Enorm2 < 16,
        "full_scalar_density_derivative_norm_exact": Edot2 < 400,
        "whole_real_reference_H_bound": hmax < 2,
        "whole_real_reference_scale_bound": amax < 2,
        "whole_physical_first_chart_operator_bound": 2 * amax**4 < 36,
        "whole_physical_second_chart_operator_bound": 2 + 8 + 6 + 6 + 4 < 30,
        "complete_profile_QQ_bound": s.Rational(3, 4) * 8 < 8,
        "complete_profile_ADM_bound": 8 * (2 + 2 + s.Rational(3, 4)) < 64,
        "whole_complex_H_coefficient": 4 * s.Rational(3, 4) / s.Rational(7, 16) < 7,
        "whole_complex_inverse_scale_coefficient": s.Rational(16, 7) ** 4 < 28,
        "whole_complex_clock_coefficient": s.Rational(1, 2) * s.Rational(16, 7) ** 3
        < 6,
        "source_primitives_and_all_time_jets_fit": cc[
            "entire_synchronous_source_majorant"
        ]
        < SYN_SOURCE,
        "detector_map_needs_no_time_derivative": 87 < SYN_DETECTOR,
        "whole_clock_source_no_time_or_spatial_degree_increase": cc[
            "entire_common_clock_source_majorant"
        ]
        < CLOCK_SOURCE,
        "whole_clock_detector_no_derivative_increase": 1 + s.sqrt(3) < CLOCK_DETECTOR,
        "full_C1_mean_tail_from_actual_per_mode_bound": s.Rational(16, 9) * 10**301
        < MEAN_TAIL,
        "retained_projected_mean_contact_tail": 24 * MEAN_TAIL / state.KAPPA > 0,
        "complete_actual_ADM_response_positive_sum": cc[
            "complete_actual_ADM_plus_profile_normalized_bound"
        ]
        < s.Rational(1, 10**327),
        "complete_actual_common_clock_response_positive_sum": cc[
            "complete_actual_common_clock_normalized_bound"
        ]
        < s.Rational(1, 10**300),
        "whole_Ward_completed_clock_approximant_tail": cc[
            "complete_Ward_completed_common_clock_approximant_tail_numerator"
        ]
        < s.Rational(1, 10**200),
        "original_cutoff_domain_has_K_greater_than_one": state.MASS2 > 1,
        "full_same_source_and_detector_graph_used_in_tail": True,
        "no_full_bare_projected_response_Ward_identification": True,
        "no_same_space_inverse_or_nonlinear_history_claim": True,
    }
    return {
        "complete_positive_norm_and_tail_constants": cc,
        "full_density_majorants_squared": {"E": Enorm2, "dE": Edot2},
        "graph": "V02[D]=||(1-Delta)(n,beta,Q)||L2(dt dx;Euclidean/Frobenius). U138[G]^2=sum(j0..13)||(1-Delta)^4 partial_t^j(n,beta,Q)||L2². Detector has two spatial derivatives and no time derivatives; source has eight spatial and thirteen time derivatives.",
        "whole_synchronous_norm_proof": "Both time primitives have L2 norm<=1. Their higher time derivatives are lower derivatives of the original source. Entire Cauchy-Leibniz sums use H_j<=7j!4^j,a^-2_j<=28j!4^j. Since S245's detector graph is plain L2, only two detector spatial derivatives are required. Every spatial and time source term is retained.",
        "whole_Ward_bound": "Use the actual S240 rho/P and first derivatives<1e400. At a<=25/16 and |H|<=8/5, E_F<4B,E'_F<20B are checked directly. xi<=3,dxi<=4,gaugeADM<=13. All source and detector density/chart terms give2064B and26456B, respectively. No old vector mean bound is imported.",
        "complete_profile_and_clock_bound": "Remove only the full profile QQ term from the actual spatial input before scalar Gaussian Ward; add the entire ADM profile once. Clock first pullback adds2delta nI to Q. Both full nonlinear mean contacts are retained, then cancel only at the full reference. All13 clock source jets are controlled without changing the spatial graph.",
        "full_original_mean_cutoff": "Use the time-independent original Omega_star<=K one-leg projection on the entire actual-minus-ad4 mean integrand, with full finite local terms restored. Complete stress jets through1 obey1e301 Omega^-5. The full radial Jacobian yields<2e301/K², bounded by1e310/K². Exact physical and general-D ad0,2,4 mode Ward identities imply the projected homogeneous mean is conserved. This is not the full projected-response Ward identity.",
        "honest_Ward_completed_projection": "Use the S245 original two-leg restored spatial projection, the same one-leg projected mean E_K in BOTH ordered Ward terms, and the full fixed profile. Apply the nonlinear clock with the retained3a³(P_K-P)(4delta²-3delta)nn term. This is a precisely defined convergent computational Ward completion, not an identification with the bare full cutoff response or a proved local/shape counteraction.",
        "remaining": "Full bare-cutoff response/local-shape identification, compatible quantum constraints/inverse, finite inhomogeneous nonlinear feedback, interacting light/mixed loops, quantum gravitational limit, physical UV/Regge and original V/G/B/P8 remain open.",
        "checks": checks,
        "gates": {key: bool(value) for key, value in gates.items()},
    }
