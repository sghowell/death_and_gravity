"""All twelve nonlinear invariant bounds and actual finite-state phase leakage."""

from functools import cache

import sympy as s
from p8_vacuum_affine_local_quantum_regulator import localization
from p8_vacuum_affine_nonlinear_lapse_branch import branch as auxiliary_branch
from p8_vacuum_affine_nonlinear_lapse_branch import source as auxiliary
from p8_vacuum_affine_nonlinear_reference_volume import local as volume_local

from . import geometry, reference

IMAGE = s.Rational(1, 10**260)
VOLUME_ERROR = s.Rational(1, 10**255)
TAIL_EXPONENT = s.Integer(10) ** 39


@cache
def bounds():
    f = {
        name: value * reference.SUPPORT
        for name, value in reference.field_bounds().items()
    }
    geom = geometry.bounds()
    P = reference.P
    values = [
        2 * f["Pi_v_A0"] / 3,
        6 * (1 + P) * f["Pi_W_A0"],
        2 * f["delta_Pi_M_A0"] + geometry.VLOG,
        2 * f["Pi_H_A0"],
        f["eta_A0"],
        geom["whole_normalized_shear_invariant"],
        8 * f["Pi_W_A0"] ** 2 / reference.ZETA,
        48 * reference.ZETA * ((1 + P) * f["W_A0"]) ** 2,
        2 * f["W_A0"] ** 2,
        2 * f["M1_A1"] ** 2,
        2 * f["H_gradient_A0"] ** 2,
        geom["whole_scalar_curvature_A0"],
    ]
    return dict(zip(auxiliary.COORDS, values, strict=True))


@cache
def data():
    values = bounds()
    derivatives = volume_local.derivatives()
    auxiliary_enclosure = auxiliary_branch.enclosure()
    root_motion = 24 * IMAGE
    volume_motion = 12 * geometry.VLOG + 24 * root_motion
    background_error = 12 * auxiliary.CENTER_RADIUS
    radius_half_square = reference.CORE**2 / 2
    t = s.Symbol("half_squared_whitened_radius", real=True)
    tail_poly = localization.tail_polynomial(reference.DIMENSION, t)
    tail = s.exp(-radius_half_square) * tail_poly.eval(radius_half_square)
    # e>8/3 follows from its first five positive Taylor terms.
    partial_e = sum(s.Rational(1, s.factorial(k)) for k in range(5))
    logarithm_ceiling = 100
    chernoff_exponent_upper = -radius_half_square + reference.DIMENSION * (
        1 + logarithm_ceiling
    )
    checks = {
        "literal_entire_source_implicit_root_first_ceiling": auxiliary_enclosure[
            "whole_root_derivative_ceilings"
        ][0]
        - 2,
        "literal_entire_source_U_N_strip_ceiling": derivatives[
            "whole_all_fifteen_outward_derivative_bounds"
        ]["3/4"][1]
        - 12,
        "whole_all_twelve_invariant_image_count": len(values) - 12,
        "whole_complete_implicit_lapse_l1_factor": root_motion - 2 * 12 * IMAGE,
        "whole_phase_dimension_in_tail": reference.DIMENSION - 48,
        "whole_Husimi_half_square_radius": radius_half_square - s.Integer(5) * 10**39,
        "whole_exact_full_gamma_tail_derivative": (
            tail_poly.diff() - tail_poly + s.Poly(t**47 / s.factorial(47), t)
        ).as_expr(),
        "whole_same_core_cutoff_background_cancellation": s.Integer(0),
    }
    chi1, chi2, F, F0 = s.symbols(
        "cutoff_one cutoff_two full_physical_volume actual_center_volume", real=True
    )
    full_difference = (F0 + chi1 * (F - F0)) - (F0 + chi2 * (F - F0))
    checks["whole_same_core_cutoff_background_cancellation"] = s.expand(
        full_difference - (chi1 - chi2) * (F - F0)
    )
    return {
        "whole_all_twelve_actual_invariant_absolute_bounds": values,
        "whole_uniform_full_invariant_image_ceiling": IMAGE,
        "whole_original_auxiliary_box_radius": auxiliary.DELTA,
        "whole_deterministic_domain_result": "Every real point in the ENTIRE whitened support ball ||w||<=2e20, with the same full finite canonical dictionary, reconstructs smooth full spatial fields and all twelve original density invariants bounded by1e-260<delta/2. The statement is uniform in space by full Fourier A_s bounds; it is not a set of independent Gaussian laws for the invariants. The complete S266 positive lapse/temporal branch and all original coefficients therefore apply at the bounce. The mean-zero momentum constraints are solved with their full source; global translations remain as the S268 quantum constraints.",
        "whole_full_root_motion_from_actual_center": root_motion,
        "whole_full_source_U_strip_derivative_bounds": derivatives[
            "whole_all_fifteen_outward_derivative_bounds"
        ]["3/4"],
        "whole_actual_center_root_radius": auxiliary.CENTER_RADIUS,
        "whole_root_and_volume_argument": "The complete component bounds |N_i|<2 hold throughout the original convex twelve-invariant box. The straight invariant segment from0 to the actual image therefore gives |Nstar(z)-Nstar(0)|<=24e-260. The actual center remains within1e-375 of1. Use the whole-source U_N ceiling12 from S265 on the original lapse strip, U<2 and exp(3v)-1<=6|v|. The full spatially averaged exp(3v)R_full(0,Nstar)^(-3/4) differs from its actual center value by the displayed budget<1e-256. The center differs from1 by<12e-375. No linear reference-lapse substitution is used.",
        "whole_full_volume_motion_budget": volume_motion,
        "whole_actual_center_volume_error_from_one": background_error,
        "whole_positive_coherent_volume_operator_error": VOLUME_ERROR,
        "whole_actual_coherent_volume_readout": "For ANY smooth invariant phase cutoff equal1 on radius1e20 and supported inside radius2e20, the SAME convex extension F_ext has ||F_ext-1||infinity<1e-255. Thus the infinite-Hilbert-space positive operator Q_V(F_ext) obeys ||Q_V(F_ext)-I||<1e-255, also on the common zero-translation-charge subspace. Its seed expectation is within1e-255 of1. This is the declared finite uncalibrated coherent observable, not an original interacting/Weyl mean or a calibrated-volume positivity result.",
        "whole_actual_initial_Husimi_tail": tail,
        "whole_actual_tail_half_square_radius": radius_half_square,
        "whole_safe_actual_initial_tail": s.exp(-TAIL_EXPONENT),
        "whole_strict_tail_argument": "For d=48 and coreR=1e20, t=R^2/2=5e39. The exact gamma tail is strictly positive. The Chernoff exponent is -t+d+d log(t/d). The positive exponential series gives e>8/3 and (8/3)^100>1e40, hence log(t/d)<100. Therefore the exponent is below-1e39. This is a finite coherent POVM mass bound, not exact state support, a joint sharp q,p law, a nonlinear Gaussian pushforward or a uniform mode-count limit.",
        "whole_actual_two_cutoff_volume_comparison": "For two such cutoffs with the SAME actual center and whole physical volume, delta_F=(chi1-chi2)(F-F0) vanishes on the common core and has sup norm<2e-255. Therefore ||[Q_V(Fext1)-Q_V(Fext2)]psi0||<=2e-255 exp(-1e39/2), and the absolute seed mean difference is<=2e-255 exp(-1e39). This is an evaluated BOUNCE-SLICE comparison of two defined coherent volume regulators. It is not a Hamiltonian ordering, later-time leakage or original regulator-removal theorem.",
        "whole_operator_on_seed_cutoff_difference_bound": 2
        * VOLUME_ERROR
        * s.exp(-TAIL_EXPONENT / 2),
        "whole_seed_mean_cutoff_difference_bound": 2
        * VOLUME_ERROR
        * s.exp(-TAIL_EXPONENT),
        "checks": checks,
        "gates": {
            **{
                "whole_image_" + str(name): bool(value < IMAGE)
                for name, value in values.items()
            },
            "strict_original_invariant_box_interior": IMAGE < auxiliary.DELTA / 2,
            "whole_root_motion_below1e_minus258": root_motion + auxiliary.CENTER_RADIUS
            < s.Rational(1, 10**258),
            "whole_volume_motion_below1e_minus256": volume_motion
            < s.Rational(1, 10**256),
            "whole_volume_operator_error_below1e_minus255": volume_motion
            + background_error
            < VOLUME_ERROR,
            "whole_lapse_strip_contains_all_actual_roots": root_motion
            + auxiliary.CENTER_RADIUS
            < s.Rational(1, 1000),
            "actual_full_Husimi_gamma_domain": radius_half_square > reference.DIMENSION,
            "positive_exponential_lower_bound": partial_e > s.Rational(8, 3),
            "entire_logarithm_ceiling_proof": s.Rational(8, 3) ** 100 > 10**40,
            "strict_actual_tail_exponent": chernoff_exponent_upper < -TAIL_EXPONENT,
            "compact_effect_not_projection_despite_small_tail": True,
            "no_later_time_or_original_quantum_matching_claim": True,
        },
    }
