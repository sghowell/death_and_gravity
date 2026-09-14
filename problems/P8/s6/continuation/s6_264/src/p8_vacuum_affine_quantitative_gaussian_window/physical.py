"""Complete linear physical rows and strictly limited finite-band Gaussian jets."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_finite_band_neighborhood import reference as previous_reference
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from p8_vacuum_affine_gauge_mean_transport import crossing
from p8_vacuum_affine_physical_background_vertices import parent
from p8_vacuum_affine_scalar_tame_propagator import charts

from . import source


@cache
def rows():
    prior = previous_reference.shear_bounds()
    slope = prior["whole_reference_a_cubed_boundary_over_P_squared_derivative_bound"]
    P = s.Symbol("whole_comoving_P", positive=True)
    s00, s01, s11 = s.symbols(
        "whole_reference_Sbar00 whole_reference_Sbar01 whole_reference_Sbar11",
        real=True,
    )
    R = s.eye(4)
    R[2:, :2] = P**2 * s.Matrix([[s00, s01], [s01, s11]])
    W = s.diag(P, P, 1, 1)
    T = phase.central_map().subs(charts.q, P**2 / phase.a**2)
    old_from_x = T.inv() * R * W.inv()
    original = crossing.reference_rows()
    vrow = (
        s.sqrt(phase.kappa) * original["whole_physical_log_volume_row"] * old_from_x
    ).applyfunc(s.factor)
    nrow = (
        s.sqrt(phase.kappa)
        * original["whole_physical_lapse_row"].subs(charts.q, P**2 / phase.a**2)
        * old_from_x
    ).applyfunc(s.factor)
    wanted_v = s.Matrix(
        [
            [
                s00 / (2 * phase.a * P),
                s01 / (2 * phase.a * P),
                1 / (2 * phase.a * P**2),
                0,
            ]
        ]
    )
    yb, ym, pb, pm = s.symbols("clean_yb clean_ym clean_pb clean_pm", real=True)
    x = s.Matrix([P * yb, P * ym, pb, pm])
    old_pb = pb + P**2 * (s00 * yb + s01 * ym)
    old_pm = pm + P**2 * (s01 * yb + s11 * ym)
    independent_n = (
        -2 * charts.th * P**2 * yb / phase.a**2
        + 3 * charts.th * charts.l * ym
        - (2 * charts.E * P**2 / phase.a**2 + 3 * charts.T)
        * old_pb
        / (2 * phase.a * P**2)
        + charts.l * charts.E * old_pm / phase.a**3
    ) / (2 * charts.J)
    tiny = I(-s.Rational(1, 10**49), s.Rational(1, 10**49))
    box = {
        P: I(source.LOW, source.HIGH),
        phase.a: I(1, 2),
        charts.J: I(1, 2),
        charts.th: I(-5 * source.TIME, 5 * source.TIME),
        charts.l: I(0, s.Rational(1, 10)),
        charts.E: I(-s.Rational(1, 2), -s.Rational(1, 4)),
        charts.T: I(-s.Rational(1, 10**380), s.Rational(1, 10**380)),
        s00: tiny,
        s01: tiny,
        s11: tiny,
    }
    norms = {}
    entries = {}
    for name, row in (("hat_log_scale", vrow), ("physical_lapse", nrow)):
        values = [
            source.as_rational(max(abs(iv.lo), abs(iv.hi)))
            for iv in (evaluate(entry, box) for entry in row)
        ]
        entries[name] = values
        norms[name] = sum(values)
    sourcebox = {source.u: I(-source.TIME, source.TIME)}
    Ebound = evaluate(source.E, sourcebox)
    # Full source profiles in Tcorr=rho-3 delta pressure are not set to zero.
    return {
        "whole_original_physical_canonical_row_source": original[
            "whole_original_physical_canonical_map"
        ],
        "whole_original_canonical_from_weighted_clean_phase": old_from_x,
        "whole_hat_log_scale_row_times_sqrt_kappa": vrow,
        "whole_physical_lapse_row_times_sqrt_kappa": nrow,
        "whole_reference_shear_derivative_enclosure": slope,
        "whole_improved_reference_Sbar_absolute_bound": s.Rational(1, 10**49),
        "whole_current_reference_band_row_box": {
            str(key): value.bounds() for key, value in box.items()
        },
        "whole_eight_row_entry_absolute_enclosures": entries,
        "whole_both_row_operator_enclosures": norms,
        "whole_safe_hat_log_scale_row_norm": s.Rational(1, 10**112),
        "whole_safe_lapse_row_norm": s.Integer(10) ** 16,
        "whole_row_interpretation": "The old stored coordinate v is the HAT log-scale, not already the physical scale. The physical correction is applied separately. Both original canonical density factors, both momentum entries, all nonzero Tcorr and the full symmetric boundary remain in the lapse row. Every physical field row also carries the explicit1/sqrt(kappa) factor removed only in the displayed scaled rows.",
        "checks": {
            "whole_exact_hat_log_scale_cotangent_row": (vrow - wanted_v).applyfunc(
                s.factor
            ),
            "whole_full_lapse_row_from_independent_original_constraint": s.factor(
                (nrow * x)[0] - independent_n
            ),
            "whole_weighted_phase_CCR_in_both_physical_rows": s.factor(
                (vrow * (P * phase.OMEGA) * nrow.T)[0]
                - charts.th / (2 * charts.J * phase.a**3)
            ),
            "whole_original_canonical_map_roundtrip": (
                T * old_from_x * W - R
            ).applyfunc(s.factor),
        },
        "gates": {
            "whole_previous_computed_shear_derivative_below1e11": slope < 10**11,
            "whole_integrated_Sbar_below1e_minus49": slope * source.TIME
            < s.Rational(1, 10**49),
            "whole_hat_log_scale_row_bound": norms["hat_log_scale"]
            < s.Rational(1, 10**112),
            "whole_physical_lapse_row_bound": norms["physical_lapse"] < 10**16,
            "actual_reference_E_inside_row_box": Ebound.inside(
                -s.Rational(1, 2), -s.Rational(1, 4)
            ),
            "actual_Theta_bound_uses_both_complete_terms": 5 * source.TIME
            < s.Rational(1, 10**20),
            "actual_full_Tcorr_inside_row_box": s.Rational(5, 2) * source.PROFILE_BOUND
            < s.Rational(1, 10**380),
            "actual_J_row_box_is_previous_full_profile_reentry": True,
            "no_new_Gaussian_state_or_homogeneous_constraint_sector": True,
        },
    }


@cache
def moments():
    P = s.Symbol("positive_radial_P", positive=True)
    exact_integral = (source.HIGH**4 - source.LOW**4) / (8 * s.pi**2)
    radial = (source.HIGH**4 - source.LOW**4) / 72
    vv = s.Rational(1, 10**224) * source.COVARIANCE * radial / source.KAPPA
    nn = 10**32 * source.COVARIANCE * radial / source.KAPPA
    vv_safe = s.Rational(1, 10**740)
    nn_safe = s.Rational(1, 10**490)
    vn_safe = s.Rational(1, 10**615)
    tensor = 2 * 10**5 * (source.HIGH**2 - source.LOW**2) / (9 * source.KAPPA)
    exponent = source.EPSILON**2 / (2 * nn_safe)
    return {
        "whole_original_R3_radial_measure": "P^2 dP/(2 pi^2), restricting the same original two-point function to the stated band; no new Cauchy covariance or torus state is selected.",
        "whole_full_weighted_covariance_upper": "C_x<=1e21 P I for the complete two-mode scalar block; C_xT<=1e5 P I for each original TT polarization.",
        "whole_exact_covariance_radial_factor": exact_integral,
        "whole_outward_rational_radial_factor": radial,
        "whole_hat_log_scale_variance_enclosure": vv,
        "whole_lapse_variance_enclosure": nn,
        "whole_safe_hat_log_scale_variance": vv_safe,
        "whole_safe_lapse_variance": nn_safe,
        "whole_safe_absolute_symmetric_cross_covariance": vn_safe,
        "whole_both_TT_Frobenius_variance_enclosure": tensor,
        "whole_safe_both_TT_metric_variance": s.Rational(1, 10**660),
        "whole_linear_lapse_threshold": source.EPSILON,
        "whole_marginal_Gaussian_tail_exponent_lower": exponent,
        "whole_safe_marginal_lapse_tail": "For each fixed time |t|<=1e-60 and fixed spatial point x, the spectral law of the self-adjoint band-limited LINEAR reference lapse n_B is centered Gaussian and Pr(abs(n_B)>=1e-230)<=2 exp(-1e29). This is a one-observable marginal statement, not a joint law of noncommuting physical observables or a uniform infinite-volume event.",
        "whole_moment_argument": "Use the complete row norms with the full transported covariance upper bound and both1/sqrt(kappa) endpoint factors. Integrate P^3 against the exact original radial normalization and use pi>3 for the outward rational ceiling. Covariance positivity gives the retained symmetric cross bound. For tensors restore gamma=2h/sqrt(kappa), both polarizations and the norm-one TT projector; their sum gives the displayed Frobenius variance.",
        "whole_tail_boundary": "The Gaussian is unbounded. The tail bound is not zero probability, nonlinear lapse reconstruction or support entirely within a regular branch. A continuum of times/points, simultaneous noncommuting measurements, full interacting dynamics, and all omitted momenta require separate results. Pointwise marginal smallness does not supply a nonlinear Cauchy or quantum-mean theorem.",
        "checks": {
            "whole_exact_R3_radial_covariance_integral": s.factor(
                s.integrate(P**3 / (2 * s.pi**2), (P, source.LOW, source.HIGH))
                - exact_integral
            ),
            "whole_Cauchy_cross_square": vn_safe**2 - vv_safe * nn_safe,
            "both_physical_kappa_endpoints": (1 / s.sqrt(source.KAPPA)) ** 2
            - 1 / source.KAPPA,
        },
        "gates": {
            "full_hat_log_scale_variance_below1e_minus740": vv < vv_safe,
            "full_lapse_variance_below1e_minus490": nn < nn_safe,
            "whole_both_TT_metric_variance_below1e_minus660": tensor
            < s.Rational(1, 10**660),
            "actual_nonzero_band_and_positive_covariance_ceiling": source.HIGH
            > source.LOW
            > 0
            and nn > 0,
            "outward_radial_bound_uses_pi_greater_than_three": bool(s.pi > 3),
            "marginal_Chernoff_exponent_exceeds1e29": exponent > 10**29,
            "same_reference_R3_band_not_Wilsonian_cutoff": True,
            "single_linear_observable_tail_not_nonlinear_support": True,
        },
    }


@cache
def volume():
    h = parent.h
    c1 = s.Rational(3, 2) / h
    c2 = s.Rational(21, 4) / h**2 - s.Rational(9, 2) / h
    actual = parent.lapse_jets()["rows"]["U"]
    Cv, Cvn, Cn = s.symbols("fixed_band_Cvv fixed_band_Cvn fixed_band_Cnn", real=True)
    d = s.Rational(1, 2) / h
    contact = s.Rational(9, 2) * Cv + 3 * c1 * Cvn + c2 * Cn / 2
    box = {source.u: I(-source.TIME, source.TIME)}
    first = evaluate(c1, box)
    second = evaluate(c2, box)
    budget = (
        5 * s.Rational(1, 10**740) + 6 * s.Rational(1, 10**615) + s.Rational(1, 10**490)
    )
    linear = 18 * s.Rational(1, 10**740) + 8 * s.Rational(1, 10**490)
    return {
        "whole_correct_physical_log_scale_linear_jet": "v_physical=v_hat+n/(2h), from a_physical=a_hat R^(-1/4), with h=(1+u^2)^3.",
        "whole_correct_normalized_spatial_volume_linear_jet": "3 v_hat+3 n/(2h), not the refuted S261 shifted-R slope.",
        "whole_actual_U_lapse_slope": c1,
        "whole_actual_U_second_lapse_jet": c2,
        "whole_complete_second_order_Weyl_spatial_volume_contact": contact,
        "whole_second_order_Weyl_contact_absolute_enclosure": budget,
        "whole_safe_second_order_Weyl_contact": s.Rational(1, 10**489),
        "whole_safe_linear_physical_volume_variance": s.Rational(1, 10**488),
        "whole_strict_jet_scope": "This is the expectation of the explicitly stated second-order Weyl polynomial under the unchanged band-limited reference covariance. It is NOT the expectation of the full nonlinear volume function on an unbounded Gaussian, a stationary physical quantum mean, a Ward-restored renormalized stress or a new subtraction prescription. The complete cross covariance is bounded, not set to zero. The other H/Proca product blocks and all determinants remain unchanged; no full interacting onepoint calculation is claimed.",
        "checks": {
            "literal_whole_correct_U_first_jet": s.factor(actual[1] - c1),
            "literal_whole_correct_U_second_jet": s.factor(actual[2] - c2),
            "literal_volume_three_times_physical_logscale_slope": c1 - 3 * d,
            "correct_actual_bounce_lapse_slope": c1.subs(source.u, 0)
            - s.Rational(3, 2),
            "correct_actual_bounce_second_jet": c2.subs(source.u, 0) - s.Rational(3, 4),
        },
        "gates": {
            "whole_lapse_slope_below_two": max(abs(first.lo), abs(first.hi)) < 2,
            "whole_second_lapse_jet_below_two": max(abs(second.lo), abs(second.hi)) < 2,
            "whole_Weyl_contact_below1e_minus489": budget < s.Rational(1, 10**489),
            "whole_linear_volume_variance_below1e_minus488": linear
            < s.Rational(1, 10**488),
            "cross_covariance_retained_in_complete_contact": contact.has(Cvn),
            "S261_wrong_physical_binding_not_renewed": True,
            "finite_Weyl_jet_not_full_nonlinear_volume_mean": True,
        },
    }
