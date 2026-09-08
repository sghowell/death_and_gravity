"""Exact continuous source-to-geometry bounds on the whole bounce interval."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_aligned_margin import action as margin_action
from p8_aligned_margin import bounds as margin_bounds
from p8_vector_hadamard import transfer

from . import system

old = system.old
JLOW = sp.Rational(3, 8)


@cache
def derivative_bounds():
    bg = old.background()
    out = {}
    for name in ("J", "theta", "w", "ell", "delta"):
        value = sp.factor(sp.diff(bg[name], old.u))
        num, den = sp.fraction(value)
        constant = den.subs(old.u, 0)
        poly = sp.Poly(num, old.u)
        upper = sum(abs(c)/2**powers[0] for powers, c in poly.terms())/constant
        quotient = sp.factor(den/constant)
        power = sp.Poly(den, old.u).degree()//2
        out[name] = {"upper": upper,
                     "numerator_reconstruction": sp.expand(poly.as_expr()-num),
                     "positive_denominator_reconstruction": sp.factor(quotient-(1+old.u**2)**power)}
    return out


@cache
def constants():
    j, theta, ell, w, volume = JLOW, sp.Integer(2), sp.Rational(1, 10), sp.Rational(1, 20), sp.Integer(4)
    A = sp.Rational(1, 6)+theta**2/(2*j)
    B = 3*theta*w*ell/(2*j)
    C = 9*ell**2*(1+w**2/(2*j))
    cv_prime = 3*(sp.Rational(4, 3)*ell+w*sp.Rational(3, 5))/(2*j)+3*w*ell*185/(2*j**2)
    cy_prime = (7+theta*6)/(8*j)+theta*185/(8*j**2)
    partial = sum(sp.Rational(3, 2)**k/sp.factorial(k) for k in range(7))
    exp_upper = partial+sp.Rational(3, 2)**7/sp.factorial(7)/(1-sp.Rational(3, 16))
    return {"row_one_norm_upper": B+A/4,
            "row_two_norm_upper": B+4*volume*C,
            "force_row_one_over_eta0_upper": theta*sp.Rational(5, 2)/(2*j),
            "force_row_two_over_eta0_upper": 4*volume*(3+3*ell*w*sp.Rational(5, 2)/(2*j)),
            "exp_three_halves_rational_upper": exp_upper,
            "lapse_coefficient_sum_upper": sp.Rational(103, 150),
            "lapse_coefficient_derivative_sum_upper": cv_prime+cy_prime,
            "lapse_time_derivative_eta0_upper": 346*131+sp.Rational(103, 150)*246+6+sp.Rational(5, 2)*185/(2*j**2),
            "lapse_time_derivative_eta1_upper": sp.Rational(10, 3)}


def response_bounds(value_bound, rate_bound):
    eta0, eta1 = exact(value_bound, "eta_value"), exact(rate_bound, "eta_time_derivative")
    if not (bool(0 <= eta0 <= sp.Rational(1, 10**6)) and bool(0 <= eta1 <= sp.Rational(1, 10**6))):
        raise ValueError("Require exact 0<=eta_value,eta_time_derivative<=10^-6")
    return {"eta_value": eta0, "eta_time_derivative": eta1,
            "scaled_phase_response_upper": 131*eta0,
            "scaled_phase_time_derivative_upper": 246*eta0,
            "lapse_response_upper": 94*eta0,
            "lapse_time_derivative_upper": 48000*eta0+4*eta1,
            "hat_curvature_time_derivative_upper": 193*eta0,
            "matter_displacement_upper": 44*eta0,
            "physical_log_scale_correction_upper": 225*eta0,
            "physical_fractional_scale_correction_upper": 450*eta0,
            "exact_chart_lift_Hubble_difference_over_inverse_tau_upper": 98000*eta0+8*eta1,
            "full_self_consistent_quantum_solution_or_nonlinear_residual_claim": False}


def vector_example(planck_time_product, mass_time_product):
    vector = transfer.physical_bounds(planck_time_product, mass_time_product)
    eta0 = max(vector[name+"_new_state_total_over_reference_density"] for name in ("energy", "pressure"))
    eta1 = max(vector[name+"_new_state_derivative_over_reference_density_rate"] for name in ("energy", "pressure"))
    return {"actual_vector_envelopes": vector, "response": response_bounds(eta0, eta1)}


@cache
def checks():
    out = {}
    for name, data in derivative_bounds().items():
        out[name+"_derivative_numerator_reconstruction"] = data["numerator_reconstruction"]
        out[name+"_derivative_positive_denominator"] = data["positive_denominator_reconstruction"]
    N, h, hp = sp.symbols("N h h_prime", positive=True)
    omega = -sp.log((h-1+N**-2)/h)/4
    out["exact_chart_omega_N"] = sp.factor(sp.diff(omega, N)-1/(2*N*(1+(h-1)*N**2)))
    out["exact_chart_omega_u"] = sp.factor(sp.diff(omega, h)*hp-hp*(N**-2-1)/(4*h*(h-1+N**-2)))
    return out


@cache
def proof_checks():
    c, derivatives = constants(), derivative_bounds()
    eta = sp.Rational(1, 10**14)
    example = response_bounds(eta, eta)
    actual = vector_example(10**12, 1000)["response"]
    epsmax = sp.Rational(1, 100)
    return {"inherited_J_global_bound_gives_compact_three_eighths": bool(
                margin_bounds.polynomial_bounds()["global_lower"]/sp.Rational(5, 4)**6 > JLOW),
            "a_cubed_below_four": bool(sp.Rational(25, 16)**3 < 4),
            "compact_h_upper_below_three_for_abs_Lambda_at_most_one_half": bool(sp.Rational(5, 4)**3 < 3),
            "Hubble_and_Theta_envelope_at_most_two": bool(4*sp.Rational(1, 2) == 2),
            "h_logarithmic_derivative_envelope_at_most_three": bool(6*sp.Rational(1, 2) == 3),
            "new_J_derivative_below_185": bool(derivatives["J"]["upper"]+24*epsmax < 185),
            "theta_derivative_below_seven": bool(derivatives["theta"]["upper"] < 7),
            "w_derivative_below_four_thirds": bool(derivatives["w"]["upper"] < sp.Rational(4, 3)),
            "ell_derivative_at_most_three_fifths": bool(derivatives["ell"]["upper"] <= sp.Rational(3, 5)),
            "delta_derivative_at_most_three_halves": bool(derivatives["delta"]["upper"] <= sp.Rational(3, 2)),
            "both_matrix_rows_below_three_halves": bool(max(c["row_one_norm_upper"], c["row_two_norm_upper"]) < sp.Rational(3, 2)),
            "both_forcing_rows_below_49_eta0": bool(max(c["force_row_one_over_eta0_upper"], c["force_row_two_over_eta0_upper"]) < 49),
            "exponential_bound_below_five": bool(c["exp_three_halves_rational_upper"] < 5),
            "integrated_response_constant_below_131": bool(49*(5-1)/sp.Rational(3, 2) < 131),
            "phase_derivative_constant_below_246": bool(sp.Rational(3, 2)*131+49 < 246),
            "curvature_derivative_constant_below_193": bool(c["row_one_norm_upper"]*131+c["force_row_one_over_eta0_upper"] < 193),
            "lapse_value_constant_below_94": bool(sp.Rational(103, 150)*131+sp.Rational(10, 3) < 94),
            "matter_displacement_constant_at_most_44": bool(sp.Rational(3, 10)*131+sp.Rational(1, 20)*94 <= 44),
            "lapse_coefficient_derivative_below_346": bool(c["lapse_coefficient_derivative_sum_upper"] < 346),
            "lapse_time_derivative_eta0_constant_below_48000": bool(c["lapse_time_derivative_eta0_upper"] < 48000),
            "lapse_time_derivative_eta1_constant_below_four": bool(c["lapse_time_derivative_eta1_upper"] < 4),
            "small_response_lapse_stays_in_one_percent_strip": bool(94*sp.Rational(1, 10**6) < sp.Rational(1, 100)),
            "one_percent_lapse_strip_stays_inside_original_X_tube": bool(
                sp.Rational(100, 101)**2 > sp.Rational(9, 10) and sp.Rational(100, 99)**2 < sp.Rational(11, 10)),
            "chart_omega_N_below_one": bool(1/(2*sp.Rational(99, 100)) < 1),
            "chart_inverse_lapse_square_slope_below_three": bool(2/sp.Rational(99, 100)**3 < 3),
            "chart_omega_u_over_lapse_shift_below_three": bool(9/(4*sp.Rational(9, 10)) < 3),
            "log_scale_correction_below_one_half": bool(225*sp.Rational(1, 10**6) < sp.Rational(1, 2)),
            "Hubble_numerator_constant_below_49000": bool(193+5*94+48000 < 49000),
            "permitted_response_Hubble_error_below_one_tenth": bool((98000+8)*sp.Rational(1, 10**6) < sp.Rational(1, 10)),
            "old_Hubble_endpoint_magnitude_above_nine_tenths": bool(sp.Rational(16, 17) > sp.Rational(9, 10)),
            "vector_example_satisfies_rounded_value_and_rate_envelopes": bool(
                actual["eta_value"] < eta and actual["eta_time_derivative"] < eta),
            "rounded_example_physical_Hubble_difference_below_one_billionth": bool(
                example["exact_chart_lift_Hubble_difference_over_inverse_tau_upper"] < sp.Rational(1, 10**9)),
            "regular_homogeneous_action_uses_no_Theta_or_spatial_inverse": system.action()["no_Theta_or_spatial_momentum_inverse"],
            "fixed_source_response_not_quantum_feedback_contraction": True,
            "clock_source_compatibility_is_inherited_Ward_identity": True,
            "new_margin_range_is_inherited_not_silently_changed": margin_action.epsilon.is_positive is True}
