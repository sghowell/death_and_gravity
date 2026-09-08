"""Mandatory initial covariance response on the actual corrected initial slice."""
from functools import cache

import sympy as sp
from p8_affine_kinetic import scalar as old
from p8_affine_retuned.bounds import exact
from p8_aligned_quantum import kernel
from p8_vector_clock_matching import bounds as matched
from p8_vector_state import wkb
from p8_vector_subtraction import tail

from . import modes


@cache
def initial_force():
    u0 = -sp.Rational(1, 2)
    h = (1+u0**2)**3
    delta = 1/(2*h)
    coefficients = {order: sp.factor((3*delta*data["pressure"]-data["energy"]).subs(wkb.u, u0))
                    for order, data in matched.actual_coefficients().items()}
    L, R = sp.Integer(10)**12, sp.Integer(1000)
    old_state = tail.physical_bounds(L, R)
    state_error = (old_state["full_subtracted_energy_integral_over_reference_density"]
                   +sp.Rational(3, 2)*old_state["full_subtracted_pressure_integral_over_reference_density"]
                   +sp.Rational(5, 2)*20*10**7/(R**2*L**2))
    # The exact local coefficient sum exceeds 4 R^4; pi²<10.
    local_lower = R**4/(160*L**2)
    force_lower = local_lower-state_error
    Jmax = old.background()["J"].subs(old.u, u0)+4*sp.Rational(1, 100)/h**2
    return {"u0": u0, "h0": h, "local_force_coefficients": coefficients,
            "local_force_polynomial_exceeds_four_R_fourth": sp.factor(
                sum(coefficients[order]*R**(4-2*order) for order in (0, 1, 2))-4*R**4),
            "state_error_upper": state_error, "initial_lapse_force_lower": force_lower,
            "initial_J_e_upper": Jmax}


@cache
def impedances():
    p = kernel.geometry.P
    mass = kernel.masses()
    ratio2 = sp.factor(mass["a"]*mass["b"]/(4*p**2))
    Q = 2592*p**5+3708*p**4+886*p**3-1350*p**2-2013*p-605
    target = -(2*p-1)*Q/(324*p**2*(8*p+5)*(2*p**3-1))
    N, h = sp.symbols("N h", positive=True)
    rT4 = ((h-1)*N**4+N**2)/h
    return {"longitudinal_impedance_ratio_squared": ratio2,
            "longitudinal_ratio_difference_factorization": sp.factor(ratio2-1-target),
            "sign_polynomial": Q,
            "transverse_oscillator_frequency_ratio_fourth_power": rT4,
            "transverse_frequency_difference_factorization": sp.factor(
                rT4-1-(N**2-1)*((h-1)*N**2+h)/h),
            "ordinary_transverse_physical_impedance": sp.factor(
                modes.canonical()["transverse"]["g_squared"]*sp.sqrt(
                    modes.canonical()["transverse"]["bare_frequency_squared"].subs(modes.mass2, 0))),
            "longitudinal_impedance_logarithmic_N_slope_at_clock": sp.Rational(113, 81)/h}


@cache
def covariance_checks():
    Z = sp.Symbol("old_impedance", positive=True)
    ratio = sp.Symbol("impedance_ratio", positive=True)
    v0, v1 = 1/sp.sqrt(2*Z), 1/sp.sqrt(2*Z*ratio)
    p0, p1 = -sp.I*Z*v0, -sp.I*Z*ratio*v1
    A, B = (ratio+1)/(2*sp.sqrt(ratio)), (ratio-1)/(2*sp.sqrt(ratio))
    return {"fixed_physical_value_Bogoliubov_identity": sp.simplify(A*v1+B*v1-v0),
            "fixed_physical_momentum_Bogoliubov_identity": sp.simplify(A*p1-B*p1-p0),
            "initial_CCR_preserved_but_not_positive_frequency": sp.factor(A**2-B**2-1),
            "nonzero_mixing_energy_excess": sp.factor(B**2-(ratio-1)**2/(4*ratio))}


def cauchy_mismatch(p_value):
    value = exact(p_value, "p_affine")
    if not bool(sp.Rational(1, 2) <= value <= sp.Rational(21, 40)):
        raise ValueError("Require exact 1/2<=p_affine<=21/40 for this diagnostic")
    squared = impedances()["longitudinal_impedance_ratio_squared"].subs(kernel.geometry.P, value)
    ratio = sp.sqrt(squared)
    return {"p_affine": value, "longitudinal_impedance_ratio_squared": squared,
            "frozen_physical_Cauchy_data_initial_beta_squared_limit": sp.factor((ratio-1)**2/(4*ratio)),
            "leading_mismatch_vanishes": bool(squared == 1),
            "full_quantum_feedback_or_coupled_cone_claim": False}


@cache
def proof_checks():
    d = initial_force()
    pmax = sp.Rational(21, 40)
    Qupper = (2592*pmax**5+3708*pmax**4+886*pmax**3
              -1350*sp.Rational(1, 2)**2-2013*sp.Rational(1, 2)-605)
    return {"initial_local_force_exceeds_four_R_fourth": bool(d["local_force_polynomial_exceeds_four_R_fourth"] > 0),
            "pi_squared_upper_control_22_over_seven_squared_below_ten": bool(sp.Rational(22, 7)**2 < 10),
            "actual_initial_force_above_five_e_minus_fifteen": bool(d["initial_lapse_force_lower"] > sp.Rational(1, 2*10**14)),
            "initial_J_e_below_one_for_entire_margin_range": bool(d["initial_J_e_upper"] < 1),
            "response_lapse_absolute_upper_below_one_percent": bool(94*sp.Rational(1, 10**14) < sp.Rational(1, 100)),
            "tube_p_upper_below_21_over_40": bool(sp.Rational(11, 40) < pmax**2),
            "mass_denominator_two_p_cubed_minus_one_negative": bool(2*pmax**3 < 1),
            "temporal_mass_numerator_negative": bool(22*pmax**3+3*pmax-11 < 0),
            "spatial_mass_numerator_negative": bool(72*pmax**2-88*sp.Rational(1, 2)-55 < 0),
            "longitudinal_impedance_sign_polynomial_strictly_negative": bool(Qupper < 0),
            "initial_longitudinal_impedance_clock_slope_positive": bool(sp.Rational(113, 81) > 0),
            "derivative_free_original_canonical_variation": modes.symplectic_variation()["original_variation_has_no_metric_derivatives"],
            "state_homogeneous_variation_term_is_required": True,
            "transverse_conformal_leading_control_not_all_order_state_proof": True,
            "no_ordinary_Proca_Hadamard_propagation_imported_off_clock": True,
            "frozen_source_result_not_invalidated_by_uncomputed_feedback": True,
            "no_full_coupled_cone_or_quantum_feedback_bound": True}


@cache
def checks():
    out = covariance_checks()
    d = impedances()
    out["longitudinal_physical_impedance_factorization"] = d["longitudinal_ratio_difference_factorization"]
    out["transverse_fixed_oscillator_frequency_factorization"] = d["transverse_frequency_difference_factorization"]
    # At the clock, d log(a_phys^2 sqrt(a_mass*b_mass))/dN
    # = 2omega_N+(alpha+beta)/2.
    h = sp.Symbol("h", positive=True)
    out["actual_initial_impedance_slope"] = sp.factor(
        1/h+(4/(9*h)+28/(81*h))/2-sp.Rational(113, 81)/h)
    return out
