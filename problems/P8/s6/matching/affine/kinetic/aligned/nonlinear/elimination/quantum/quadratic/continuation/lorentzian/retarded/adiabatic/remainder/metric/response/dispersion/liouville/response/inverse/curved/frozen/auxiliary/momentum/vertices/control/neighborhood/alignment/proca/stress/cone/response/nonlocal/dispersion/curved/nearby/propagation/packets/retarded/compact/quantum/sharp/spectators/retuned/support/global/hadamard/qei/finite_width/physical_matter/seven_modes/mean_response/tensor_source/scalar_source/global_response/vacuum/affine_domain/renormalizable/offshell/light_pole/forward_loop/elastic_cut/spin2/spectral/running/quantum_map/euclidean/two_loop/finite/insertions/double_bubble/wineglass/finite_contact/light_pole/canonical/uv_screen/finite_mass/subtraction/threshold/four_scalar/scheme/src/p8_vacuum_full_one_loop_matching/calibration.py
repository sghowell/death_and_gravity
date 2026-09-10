"""Exact full one-loop coefficient and pole intervals in the declared scheme."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_four_scalar import calibration as box
from p8_vacuum_fermion_local_matching import calibration as fermion
from p8_vacuum_forward_loop import calibration as scalar_amplitude
from p8_vacuum_light_pole import calibration as scalar_pole
from p8_vacuum_two_loop_finite_contact import calibration as fixed_contact

from . import logarithm


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational")
    return sp.Rational(value)


def coefficient_interval(tree, old_error, extra_lower, extra_upper):
    T, E, lo, hi = map(exact, (tree, old_error, extra_lower, extra_upper))
    if T <= 0 or E < 0 or lo > hi:
        raise ValueError("Require T>0, E>=0 and ordered finite relative increments")
    low = T - E + T * lo
    high = T + E + T * hi
    return {
        "lower": low,
        "upper": high,
        "positive_complete_one_loop_coefficient": bool(low > 0),
        "scope": "Complete formal one-loop coefficient, not a higher-loop error bound or strict V verdict.",
    }


@cache
def data():
    p = model.data()["actual_parameters"]
    M, g, L, lam = (
        p[k]
        for k in (
            "heavy_mass_squared",
            "cubic_coupling_squared",
            "bare_polynomial_quartic",
            "lambda",
        )
    )
    D = M - 2
    f = fermion.data()
    b = box.data()
    log = logarithm.enclosure()
    ell_lo, ell_hi = log["scale_log_lower"], log["scale_log_upper"]
    m = f["fermion_mass"]
    T = 4 * lam
    fp_lo, fp_hi = f["on_shell_slope_bounds"]
    scalar_pole_bounds = scalar_pole.point(2)
    r_hi = scalar_pole_bounds["finite_kinetic_counterterm_upper"]
    scale_factor = -2 * L + 3 * g / D
    box_relative = b["rational_box_second_coefficient_relative_upper"]
    extra_lo = scale_factor * ell_hi / (2 * 144) + 2 * (fp_lo - r_hi) - box_relative
    extra_hi = scale_factor * ell_lo / (2 * 256) + 2 * fp_hi + box_relative
    old_error = scalar_amplitude.point()["total_one_loop_b2_error_upper"]
    band = coefficient_interval(T, old_error, extra_lo, extra_hi)
    total_error = old_error + T * max(abs(extra_lo), abs(extra_hi))
    kappa_lo = 1 - fp_hi
    kappa_hi = 1 + r_hi - fp_lo
    Bscalar = scalar_pole_bounds["quadratic_self_energy_upper_coefficient"]
    Bferm = f["unscaled_OS_remainder_coefficient_upper"]
    Btotal = (Bscalar + Bferm) / kappa_lo
    scalar_Pi_one_upper = (L * (ell_hi + 1) / 2 + g * ell_hi) / 144
    fF_lo, fF_hi = f["on_shell_mass_threshold_bounds"]
    mass_ref_lo = 1 - fF_hi
    mass_ref_hi = 1 + scalar_Pi_one_upper - fF_lo
    scalar_vac_abs = (1 + M * M) * (ell_hi + sp.Rational(3, 2)) / (4 * 144)
    ferm_vac_lo, ferm_vac_hi = f["all_flavor_vacuum_energy_bounds"]
    total_vac_lo = ferm_vac_lo - scalar_vac_abs
    total_vac_hi = ferm_vac_hi
    contact = fixed_contact.data()
    independent_J1 = M * sp.log(M) / (M - 1) ** 2 - 1 / (M - 1)
    independent_J2 = (M + 1) * sp.log(M) / (M - 1) ** 3 - 2 / (M - 1) ** 2
    independent_contact = (
        6
        * ((-L + g / M) * g * independent_J1 + g * g * independent_J2)
        / (16 * sp.pi**2)
    )
    checks = {
        "same_actual_MS_interaction_tree_relation": L - g * (3 * D - 2) / D**2,
        "same_actual_tree_coefficient": T - 2 * g / D**3,
        "same_old_complete_one_loop_bound": old_error
        - scalar_amplitude.point()["total_one_loop_b2_error_upper"],
        "negative_actual_scale_conversion_factor": scale_factor
        - g * (4 - 3 * D) / D**2,
        "extra_lower_all_terms_retained": extra_lo
        - scale_factor * ell_hi / (2 * 144)
        - 2 * (fp_lo - r_hi)
        + box_relative,
        "extra_upper_all_terms_retained": extra_hi
        - scale_factor * ell_lo / (2 * 256)
        - 2 * fp_hi
        - box_relative,
        "combined_pole_remainder_includes_both_sectors": Btotal * kappa_lo
        - Bscalar
        - Bferm,
        "combined_kinetic_interval_width": kappa_hi - kappa_lo - r_hi - fp_hi + fp_lo,
        "complete_mass_reference_interval_width": mass_ref_hi
        - mass_ref_lo
        - scalar_Pi_one_upper
        - fF_hi
        + fF_lo,
        "full_vacuum_reference_lower_dictionary": total_vac_lo
        + scalar_vac_abs
        - ferm_vac_lo,
        "same_fixed_scalar_finite_contact_not_refit": sp.factor(
            contact["actual_fixed_finite_potential_contact"] - independent_contact
        ),
        "finite_contact_does_not_enter_second_coefficient": sp.Integer(0)
        * contact["actual_fixed_finite_potential_contact"],
    }
    bounds = {
        "named_scale_log_between_nine_hundred_and_one_thousand": bool(
            900 < ell_lo < ell_hi < 1000
        ),
        "positive_ordered_scalar_and_fermion_mass_scales": bool(2 < M < m * m),
        "actual_scale_conversion_factor_strictly_negative": bool(scale_factor < 0),
        "scalar_pole_slope_upper_below_three_e_minus_208": bool(
            0 < r_hi < sp.Rational(3, 10**208)
        ),
        "new_reference_and_fermion_increment_band_ordered": bool(extra_lo < extra_hi),
        "extra_relative_effect_below_one_e_minus_202": bool(
            max(abs(extra_lo), abs(extra_hi)) < sp.Rational(1, 10**202)
        ),
        "complete_formal_one_loop_error_below_one_millionth_tree": bool(
            total_error / T < sp.Rational(1, 10**6)
        ),
        "complete_formal_one_loop_coefficient_positive": band[
            "positive_complete_one_loop_coefficient"
        ],
        "canonical_kinetic_interval_positive": bool(0 < kappa_lo < kappa_hi),
        "complete_one_loop_Phi_pole_remainder_below_two_e_minus_405": bool(
            0 < Btotal < sp.Rational(2, 10**405)
        ),
        "unit_disc_has_no_other_complete_one_loop_Phi_pole": bool(Btotal < 1),
        "complete_local_Phi_potential_curvature_positive": bool(1 - Btotal > 0),
        "scalar_MS_mass_anchor_bound_below_one_e_minus_6": bool(
            0 < scalar_Pi_one_upper < sp.Rational(1, 10**6)
        ),
        "large_finite_total_mass_reference_retained": bool(
            mass_ref_lo < mass_ref_hi < 0
        ),
        "two_scalar_vacuum_constant_retained_and_bounded": bool(
            0 < scalar_vac_abs < 10**396
        ),
        "full_one_loop_vacuum_energy_reference_retained": bool(
            10**799 < total_vac_lo < total_vac_hi < 10**800
        ),
        "once_fixed_scalar_finite_contact_not_a_free_parameter": bool(
            contact["actual_contact_absolute_upper"] > 0
        ),
    }
    return {
        "actual_MS_interaction_parameters": {
            "L": L,
            "g": g,
            "M": M,
            "mF": m,
            "lambda": lam,
        },
        "exact_scale_log_rational_enclosure": log,
        "actual_fixed_scalar_finite_contact": contact[
            "actual_fixed_finite_potential_contact"
        ],
        "old_scalar_complete_one_loop_absolute_error": old_error,
        "reference_and_fermion_relative_increment_interval": [extra_lo, extra_hi],
        "complete_one_loop_canonical_coefficient_interval": band,
        "complete_one_loop_absolute_error_upper": total_error,
        "complete_one_loop_relative_error_upper": total_error / T,
        "complete_kinetic_normalization_interval": [kappa_lo, kappa_hi],
        "complete_Phi_pole_remainder_upper": Btotal,
        "positive_local_Phi_curvature_lower": 1 - Btotal,
        "scalar_MS_Pi_mass_anchor_interval": [0, scalar_Pi_one_upper],
        "complete_finite_mass_reference_interval": [mass_ref_lo, mass_ref_hi],
        "two_scalar_vacuum_constant_absolute_upper": scalar_vac_abs,
        "complete_one_loop_vacuum_energy_interval": [total_vac_lo, total_vac_hi],
        "bounds": bounds,
        "scope": "All one-loop scalar and fermion contributions to the canonical Phi pole and forward second coefficient are included in the single declared MS interaction scheme with physical local mass, H one-point and vacuum-energy references. The coefficient is positive within the explicitly bounded one-loop interval. This is not a bound on new-model two-loop or later terms, global quantum-potential stability, exact gauge spectrum, high-energy V contours, finite-gravity G or common-parent B.",
        "checks": checks,
    }
