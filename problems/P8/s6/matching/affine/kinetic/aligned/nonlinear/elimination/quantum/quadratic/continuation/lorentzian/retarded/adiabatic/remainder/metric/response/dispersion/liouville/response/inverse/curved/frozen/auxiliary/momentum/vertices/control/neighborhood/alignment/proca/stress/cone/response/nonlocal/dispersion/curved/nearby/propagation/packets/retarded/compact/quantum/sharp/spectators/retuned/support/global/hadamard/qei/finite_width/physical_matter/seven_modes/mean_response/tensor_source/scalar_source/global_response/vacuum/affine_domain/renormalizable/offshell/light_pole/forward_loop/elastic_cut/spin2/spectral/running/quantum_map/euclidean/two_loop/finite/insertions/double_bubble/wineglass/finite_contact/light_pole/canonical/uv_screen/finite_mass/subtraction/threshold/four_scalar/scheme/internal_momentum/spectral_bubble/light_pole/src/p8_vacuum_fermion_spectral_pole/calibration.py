"""Actual spectral quadratic-family and selected-functional interval bounds."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_local_matching import calibration as fermion
from p8_vacuum_full_one_loop_matching import calibration as combined
from p8_vacuum_global_one_loop_insertions import calibration as global_one_loop
from p8_vacuum_global_one_loop_insertions import window

from . import enclosure


@cache
def data():
    p = model.data()["actual_parameters"]
    f = fermion.data()
    one = combined.data()
    global_one = global_one_loop.data()
    m, Y, g, M = (
        f["fermion_mass"],
        f["rational_Yukawa_squared_upper"],
        p["cubic_coupling_squared"],
        p["heavy_mass_squared"],
    )
    d = enclosure.bound(m, Y, g, M, 144)
    r2 = d["finite_positive_outer_slope_upper"]
    B2 = d["uniform_outer_OS_remainder_coefficient_upper"]
    klo, khi = one["complete_kinetic_normalization_interval"]
    selected_kinetic = [klo, khi + r2]
    Bselected = one["complete_Phi_pole_remainder_upper"] + B2 / klo
    window_bound = window.defect_bound(
        global_one["scalar_slope_upper"] + r2,
        global_one["fermion_logarithm_factor_upper"],
        global_one["inherited_rational_scale_log_enclosure"]["scale_log_upper"],
        klo,
    )
    old_defect = global_one["complete_one_loop_reference_window_inverse_enclosure"][
        "fractional_inverse_defect_upper"
    ]
    new_defect = window_bound["fractional_inverse_defect_upper"]
    checks = {
        "same_actual_fermion_mass": m - 10**200,
        "same_actual_cubic_squared": g - sp.Rational(1, 67108864),
        "actual_family_remainder_normalization": sp.factor(
            B2 - 6 * Y * g / (6 * 144**2 * m**4)
        ),
        "selected_kinetic_interval_contains_named_second_slope": selected_kinetic[1]
        - khi
        - r2,
        "selected_complex_pole_bound_retains_new_family": (
            Bselected - one["complete_Phi_pole_remainder_upper"]
        )
        * klo
        - B2,
        "selected_Euclidean_window_retains_new_family": (new_defect - old_defect) * klo
        - r2,
    }
    bounds = {
        "actual_mixed_mass_reference_gap_positive": bool(M > 1 and m >= 2),
        "strict_actual_positive_spectral_family": d[
            "strict_family_slope_and_spacelike_remainder_positive"
        ],
        "finite_second_order_slope_below_one_e_minus_616": bool(
            0 < r2 < sp.Rational(1, 10**616)
        ),
        "uniform_second_order_remainder_below_one_e_minus_1017": bool(
            0 < B2 < sp.Rational(1, 10**1017)
        ),
        "selected_kinetic_interval_strictly_positive": bool(
            0 < selected_kinetic[0] < selected_kinetic[1]
        ),
        "selected_radius_two_remainder_below_two_e_minus_405": bool(
            0 < Bselected < sp.Rational(2, 10**405)
        ),
        "selected_unit_disc_has_no_other_light_pole": bool(Bselected < 1),
        "selected_local_Phi_curvature_positive": bool(1 - Bselected > 0),
        "selected_finite_Euclidean_window_defect_below_one_e_minus_202": bool(
            0 < new_defect < sp.Rational(1, 10**202)
        ),
        "selected_finite_Euclidean_window_inverse_positive": window_bound[
            "strictly_positive_on_declared_log_window"
        ],
    }
    return {
        "actual_reference_parameters": {"mF": m, "Y_upper": Y, "g": g, "M": M},
        "quadratic_covariance_insertion_family_enclosure": d,
        "selected_kinetic_normalization_interval": selected_kinetic,
        "selected_radius_two_OS_remainder_coefficient_upper": Bselected,
        "selected_local_Phi_curvature_lower": 1 - Bselected,
        "selected_reference_window_energy_not_a_cutoff": global_one[
            "reference_energy_not_a_derived_cutoff"
        ],
        "selected_reference_window_inverse_enclosure": window_bound,
        "bounds": bounds,
        "scope": "The complete one-loop quadratic functional plus this one first covariance-insertion sector is controlled after its explicit local references. This does not complete the enlarged-model two-loop pole or its MS interaction conversion, let alone higher-loop errors, cutoff, V/G/B or P8.",
        "checks": checks,
    }
