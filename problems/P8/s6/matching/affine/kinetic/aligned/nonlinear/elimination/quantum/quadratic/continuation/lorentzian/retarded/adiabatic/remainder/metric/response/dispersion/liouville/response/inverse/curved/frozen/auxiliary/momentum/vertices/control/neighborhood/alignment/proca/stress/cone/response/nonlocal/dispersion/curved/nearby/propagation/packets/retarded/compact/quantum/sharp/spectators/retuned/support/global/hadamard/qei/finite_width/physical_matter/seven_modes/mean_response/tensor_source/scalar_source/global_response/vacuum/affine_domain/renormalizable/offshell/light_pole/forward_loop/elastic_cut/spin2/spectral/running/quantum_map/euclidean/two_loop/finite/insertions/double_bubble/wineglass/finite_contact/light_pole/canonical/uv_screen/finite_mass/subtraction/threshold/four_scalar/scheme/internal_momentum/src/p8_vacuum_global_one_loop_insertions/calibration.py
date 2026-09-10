"""Exact actual half-plane insertion and finite reference-window bounds."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_local_matching import calibration as fermion
from p8_vacuum_full_one_loop_matching import calibration as combined
from p8_vacuum_full_one_loop_matching import logarithm
from p8_vacuum_light_pole import calibration as scalar

from . import halfplane, routing, window


@cache
def data():
    f = fermion.data()
    c = combined.data()
    m = f["fermion_mass"]
    Yhi = f["rational_Yukawa_squared_upper"]
    kappa_lo, kappa_hi = c["complete_kinetic_normalization_interval"]
    rhi = scalar.point(2)["finite_kinetic_counterterm_upper"]
    ell = logarithm.enclosure()
    ell_hi = ell["scale_log_upper"]
    C_hi = 2 * halfplane.N * Yhi / 144
    hp = halfplane.enclosure(m, Yhi, 144, 1)
    B = hp["unscaled_OS_insertion_uniform_upper"]
    Bcan = B / kappa_lo
    K = 4 * m * m - 1
    energy = sp.Integer(10) ** 400
    tmax = energy * energy
    log_argument = 1 + (tmax + 1) / K
    w = window.defect_bound(rhi, C_hi, ell_hi, kappa_lo)
    epsilon = w["fractional_inverse_defect_upper"]
    route = routing.channel_margin(2, 1, 1)
    coarse = C_hi / K
    checks = {
        "same_actual_fermion_mass": m - 10**200,
        "same_canonical_kinetic_lower_includes_fermion_slope": kappa_lo
        - 1
        + f["on_shell_slope_bounds"][1],
        "reference_energy_equal_to_fermion_mass_squared": energy - m * m,
        "reference_log_argument_strict_gap_identity": sp.factor(
            m * m - log_argument - (3 * m**4 - 5 * m * m) / K
        ),
        "canonical_halfplane_bound_not_unscaled": Bcan * kappa_lo - B,
        "global_coarse_insertion_bound": coarse * K - C_hi,
        "finite_window_both_sectors_retained": epsilon * kappa_lo - rhi - C_hi * ell_hi,
        "reference_window_normalized_inverse_lower": w["normalized_inverse_lower"]
        + epsilon
        - 1,
        "same_selected_channel_routing_margin": route["halfplane_margin"]
        - sp.Rational(1, 4),
    }
    bounds = {
        "actual_halfplane_parameter_gap_strictly_positive": bool(
            hp["parameter_gap_lower"] > 0
        ),
        "actual_Y_upper_positive": bool(Yhi > 0),
        "actual_canonical_kinetic_lower_positive": bool(0 < kappa_lo < kappa_hi),
        "unscaled_global_halfplane_insertion_below_one_e_minus_606": bool(
            0 < B < sp.Rational(1, 10**606)
        ),
        "canonical_global_halfplane_insertion_below_one_e_minus_606": bool(
            0 < Bcan < sp.Rational(1, 10**606)
        ),
        "halfplane_bound_sharper_than_old_small_disc_coefficient": bool(
            B < f["unscaled_OS_remainder_coefficient_upper"]
        ),
        "halfplane_bound_sharper_than_global_coarse_log_bound": bool(B < coarse),
        "selected_bubble_all_internal_momenta_covered": route[
            "selected_symmetric_bubble_routing_covered"
        ],
        "selected_bubble_has_strict_halfplane_margin": bool(
            route["halfplane_margin"] > 0
        ),
        "finite_reference_log_argument_positive_and_below_m_squared": bool(
            1 < log_argument < m * m
        ),
        "reference_scale_log_enclosure_below_one_thousand": bool(
            900 < ell["scale_log_lower"] < ell_hi < 1000
        ),
        "scalar_slope_bound_retained_and_positive": bool(rhi > 0),
        "finite_reference_window_fractional_defect_below_one_e_minus_202": bool(
            0 < epsilon < sp.Rational(1, 10**202)
        ),
        "complete_selected_one_loop_inverse_positive_on_reference_window": w[
            "strictly_positive_on_declared_log_window"
        ],
    }
    return {
        "fermion_mass": m,
        "reference_energy_not_a_derived_cutoff": energy,
        "Euclidean_momentum_squared_window": [0, tmax],
        "actual_rational_Y_upper": Yhi,
        "complete_canonical_kinetic_interval": [kappa_lo, kappa_hi],
        "all_complex_halfplane_insertion_enclosure": hp,
        "canonical_halfplane_insertion_upper": Bcan,
        "spacelike_coarse_uniform_insertion_upper": coarse,
        "selected_symmetric_bubble_routing": route,
        "finite_reference_window_maximum_log_argument": log_argument,
        "inherited_rational_scale_log_enclosure": ell,
        "scalar_slope_upper": rhi,
        "fermion_logarithm_factor_upper": C_hi,
        "complete_one_loop_reference_window_inverse_enclosure": w,
        "bounds": bounds,
        "scope": "The same full scalar-fermion one-loop reference functional is controlled on this explicit finite Euclidean momentum window. This is neither a Wilsonian cutoff determination nor a bound on new-model two-loop, higher-loop, gauge-spectrum, V-contour, G or B errors.",
        "checks": checks,
    }
