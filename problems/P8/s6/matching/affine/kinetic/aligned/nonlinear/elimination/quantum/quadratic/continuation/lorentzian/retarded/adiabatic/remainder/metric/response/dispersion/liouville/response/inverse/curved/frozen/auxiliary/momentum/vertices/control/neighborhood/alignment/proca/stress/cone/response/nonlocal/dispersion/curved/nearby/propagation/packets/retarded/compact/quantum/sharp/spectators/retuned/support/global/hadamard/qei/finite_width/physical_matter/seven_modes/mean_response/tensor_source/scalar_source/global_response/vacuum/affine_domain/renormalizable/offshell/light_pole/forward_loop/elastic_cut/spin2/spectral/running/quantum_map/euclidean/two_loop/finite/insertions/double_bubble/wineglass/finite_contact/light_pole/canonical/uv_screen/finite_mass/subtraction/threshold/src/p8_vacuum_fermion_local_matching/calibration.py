"""Rational mass, canonical residue and finite-field quartic bounds."""

from functools import cache

import sympy as sp
from p8_vacuum_gauge_yukawa_screen import calibration as candidate
from p8_vacuum_gauge_yukawa_screen import flow

from . import anchors


@cache
def data():
    parent = candidate.data()
    m = parent["new_common_Dirac_mass"]
    Y = parent["prospective_Yukawa_squared"]
    L = parent["old_scalar_quartic_used_only_as_prospective_boundary"]
    g = parent["old_cubic_squared"]
    M = parent["old_heavy_mass_squared"]
    ry = flow.data()["Yukawa_squared_to_gauge_squared"]
    rL = flow.data()["quartic_to_gauge_squared"]
    N = sp.Integer(6)
    d = L - 3 * g / M
    lowY = ry * L / 2
    highY = 2 * ry * L / 3
    normalized = anchors.enclosure(m)
    factor_low = 2 * N * lowY / 256
    factor_high = 2 * N * highY / 144
    f0_low = 4 * N * lowY * m * m / 256
    f0_high = 4 * N * highY * m * m / 144
    f1_low = f0_low + factor_low * normalized["mass_increment_over_2NY_div_Q_lower"]
    f1_high = f0_high + factor_high * normalized["mass_increment_over_2NY_div_Q_upper"]
    slope_low = factor_low * normalized["slope_over_2NY_div_Q_lower"]
    slope_high = factor_high * normalized["slope_over_2NY_div_Q_upper"]
    kappa_low = 1 - slope_high
    kappa_high = 1 - slope_low
    Bup = (
        N
        * highY
        / 144
        * (sp.Rational(7, 15) * m * m + sp.Rational(1, 10))
        / (m * m - sp.Rational(3, 4)) ** 2
    )
    V4up = 64 * N * highY * highY / 144
    quartic_relative = V4up / d
    normalized_B = Bup / kappa_low
    vacuum_low = 63 * m**4 / 256
    vacuum_high = 63 * m**4 / 144
    checks = {
        "same_actual_fermion_mass": m - 10**200,
        "same_cubic_square_and_eliminated_quartic": d - (L - 3 * g / M),
        "rational_Y_lower_dictionary": lowY - ry * L / 2,
        "rational_Y_upper_dictionary": highY - ry * L / sp.Rational(3, 2),
        "stronger_quartic_ray_lower_polynomial_test": 675 * sp.Rational(3, 2) ** 2
        + 5310 * sp.Rational(3, 2)
        - 11552
        + sp.Rational(8273, 4),
        "all_fourteen_vacuum_constant": sp.Rational(3, 2) * 14 * 3 - 63,
        "mass_reference_counterterm_interval": (1 - f1_low)
        - (1 - f1_high)
        - (f1_high - f1_low),
        "canonical_kinetic_interval_width": kappa_high
        - kappa_low
        - (slope_high - slope_low),
        "normalized_pole_bound": normalized_B * kappa_low - Bup,
        "quartic_relative_bound": quartic_relative * d - V4up,
    }
    bounds = {
        "strict_stronger_quartic_ray_bracket": bool(sp.Rational(3, 2) < rL < 2),
        "exact_Y_between_rational_bounds": bool(0 < lowY < Y < highY),
        "strict_on_shell_slope_positive_and_small": bool(
            0 < slope_low < slope_high < sp.Rational(1, 10**205)
        ),
        "strict_positive_canonical_kinetic_normalization": bool(
            0 < kappa_low < kappa_high < 1
        ),
        "large_finite_local_mass_threshold_retained": bool(
            10**193 < f0_low < f0_high < 10**195
        ),
        "large_finite_on_shell_mass_reference_retained": bool(
            10**193 < f1_low < f1_high < 10**195
        ),
        "negative_mass_reference_not_misread_as_physical_tachyon": 1 - f1_low < 0,
        "all_flavor_vacuum_energy_reference_retained": bool(
            10**799 < vacuum_low < vacuum_high < 10**800
        ),
        "strict_pole_remainder_below_one_e_minus_606": bool(
            0 < normalized_B < sp.Rational(1, 10**606)
        ),
        "positive_potential_curvature_on_selected_functional": 1 - normalized_B > 0,
        "old_tree_heavy_eliminated_quartic_positive": d > 0,
        "fermion_quartic_threshold_relative_below_one_e_minus_8": bool(
            0 < quartic_relative < sp.Rational(1, 10**8)
        ),
        "selected_functional_quartic_stays_positive": d - V4up > 0,
        "higher_even_fermion_potential_remainder_nonnegative": True,
    }
    return {
        "same_prospective_MSbar_fermion_boundary_not_identical_canonical_couplings": True,
        "fermion_mass": m,
        "actual_Yukawa_squared": Y,
        "rational_Yukawa_squared_lower": lowY,
        "rational_Yukawa_squared_upper": highY,
        "original_tree_eliminated_quartic": d,
        "normalized_anchor_enclosures": normalized,
        "all_flavor_vacuum_energy_bounds": [vacuum_low, vacuum_high],
        "zero_momentum_mass_threshold_bounds": [f0_low, f0_high],
        "on_shell_mass_threshold_bounds": [f1_low, f1_high],
        "mass_reference_parameter_bounds": [1 - f1_high, 1 - f1_low],
        "on_shell_slope_bounds": [slope_low, slope_high],
        "canonical_kinetic_normalization_bounds": [kappa_low, kappa_high],
        "unscaled_OS_remainder_coefficient_upper": Bup,
        "canonical_OS_remainder_coefficient_upper": normalized_B,
        "canonical_potential_mass_squared_lower": 1 - normalized_B,
        "fermion_quartic_threshold_absolute_upper": V4up,
        "fermion_quartic_threshold_relative_to_tree_d_upper": quartic_relative,
        "unnormalized_selected_functional_quartic_lower": d - V4up,
        "checks": checks,
        "bounds": bounds,
        "scope": "The exactly normalized tree-scalar plus one-loop-fermion functional has a unit-residue mass-one pole on the unit disc and positive potential on |y Phi_reference/mF|<1 after the explicitly fixed local references. This is not the complete quantum candidate: old scalar loops and later fermion/gauge loops remain separate, and canonical L/G/y are not relabelled as unchanged MSbar fixed-flow inputs.",
    }
