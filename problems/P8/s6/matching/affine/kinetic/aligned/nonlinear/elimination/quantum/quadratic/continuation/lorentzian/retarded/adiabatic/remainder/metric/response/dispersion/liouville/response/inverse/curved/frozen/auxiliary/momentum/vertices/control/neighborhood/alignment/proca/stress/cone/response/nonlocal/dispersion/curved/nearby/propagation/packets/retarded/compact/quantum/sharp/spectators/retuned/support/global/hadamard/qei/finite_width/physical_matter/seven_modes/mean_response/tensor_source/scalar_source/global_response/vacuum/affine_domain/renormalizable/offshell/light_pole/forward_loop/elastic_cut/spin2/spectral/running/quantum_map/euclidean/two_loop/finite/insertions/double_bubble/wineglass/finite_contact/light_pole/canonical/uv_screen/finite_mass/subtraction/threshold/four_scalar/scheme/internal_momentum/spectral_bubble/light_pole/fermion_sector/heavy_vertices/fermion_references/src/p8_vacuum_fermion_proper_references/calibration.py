"""Actual shared-reference bounds for the proper one-loop fermion subgraphs."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_local_matching import calibration as fermion
from p8_vacuum_full_one_loop_matching import calibration as complete
from p8_vacuum_gauge_yukawa_screen import calibration as candidate

from . import bounds


@cache
def data():
    f = fermion.data()
    c = candidate.data()
    same = complete.data()
    m = f["fermion_mass"]
    Yhi = f["rational_Yukawa_squared_upper"]
    L = c["old_scalar_quartic_used_only_as_prospective_boundary"]
    ahi = 2 * L / 3
    kappa = same["complete_kinetic_normalization_interval"]
    Phi_abs = max(abs(k - 1) for k in kappa)
    result = bounds.enclosure(m, Yhi, ahi, Phi_abs, 144)
    mass = result["selected_mass_ratio_minus_one_absolute_upper"]
    yuk = result["selected_Yukawa_ratio_minus_one_absolute_upper"]
    z = result["finite_fermion_kinetic_increment_upper"]
    return {
        "actual_reference_parameters": {
            "mF": m,
            "actual_a": c["prospective_gauge_squared"],
            "actual_Y": f["actual_Yukawa_squared"],
            "rational_a_upper": ahi,
            "rational_Y_upper": Yhi,
        },
        "same_complete_Phi_kinetic_interval": kappa,
        "actual_proper_fermion_reference_enclosure": result,
        "checks": {
            "same_fermion_mass_reference": m - 10**200,
            "rational_gauge_upper_from_same_quartic_ray": ahi - 2 * L / 3,
            "same_complete_Phi_field_interval": Phi_abs
            - max(abs(k - 1) for k in kappa),
            "exact_scalar_internal_mass_ratio": result[
                "scalar_boson_to_fermion_mass_squared_ratio"
            ]
            - sp.Rational(1, 10**400),
            "finite_mass_conversion_bound_sum": mass
            - result["finite_mass_relative_anchor_absolute_upper"]
            - z,
            "finite_Yukawa_conversion_retains_scalar_field": yuk
            - 2
            * (result["finite_Yukawa_relative_anchor_absolute_upper"] + z + Phi_abs),
        },
        "bounds": {
            "actual_gauge_coupling_below_rational_upper": bool(
                0 < c["prospective_gauge_squared"] < ahi
            ),
            "actual_Yukawa_coupling_below_rational_upper": bool(
                0 < f["actual_Yukawa_squared"] < Yhi
            ),
            "shared_scalar_field_bound_below_one_half": bool(
                0 < Phi_abs < sp.Rational(1, 2)
            ),
            "positive_selected_fermion_kinetic_reference": bool(
                1
                <= result["fermion_kinetic_interval"][0]
                <= result["fermion_kinetic_interval"][1]
            ),
            "positive_selected_mass_reference": result[
                "positive_selected_mass_reference"
            ],
            "positive_selected_Yukawa_reference": result[
                "positive_selected_Yukawa_reference"
            ],
            "actual_mass_reference_relative_shift_below_one_e_minus_204": bool(
                0 < mass < sp.Rational(1, 10**204)
            ),
            "actual_Yukawa_reference_relative_shift_below_one_e_minus_204": bool(
                0 < yuk < sp.Rational(1, 10**204)
            ),
            "actual_fermion_field_shift_below_one_e_minus_204": bool(
                0 < z < sp.Rational(1, 10**204)
            ),
        },
        "scope": "Finite one-loop local fermion references at zero external momenta, with the actual common physical Phi field. Not a charged-fermion on-shell spectrum, completed insertion integrals, the S6.138 outer-reference conversion or complete order-two matching.",
    }
