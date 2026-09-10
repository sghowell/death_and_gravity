"""Actual bounds for the first finite Phi normalization and its MS commutator."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_local_matching import calibration as fermion
from p8_vacuum_full_one_loop_matching import calibration as first

from . import bounds


@cache
def data():
    p = model.data()["actual_parameters"]
    L, g, M, lam = [
        p[k]
        for k in (
            "bare_polynomial_quartic",
            "cubic_coupling_squared",
            "heavy_mass_squared",
            "lambda",
        )
    ]
    f = fermion.data()
    Y = f["rational_Yukawa_squared_upper"]
    E = first.data()["complete_one_loop_relative_error_upper"]
    enclosure = bounds.enclosure(L, g, M, Y, E)
    k0 = enclosure["first_normalization_zero_absolute_upper"]
    k1 = enclosure["first_normalization_epsilon_absolute_upper"]
    cross = enclosure["first_normalization_second_amplitude_relative_upper"]
    comm = enclosure["coordinate_MS_commutator_b2_relative_upper"]
    return {
        "actual_parameters": p,
        "same_rational_Yukawa_upper": Y,
        "same_complete_canonical_one_loop_relative_upper": E,
        "actual_first_finite_field_enclosure": enclosure,
        "first_field_second_amplitude_absolute_upper": 4 * lam * cross,
        "coordinate_commutator_b2_absolute_upper": 4 * lam * comm,
        "checks": {
            "same_declared_fermion_scale": f["fermion_mass"] - 10**200,
            "amputated_one_loop_single_field_factor": enclosure[
                "one_loop_amputated_relative_upper"
            ]
            - E
            - 2 * k0,
            "first_field_amplitude_cross_sum": cross
            - 3 * k0**2
            - 2 * k0 * enclosure["one_loop_amputated_relative_upper"],
            "pole_commutator_tree_ratio": comm - L * k1 / 144,
            "fundamental_tree_coefficient": 4 * lam - 2 * g / (M - 2) ** 3,
            "tree_field_square_coefficient": enclosure[
                "first_normalization_tree_square_relative_upper"
            ]
            - 3 * k0**2,
        },
        "bounds": {
            "actual_scale_and_mass_domain": bool(32 < M < 10**400),
            "first_field_coefficient_below_four_e_minus_207": bool(
                k0 < s.Rational(4, 10**207)
            ),
            "first_epsilon_coefficient_below_four_e_minus_205": bool(
                k1 < s.Rational(4, 10**205)
            ),
            "same_complete_first_loop_bound_below_one_e_minus_6": bool(
                E < s.Rational(1, 10**6)
            ),
            "first_field_second_amplitude_relative_below_one_e_minus_212": bool(
                cross < s.Rational(1, 10**212)
            ),
            "coordinate_pole_commutator_relative_below_six_e_minus_412": bool(
                comm < s.Rational(6, 10**412)
            ),
            "first_field_tree_square_relative_below_five_e_minus_413": bool(
                enclosure["first_normalization_tree_square_relative_upper"]
                < s.Rational(5, 10**413)
            ),
            "first_normalization_only_factor_nonsingular": bool(1 - k0 > 0),
            "first_field_second_amplitude_absolute_below_four_e_minus_812": bool(
                4 * lam * cross < s.Rational(4, 10**812)
            ),
        },
        "scope": "Bound only the second-order terms determined by the first finite Phi normalization. The coordinate pole commutator is recorded separately and is not added twice to that physical-amplitude bound. The genuinely new second slope, other finite scheme directions and the full GY14 map remain open.",
    }
