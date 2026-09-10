"""Actual finite outer reference and the converted full-family error."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_heavy_vertices import calibration as parent

from . import remainder


@cache
def data():
    previous = parent.data()
    p = previous["actual_reference_parameters"]
    ref = remainder.enclosure(p["mF"], p["Y_upper"], 144)
    D = p["M"] - 2
    tree = 4 * p["lambda"]
    relative_coeff = -2 * p["L"] + 3 * p["g"] / D
    correction = abs(relative_coeff) * tree * ref["finite_F_absolute_upper"]
    old = previous["literal_full_family_enclosure"][
        "full_family_disc_amplitude_and_b2_upper"
    ]
    total = old + correction
    return {
        "actual_reference_parameters": p,
        "finite_outer_MS_reference_enclosure": ref,
        "b2_reference_conversion_absolute_upper": correction,
        "reference_conversion_upper_relative_to_tree": correction / tree,
        "full_family_b2_upper_in_common_MS_interaction_scheme": total,
        "converted_family_upper_relative_to_tree": total / tree,
        "checks": {
            "same_actual_threshold": ref["spectral_threshold"] - 4 * 10**400,
            "same_actual_loop_prefactor": ref["C_over_Q_upper"]
            - 12 * p["Y_upper"] / 144**2,
            "same_tree_forward_coefficient": 2 * p["g"] / D**3 - tree,
            "conversion_absolute_relative_dictionary": correction / tree
            - abs(relative_coeff) * ref["finite_F_absolute_upper"],
            "original_family_plus_conversion": total - old - correction,
            "actual_dyadic_exponent": ref["least_dyadic_exponent"] - 1331,
        },
        "bounds": {
            "actual_mass_ratio_in_proved_domain": bool(
                0 < ref["mass_ratio"] <= sp.Rational(1, 16)
            ),
            "least_dyadic_upper_covers_threshold": bool(
                ref["spectral_threshold"] <= 2 ** ref["least_dyadic_exponent"]
            ),
            "previous_dyadic_power_is_too_small": bool(
                ref["spectral_threshold"] > 2 ** (ref["least_dyadic_exponent"] - 1)
            ),
            "actual_dimensionless_correction_below_one": bool(
                0 < ref["dimensionless_finite_correction_absolute_upper"] < 1
            ),
            "finite_reference_upper_below_one_e_minus_206": bool(
                0 < ref["finite_F_absolute_upper"] < sp.Rational(1, 10**206)
            ),
            "finite_mass_ratio_correction_below_one_e_minus_600": bool(
                0 < ref["finite_F_correction_absolute_upper"] < sp.Rational(1, 10**600)
            ),
            "reference_conversion_relative_below_one_e_minus_400": bool(
                0 < correction / tree < sp.Rational(1, 10**400)
            ),
            "converted_family_b2_upper_below_one_e_minus_622": bool(
                0 < total < sp.Rational(1, 10**622)
            ),
            "converted_family_relative_upper_below_one_e_minus_22": bool(
                0 < total / tree < sp.Rational(1, 10**22)
            ),
        },
        "scope": "The full S6.138 family now has its outer-reference interaction conversion bounded; this is not the complete matched order-two amplitude or canonical-field correction.",
    }
