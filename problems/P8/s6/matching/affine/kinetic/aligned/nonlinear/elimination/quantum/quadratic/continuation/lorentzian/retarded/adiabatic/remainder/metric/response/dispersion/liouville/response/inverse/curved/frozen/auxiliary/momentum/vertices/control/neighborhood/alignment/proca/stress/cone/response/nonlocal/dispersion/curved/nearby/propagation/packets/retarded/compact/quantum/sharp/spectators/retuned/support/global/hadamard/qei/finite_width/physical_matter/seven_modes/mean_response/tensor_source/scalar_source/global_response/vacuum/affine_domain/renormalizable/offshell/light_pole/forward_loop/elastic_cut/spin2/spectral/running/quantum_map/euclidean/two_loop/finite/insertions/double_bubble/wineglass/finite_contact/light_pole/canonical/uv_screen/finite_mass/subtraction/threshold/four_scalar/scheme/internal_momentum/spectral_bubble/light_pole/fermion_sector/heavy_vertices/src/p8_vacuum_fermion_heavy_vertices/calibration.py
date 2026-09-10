"""The full nonlocal-vertex family at the actual reference parameters."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_local_matching import calibration as fermion

from . import bounds


@cache
def data():
    p = model.data()["actual_parameters"]
    f = fermion.data()
    m, Y = f["fermion_mass"], f["rational_Yukawa_squared_upper"]
    L, g, M, lam = (
        p[k]
        for k in (
            "bare_polynomial_quartic",
            "cubic_coupling_squared",
            "heavy_mass_squared",
            "lambda",
        )
    )
    result = bounds.enclosure(m, Y, L, g, M, 144)
    E = result["full_family_disc_amplitude_and_b2_upper"]
    relative = E / (4 * lam)
    return {
        "actual_reference_parameters": {
            "mF": m,
            "Y_upper": Y,
            "L": L,
            "g": g,
            "M": M,
            "lambda": lam,
        },
        "literal_full_family_enclosure": result,
        "full_family_upper_relative_to_tree": relative,
        "checks": {
            "same_fixed_fermion_mass": m - 10**200,
            "same_fixed_tree_reference": 4 * lam - 4 * sp.Rational(1, 10**600),
            "same_positive_spectral_density_factor": result["spectral_density_C_upper"]
            - 12 * Y / 144,
            "all_three_contributions_summed": E
            - sum(
                result[k]
                for k in (
                    "all_channel_bubble_amplitude_upper",
                    "all_channel_triangle_amplitude_upper",
                    "all_channel_box_amplitude_upper",
                )
            ),
            "relative_reference": relative * (4 * lam) - E,
            "dyadic_actual_argument": 4 * result["spectral_threshold"] - 16 * 10**400,
        },
        "bounds": {
            "actual_parameters_within_proved_domain": bool(
                m >= 2 and M >= 24 and Y > 0 and L > 0 and g > 0
            ),
            "actual_dyadic_upper_covers_log_argument": bool(
                16 * m * m <= 2 ** result["dyadic_exponent"]
            ),
            "actual_previous_dyadic_power_is_too_small": bool(
                16 * m * m > 2 ** (result["dyadic_exponent"] - 1)
            ),
            "full_family_absolute_upper_below_one_e_minus_622": bool(
                0 < E < sp.Rational(1, 10**622)
            ),
            "full_family_relative_upper_below_one_e_minus_22": bool(
                0 < relative < sp.Rational(1, 10**22)
            ),
        },
        "scope": "Complete paired W2-W2-F0 family in the stated local subtraction scheme, not the complete enlarged-model order-two matching correction.",
    }
