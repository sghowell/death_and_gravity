"""A nonempty common Schwartz/Fourier class for the two distinct tree errors."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import calibration as parameters

from . import norms


@cache
def data():
    n = norms.data()
    p = parameters.point()
    C = n["actual_C_R"]
    light_radius = sp.Integer(1)
    mapped_radius = 3 * light_radius
    source_radius = 2 * mapped_radius
    r = 2 + source_radius**2
    op = p["cubic_squared"] * r**4 / (8 * p["D"] ** 5 * (1 - r / p["D"]))
    combined = n["unit_jet_cube_total_field_redefinition_error"] + op * (1 + C) ** 4
    return {
        "new_field_fourier_radius": light_radius,
        "mapped_field_fourier_radius": mapped_radius,
        "mapped_quadratic_source_fourier_radius": source_radius,
        "centered_box_spectral_radius": r,
        "mapped_source_L2_squared_factor": (1 + C) ** 4,
        "operator_error_coefficient_on_common_class": op,
        "combined_tree_action_error_coefficient": combined,
        "domain": "Real Schwartz Psi with Fourier support in the Euclidean four-momentum unit ball and integral |Psi_tilde| d^4p/(2pi)^4 <=1. Let U=||Psi||_2. All Psi jets have L-infinity<=1 and L2<=U. The map Phi=Psi+R(Psi) and J=Phi^2 are evaluated before the heavy resolvent.",
        "conclusion": "|S_polynomial[Phi,H_stationary]-S_mass_one_plus_actual_vacuum_quartic[Psi]| <= combined U^2, with the full nonlocal heavy resolvent retained and all boundary currents integrable. Not an all-field/all-momentum or loop statement.",
        "checks": {
            "cubic_map_Fourier_support_radius": mapped_radius - 3,
            "quadratic_mapped_source_support_radius": source_radius - 6,
            "Minkowski_centered_box_multiplier_radius": r - 38,
        },
        "bounds": {
            "common_spectral_window_inside_actual_heavy_gap": r < p["D"],
            "combined_nonzero_tree_action_error_below_one_e_minus_800": 0
            < combined
            < sp.Rational(1, 10**800),
            "mapped_source_norm_factor_below_two": (1 + C) ** 4 < 2,
        },
    }
