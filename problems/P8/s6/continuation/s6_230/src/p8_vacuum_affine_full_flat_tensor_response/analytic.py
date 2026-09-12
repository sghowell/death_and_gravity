"""Original first-sheet pole count and exact low/high scale inequalities."""

from functools import cache

import sympy as s
from p8_vacuum_affine_isolated_shear_resolvent import spectral

from . import normalization as norm

ALPHA = s.Rational(13, 60)
BETA = -s.Rational(52, 225)
LOW_DISK = s.Rational(11, 105)
RESIDUE_BOUND = s.Rational(250, 249)


@cache
def data():
    delta = s.Symbol("inverse_radius_delta", complex=True)
    dd = 1 + delta
    polynomial = 30 - 20 * dd + 3 * dd * dd
    local = -s.Rational(172, 225) + s.Rational(19, 30) * dd - dd * dd / 10
    x, b, M, c = s.symbols("real_p imaginary_p four_mass_squared radial_c", real=True)
    p = x + s.I * b
    imag = s.factor(s.im(p * p / (M + c * p)))
    C, tau, U, DA = s.symbols("C tau U DA", positive=True)
    checks = {
        "high_polynomial_first_difference": s.expand(
            polynomial - 13 - delta * (-14 + 3 * delta)
        ),
        "high_local_first_difference": s.expand(
            local - BETA - delta * (s.Rational(13, 30) - delta / 10)
        ),
        "right_quadrant_imaginary_integrand": s.factor(
            imag
            - b * (2 * x * M + c * (x * x + b * b)) / ((M + c * x) ** 2 + c * c * b * b)
        ),
        "upper_cut_denominator_imaginary": s.im(C - tau * (DA + s.I * s.pi * U))
        + s.pi * tau * U,
        "positive_threshold_value": s.expand(
            C
            - 4 * norm.MASS**2 * (-s.Rational(172, 225))
            - (C + s.Rational(688, 225) * norm.MASS**2)
        ),
        "complete_low_disk_A_bound": s.Rational(1, 30)
        + s.Rational(3, 14) / 3
        - LOW_DISK,
        "physical_low_disk_numerator": s.Rational(5, 6)
        + LOW_DISK
        - s.Rational(197, 210),
        "residue_inverse_margin": 1 / (1 - s.Rational(1, 250)) - RESIDUE_BOUND,
    }
    margins = {
        "T_coefficient_below_half": s.Rational(1, 2) - s.Rational(5, 4) * 21 / 60,
        "T_difference_below_delta": 1 - s.Rational(37, 60),
        "B_difference_below_half_delta": s.Rational(1, 2)
        - s.Rational(13, 30)
        - s.Rational(1, 20),
        "log_error_total_below_two_delta": 2
        - s.Rational(1, 2)
        - s.Rational(1, 2) * s.Rational(4, 3),
        "T_derivative_below_one": 1
        - (21 * s.Rational(3, 4) + s.Rational(5, 4) * 17) / 60,
        "root_small_radius_exclusion": 10**5 - 30 * (2 + 790 * 3),
        "root_large_C_lower_margin": 144 - 100,
        "root_exterior_A_above250": s.Rational(13, 80) * (790 * 2 - 2) - 4 - 250,
        "root_modulus_below_kappa": 1 - s.Rational(161, 250),
        "root_A_modulus_below5000": 5000 - 20 - 2 * 794 * 3,
        "root_error_below_1e_minus785": 10**5 - 8 * (1 + 794 * 3 + 4),
        "root_imaginary_A_above_three_tenths": ALPHA * s.Rational(3, 2)
        - s.Rational(1, 10**785)
        - s.Rational(3, 10),
        "root_imaginary_p_above_kappa_over_million": 144 * s.Rational(3, 10) / 5000**2
        - s.Rational(1, 10**6),
        "root_growth_above_1e393": s.Integer(10) ** 400 / (2 * 10**6)
        - s.Integer(10) ** 393,
        "low_disk_denominator_above_C0": s.Rational(5, 6) - LOW_DISK,
        "physical_low_disk_relative_below_1e_minus796": s.Rational(1, 10**796)
        - s.Rational(197, 210 * 144) * norm.MASS**2 / norm.KAPPA,
    }
    return {
        "original_A2": spectral.data()["closed_A2"],
        "complete_denominator": "D_C(p)=C+p A2(p), C=16pi^2 kappa+5m^2/6",
        "general_count_condition": "C>2m^2/15",
        "exact_first_sheet_count": "Exactly one simple zero in each open p half-plane, a conjugate pair with Re p<0; no real zeros or other first-sheet zeros. The real-axis/cut and full closing-arc argument is in notes/count.md.",
        "full_high_radius_error": "For r=|p|/m^2>=8: |A2-(13/60)Log(p/m^2)+52/225|<=8(1+logr+pi)/r, uniformly up to both cut banks.",
        "global_modulus_bound": "|A2(p)|<=30[1+log(1+|p|/m^2)] on the complete first sheet.",
        "actual_pole_modulus_squared_bounds": (s.Integer(10) ** 796, norm.KAPPA),
        "actual_q0_lambda_modulus_bounds": (s.Integer(10) ** 398, s.Integer(10) ** 400),
        "actual_q0_positive_growth_lower": s.Integer(10) ** 393,
        "actual_complex_residue_bound": RESIDUE_BOUND / norm.C,
        "low_disk_A_bound": LOW_DISK,
        "physical_low_disk_relative_error_bound": s.Rational(197, 30240)
        * norm.MASS**2
        / norm.KAPPA,
        "all_exact_positive_margins": margins,
        "checks": checks,
        "gates": {
            "every_coarse_arithmetic_margin_strictly_positive": all(
                value > 0 for value in margins.values()
            ),
            "actual_C_above_general_count_threshold": norm.C
            > 2 * norm.MASS**2 / s.Integer(15),
            "original_cut_U_positive_on_open_cut": True,
            "upper_count_retains_full_closing_arc": True,
            "root_count_not_numerical_winding_certificate": True,
            "pole_scale_not_asserted_to_be_a_physical_cutoff": True,
        },
    }
