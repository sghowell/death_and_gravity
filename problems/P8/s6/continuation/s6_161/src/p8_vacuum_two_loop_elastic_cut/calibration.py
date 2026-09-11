"""Actual full first-amplitude and complete second elastic-cut allowances."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_elastic_cut import amplitude as tree_cut
from p8_vacuum_fermion_local_matching import calibration as fermion
from p8_vacuum_finite_field_covariance import calibration as field
from p8_vacuum_full_one_loop_matching import calibration as first
from p8_vacuum_full_two_loop_amplitude import calibration as second
from p8_vacuum_two_loop_physical_source_map import calibration as physical_source

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
    m, Y = f["fermion_mass"], f["rational_Yukawa_squared_upper"]
    ell = first.data()["exact_scale_log_rational_enclosure"]["scale_log_upper"]
    k0 = field.data()["actual_first_finite_field_enclosure"][
        "first_normalization_zero_absolute_upper"
    ]
    scalar = bounds.scalar_enclosure(L, g, M, ell, 1024)
    Bs = scalar["full_scalar_MS_amplitude_absolute_upper"]
    Bf = bounds.fermion_enclosure(m, Y)
    Bfield = 146 * k0 * lam
    B = Bs + Bf + Bfield
    cut = bounds.cut_enclosure(lam, B)
    old = second.data()["complete_canonical_two_loop_enclosure"]
    band = bounds.improved_band(lam, old["total_relative"], B)
    rel = cut["integrated_second_absolute_upper"] / (4 * lam)
    return {
        "actual_full_scalar_one_loop_amplitude_enclosure": scalar,
        "actual_full_fermion_one_loop_amplitude_upper": Bf,
        "actual_first_canonical_field_amplitude_upper": Bfield,
        "complete_one_loop_physical_amplitude_upper": B,
        "complete_second_elastic_cut_enclosure": cut,
        "complete_second_elastic_cut_tree_relative_upper": rel,
        "complete_through_two_loop_improved_formal_band": band,
        "bounds": {
            "actual_positive_reference_parameters": bool(min(L, g, M, lam, m, Y) > 0),
            "heavy_logarithm_native_power_cap_valid": bool(M < 2**1024),
            "actual_thresholds_outside_window": bool(M > 6 and 4 * m * m > 6),
            "actual_dirac_Neumann_mass_gap": bool(m >= 24),
            "complete_one_loop_amplitude_below_one_e_minus_405": bool(
                B < s.Rational(1, 10**405)
            ),
            "second_cut_relative_below_one_e_minus_407": bool(
                rel < s.Rational(1, 10**407)
            ),
            "complete_improved_formal_relative_error_below_one_e_minus_6": bool(
                band["total_relative_error"] < s.Rational(1, 10**6)
            ),
            "complete_improved_formal_uniform_lower_positive": band[
                "positive_formal_uniform_lower"
            ],
        },
        "checks": {
            "same_actual_tree_lambda_as_physical_cut": lam
            - tree_cut.data()["actual_fixed_lambda"],
            "same_actual_tree_forward_coefficient": 4 * lam - 2 * g / (M - 2) ** 3,
            "full_one_loop_amplitude_includes_all_three_groups": B - Bs - Bf - Bfield,
            "single_canonical_field_factor": Bfield - 2 * k0 * 73 * lam,
            "second_cut_relative_dictionary": rel - s.Rational(73, 5120) * B,
            "same_parent_complete_two_loop_error": old["total_relative"]
            - physical_source.data()[
                "physical_source_b2_one_plus_two_loop_relative_upper"
            ],
            "full_improved_band_includes_both_cut_orders": band["total_absolute_error"]
            - 4 * lam * old["total_relative"]
            - cut["integrated_first_upper"]
            - cut["integrated_second_absolute_upper"],
            "cut_subtracted_lower_difference": old["formal_uniform_lower"]
            - band["formal_uniform_lower"]
            - cut["integrated_first_upper"]
            - cut["integrated_second_absolute_upper"],
        },
        "scope": "The full first amplitude need not be small relative to the exceptionally canceled tree amplitude. Only its absolute bound and the resulting cut-to-b2 allowance are used. This is the complete through-two-loop Phi cut on [4,6], not the physical remainder, full low-energy gauge cut at later orders or V/G/B closure.",
    }
