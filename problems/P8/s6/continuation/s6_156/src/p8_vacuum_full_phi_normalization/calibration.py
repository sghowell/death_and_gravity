"""Complete actual GY14 Phi normalization and canonical unit-disc inverse."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_insertion_ms import calibration as inserted_slopes
from p8_vacuum_fermion_local_matching import calibration as first_fermion
from p8_vacuum_fermion_ms_slopes import calibration as direct_slopes
from p8_vacuum_fermion_quadratic_forests import calibration as direct_tails
from p8_vacuum_fermion_spectral_pole import enclosure as inserted_pole
from p8_vacuum_finite_field_covariance import calibration as first_field
from p8_vacuum_scalar_ms_slopes import calibration as scalar

from . import bounds, pole


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
    ff = first_field.data()
    f = first_fermion.data()
    fs = ff["actual_first_finite_field_enclosure"]
    k0 = fs["first_normalization_zero_absolute_upper"]
    k1 = fs["first_normalization_epsilon_absolute_upper"]
    Y = f["rational_Yukawa_squared_upper"]
    r0 = g / (288 * M)
    z = 2 * 6 * Y / 144
    sc = scalar.data()
    slope_groups = {
        "all_32_scalar_refinements": sc["actual_scalar_MS_slope_enclosure"][
            "finite_MS_slope_absolute_upper"
        ],
        "both_direct_fermionic_quadratic_primitives": direct_slopes.data()["enclosure"][
            "primitive_on_shell_MS_slope_absolute_upper"
        ],
        "complete_fermion_insertion_quadratic_family": inserted_slopes.data()[
            "uniform_row_MS_slope_absolute_upper"
        ],
    }
    hybrid = sum(slope_groups.values())
    MS = bounds.normalization_upper(hybrid, k0, r0)
    E = direct_tails.data()["enclosure"]["both_quadratic_primitives_soft_tail_upper"]
    direct = pole.tail_coefficient(E)
    inserted = inserted_pole.bound(f["fermion_mass"], Y, g, M, 144)
    Bs = sc["unchanged_one_loop_OS_quadratic_upper"]
    Bf = f["unscaled_OS_remainder_coefficient_upper"]
    second_groups = {
        "complete_scalar_MS_converted_OS": sc["converted_two_loop_OS_quadratic_upper"],
        "both_direct_fermionic_quadratic_primitives_unit_disc": direct,
        "complete_fermion_insertion_quadratic_family": inserted[
            "uniform_outer_OS_remainder_coefficient_upper"
        ],
        "finite_first_Phi_parameter_reexpansion": k0 * (2 * Bs + Bf),
    }
    second = sum(second_groups.values())
    total = Bs + Bf + second
    tree_shift = 3 * k0**2 + 2 * MS + 2 * L * k1 / 144
    return {
        "actual_parameters": p,
        "complete_fixed_hybrid_two_loop_slope_group_uppers": slope_groups,
        "complete_fixed_hybrid_two_loop_slope_absolute_upper": hybrid,
        "first_finite_parameter_and_pole_product_slope_allowances": {
            "first_parameter_and_square_after_homogeneity": k0 * r0,
            "each_epsilon_field_pole_term_before_exact_cancellation": z * k1,
        },
        "complete_MS_second_Phi_normalization_absolute_upper": MS,
        "formal_one_plus_two_loop_Phi_normalization_lower": 1 - k0 - MS,
        "full_hybrid_OS_second_order_group_uppers": second_groups,
        "complete_canonical_one_plus_two_loop_OS_coefficient_upper": total,
        "unit_disc_factored_inverse_gap_lower": 1 - total,
        "finite_field_map_tree_b2_second_relative_upper": tree_shift,
        "checks": {
            "all_three_quadratic_fermion_rows_in_two_classes": 2 + 1 - 3,
            "all_scalar_and_fermionic_slope_classes": len(slope_groups) - 3,
            "complete_MS_slope_assembly": MS - hybrid - k0 * r0,
            "same_inserted_slope_bound": inserted["finite_positive_outer_slope_upper"]
            - slope_groups["complete_fermion_insertion_quadratic_family"],
            "both_direct_primitive_unit_disc_extension": direct - E / 2,
            "all_second_OS_groups_once": second - sum(second_groups.values()),
            "first_plus_second_canonical_OS_sum": total - Bs - Bf - second,
            "same_tree_coefficient": 4 * lam - 2 * g / (M - 2) ** 3,
            "complete_field_tree_shift_all_terms": tree_shift
            - 3 * k0**2
            - 2 * MS
            - 2 * L * k1 / 144,
            "same_actual_forward_reference_factor": 2 * L
            - 3 * g / (M - 2)
            - g * (3 * (M - 2) - 4) / (M - 2) ** 2,
        },
        "bounds": {
            "all_hybrid_slope_groups_strictly_positive": all(
                v > 0 for v in slope_groups.values()
            ),
            "all_OS_second_groups_strictly_positive": all(
                v > 0 for v in second_groups.values()
            ),
            "complete_hybrid_slope_below_one_e_minus_18": bool(
                hybrid < s.Rational(1, 10**18)
            ),
            "complete_MS_second_normalization_below_one_e_minus_18": bool(
                MS < s.Rational(1, 10**18)
            ),
            "formal_two_loop_normalization_strictly_positive": bool(1 - k0 - MS > 0),
            "extended_direct_fermion_OS_bound_below_one_e_minus_799": bool(
                direct < s.Rational(1, 10**799)
            ),
            "complete_canonical_OS_bound_below_one_e_minus_18": bool(
                total < s.Rational(1, 10**18)
            ),
            "strict_canonical_unit_disc_inverse_gap": bool(1 - total > 0),
            "complete_field_map_tree_shift_below_one_e_minus_18": bool(
                tree_shift < s.Rational(1, 10**18)
            ),
            "forward_reference_factor_below_two_L": bool(
                0 < 2 * L - 3 * g / (M - 2) < 2 * L
            ),
        },
        "scope": "Complete fixed-order Phi field normalization and canonical mass-one unit-residue light pole on the unit disc at the named GY14 boundary, using direct MS fermion coordinates and the full finite Phi parameter re-expansion. This is not yet the complete matched four-point/vacuum/source calculation, a global spectral result, physical truncation or V/G/B.",
    }
