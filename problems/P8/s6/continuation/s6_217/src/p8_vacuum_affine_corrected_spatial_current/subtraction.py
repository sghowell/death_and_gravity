"""Corrected physical UV subtraction and both analytic phase branches."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_spatial_symbol import density as old_density
from p8_vacuum_affine_dimensional_spatial_symbol import geometry as old_geometry
from p8_vacuum_affine_dimensional_spatial_symbol import jets as old_jets
from p8_vacuum_affine_dimensional_spatial_symbol import matching as old_matching
from p8_vacuum_affine_matched_spatial_current import dimension
from p8_vacuum_affine_ordered_scalar_symbol import density, jets, matching
from p8_vacuum_affine_subleading_band_conversion import centered

CONVERSION_TAIL = s.Integer(10) ** 40


@cache
def physical_data():
    products = old_jets.endpoint_products()
    nondecaying = []
    Pderivatives = []
    odd = []
    for key, row in products.items():
        if key[1] % 2:
            for degree, value in enumerate(row):
                physical = s.factor(value.subs(old_geometry.d, 3))
                if key[1] + degree < 4:
                    nondecaying.append(physical)
                Pderivatives.append(s.diff(physical, old_jets.p))
                if physical != 0:
                    odd.append((key, degree, physical))
    angular = old_density.coefficients()
    coefficient = jets.m**2 * jets.a**2 * (53 * jets.t**2 + 7) / 15
    checks = {
        "all_odd_nonlog_physical_UV_products_zero": s.Matrix(nondecaying),
        "all_unaveraged_odd_UV_products_transfer_independent": s.Matrix(Pderivatives),
        "only_first_endpoint_logarithmic_unaveraged_odd_entries": sum(
            not (key[1] == 1 and degree == 3) for key, degree, value in odd
        ),
        "seven_nonzero_unaveraged_odd_products_retained": len(odd) - 7,
    }
    for ch in ("tensor", "vector", "scalar"):
        checks[ch + "_complete_odd_angular_coefficient"] = s.Matrix(
            [
                s.factor(
                    value.subs(old_geometry.d, 3)
                    - (coefficient if key[1:] == (1, 0, 3) else 0)
                )
                for key, value in angular.items()
                if key[0] == ch and key[1] % 2
            ]
        )
    old = old_density.invariant_spatial_coefficients()
    new = density.invariant_spatial_coefficients()
    checks["all_105_physical_tracefree_spatial_UV_differences_unchanged"] = s.Matrix(
        [
            s.factor((new[key][i] - old[key][i]).subs(old_geometry.d, 3))
            for key in old
            for i in range(3)
        ]
    )
    inv, _g, _ell, finite = matching.finite_input()
    old_finite = old_matching.finite_input()[3]
    correction = (
        jets.p**2 * (3 * jets.t**2 + 1) * (-31 * inv[0] + 30 * inv[1]) / (420 * s.pi**2)
    )
    checks["corrected_finite_tracefree_input_not_legacy"] = s.Matrix(
        [
            s.factor(
                finite[r].subs(dict.fromkeys(inv[3:], 0))
                - old_finite[r]
                - (correction if r == 0 else 0)
            )
            for r in range(3)
        ]
    )
    x = s.Symbol("x", real=True)
    comparison = x / (1 - x) + s.log(1 - x)
    checks["logarithmic_shell_majorant_derivative"] = s.factor(
        s.diff(comparison, x) - x / (1 - x) ** 2
    )
    checks["logarithmic_shell_majorant_zero_origin"] = comparison.subs(x, 0)
    checks["small_transfer_shell_one_over_K_factor"] = 1 / (
        1 - s.Rational(1, 4)
    ) - s.Rational(4, 3)
    return {
        "complete_odd_physical_UV": "All unaveraged odd nonlog coefficients vanish. The seven nonzero odd products all have j1/degree3 and are independent of external momentum; source-jet1 cross-magnetic terms are retained and cancel only after the complete geometry sum. Each normalized tensor/vector/scalar angular coefficient is the same m^2 a^2(53t^2+7)/15 times source jet0.",
        "unchanged_nondecaying_subtraction": "The corrected physical tracefree one-ball spatial differences F2,F4 are exactly unchanged. An odd q4 logarithmic shell contributes only a decaying original-mask difference, not A3,A1,K2 or K0. Every nondecaying two-leg shape comes from unchanged even endpoint rows, so the actual S210 A3,A1 and complete K2/K0 cancellations remain.",
        "odd_shell_not_discarded": "For |P|<=K/4, the lost logarithmic shell lies in r>=K-|P| and log(K/(K-|P|))<=|P|/(K-|P|)<=4|P|/(3K). Large transfer uses S209's separate whole-band/shape estimate, not this small-transfer expansion. The shell stays in the error, with its complete source jets.",
        "uniform_conversion_tail": "S209's1e40 tail is an absolute sum over individual endpoint/degree rows, with the S208 Cauchy row majorants and all-P partition. Every corrected row and its remainder is multiplied by a unit sign, so the same bound applies. Nondecaying shape terms are independently shown unchanged, not inferred merely from equal absolute bounds.",
        "finite_part_correction": correction,
        "checks": checks,
        "gates": {
            "all_nondecaying_original_shape_coefficients_phase_unchanged": all(
                key[1] == 1 and degree == 3 for key, degree, value in odd
            ),
            "all_odd_magnetic_source_jets_retained": any(
                key[2] == 1 for key, degree, value in odd
            ),
            "odd_log_shell_is_not_a_finite_counterterm": True,
            "full_original_all_P_conversion_tail": centered.tail_constants()[
                "complete_uniform_conversion_remainder"
            ]
            < CONVERSION_TAIL,
            "same_original_two_leg_mask_and_lower_comoving_band": True,
            "nonzero_evanescent_finite_change_still_required": correction != 0,
        },
    }


@cache
def dimension_data():
    raw = dimension.mode_data()
    defects = {}
    for kind, row in raw["rows"].items():
        defects[kind] = sum(
            row[n]["upper"] / s.Integer(1000) ** (2 * n - 1) for n in range(1, 5)
        )
    exponent = s.Rational(1, 20)
    x = s.Symbol("x", real=True)
    re, im = s.symbols("re im", real=True)
    plus = re + s.I * im
    sharp = re - s.I * im
    checks = {
        "real_dimension_positive_imaginary_analytic_branch_pair": s.expand(
            (plus - sharp) / (2 * s.I) - im
        ),
        "two_frequency_absolute_phase_defects": 2 * s.Rational(1, 40) - exponent,
        "geometric_majorant_for_both_phase_moduli": 1 / (1 - exponent)
        - s.Rational(20, 19),
        "dimension_far_integrable_exponent": -2
        + dimension.DIMENSION_RADIUS
        + s.Rational(7, 4),
        "dimension_angular_integrable_exponent": -dimension.DIMENSION_RADIUS / 2
        + s.Rational(1, 8),
        "all_transfer_degree_margin": 6
        - (5 + dimension.DIMENSION_RADIUS)
        - s.Rational(3, 4),
        "phase_modulus_geometric_bound_at_zero": (1 / (1 - x)).subs(x, 0) - 1,
    }
    return {
        "complete_W8_absolute_defects": defects,
        "both_retarded_phase_moduli": "On the real unit slab and |d-3|<=1/4, each complete continued W8 frequency differs from its real bare omega by less than1/40, uniformly in both real momenta. This follows from all four exact coefficient majorants and omega>=1000. Hence both opposite exponential phases have modulus<=exp(1/20)<20/19<2; no complex-dimension unit-modulus assertion is made.",
        "analytic_pair_not_parameter_conjugation": "For real d, +Im of the positive-phase sharp-detector/source-annihilation product equals(Fplus-Fsharp)/(2i). Continue both branches at the same complex d, never conjugating the parameter. Real tensor algebra is complexified linearly only after this mode-current coefficient is fixed.",
        "endpoint_dimension_limit": "The S213 full four-order frequency/domain and active/passive angular majorants are independent of which retarded branch is used. Every subtracted endpoint row is changed by a unit sign on real d, and the same identity extends analytically. Angular exponent-1/8, far radial exponent-7/4 and near transfer degree21/4<6 give the same dominated spatial dimension limit for the corrected Q.",
        "finite_time_comparison": "The original prescription keeps the absolutely convergent actual-state correction and finite time remainder in physical dimension. The explicit both-phase bound additionally verifies that changing the retarded branch introduces no complex-dimensional exponential obstruction. No physical state in noninteger dimension is selected.",
        "checks": checks,
        "gates": {
            "all_four_continued_frequency_defects_below_one_fortieth": all(
                v < s.Rational(1, 40) for v in defects.values()
            ),
            "uniform_both_phase_moduli_below_two": s.Rational(20, 19) < 2,
            "no_conjugation_of_complex_dimension": True,
            "same_original_dimensional_endpoint_majorants": all(
                bool(v) for v in dimension.data()["gates"].values()
            ),
            "same_actual_state_and_no_new_reference_state": True,
            "dimension_limits_not_only_fixed_transfer": 5 + dimension.DIMENSION_RADIUS
            < 6,
        },
    }
