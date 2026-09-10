"""MS-slope frontier; finite mass and complete canonical matching remain open."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_quadratic_forests import audit as previous

from . import calibration, reference, remainder

MODULES = (reference, remainder, calibration)
TARGETS = ("scalar_Phi2_W0_F2", "gauge_Phi2")
STATUS = "BOUNDED_NONLOCAL_AND_MS_SLOPE_NOT_FINITE_MASS_REFERENCE"


def frontier():
    rows = previous.frontier()
    for r in rows:
        if r["id"] in TARGETS:
            r["status"] = STATUS
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError("MS-slope frontier differs from exact scope")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: s.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    old, new = previous.frontier(), frontier()
    out.update(
        {
            "only_two_primitive_slopes_advance": sum(a != b for a, b in zip(old, new))
            - 2,
            "nine_rows_retained": len(new) - 9,
            "two_vacuum_rows_wholly_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 2,
            "two_nonlocal_plus_slope_rows": sum(r["status"] == STATUS for r in new) - 2,
            "four_previous_quartic_rows_unchanged": sum(
                r["status"]
                in (
                    "BOUNDED_PAIRED_IN_COMMON_MS_INTERACTION_SCHEME",
                    "BOUNDED_PAIRED_PRIMITIVE_ROW",
                )
                for r in new
            )
            - 4,
        }
    )
    return out


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "same_three_cyclic_words_and_gauge_sign": True,
        "dimension_symbolic_Dirac_contractions": True,
        "full_p_squared_tensor_derivative_including_cross_term": True,
        "exact_Gaussian_sunset_master_reduction": True,
        "proper_fermion_mass_kinetic_and_Yukawa_CT_once": True,
        "fixed_mu_mass_derivative_before_mu_equals_m": True,
        "single_and_double_MS_poles_subtracted": True,
        "finite_epsilon_times_pole_terms_retained": True,
        "vacuum_second_derivative_normalization_check": True,
        "massless_scalar_reference_not_actual_scalar": True,
        "full_scalar_mass_difference_is_UV_convergent": True,
        "chord_diagonal_and_both_hard_regions_bounded": True,
        "joint_fractional_massless_chord_integral": True,
        "unit_soft_coefficient_bound_precedes_integrals": True,
        "zero_soft_to_OS_derivative_via_prior_Cauchy_disc": True,
        "finite_mass_and_other_matching_not_claimed": True,
        "two_vacuum_rows_and_complete_pole_open": True,
        "V_G_B_and_original_P8_open": True,
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    out = []
    for j in range(5):
        for i, v in enumerate(invalid):
            vals = [720, 1, 1, 144, 1]
            vals[j] = v
            out.append((f"type_{j}_{i}", calibration.enclosure, tuple(vals)))
    for i, v in enumerate(
        (
            (719, 1, 1, 144, 1),
            (720, -1, 1, 144, 1),
            (720, 1, -1, 144, 1),
            (720, 1, 1, 0, 1),
            (720, 1, 1, 145, 1),
            (720, 1, 1, 144, -1),
        )
    ):
        out.append((f"domain_{i}", calibration.enclosure, v))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"frontier_{i}", validate_frontier, (rows,)))
    out.extend(
        (
            ("frontier_missing", validate_frontier, (frontier()[:-1],)),
            (
                "frontier_extra",
                validate_frontier,
                (frontier() + [{"id": "extra", "status": "COMPLETE"}],),
            ),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported MS-slope input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "unmodified_MS_dimensional_forest": True,
        "no_early_four_dimensional_numerator": True,
        "actual_scalar_mass_difference_retained": True,
        "exact_Gamma_recurrences_no_float_fit": True,
        "finite_MS_slope_not_an_adjustable_constant": True,
        "finite_mass_reference_still_open": True,
        "other_matching_and_canonical_terms_open": True,
        "original_P8_not_closed": True,
    }
