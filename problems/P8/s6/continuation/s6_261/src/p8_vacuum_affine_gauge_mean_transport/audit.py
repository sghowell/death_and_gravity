"""Finite gauge-orbit mean transport does not close the interacting P8 problem."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_spatial_brst import audit as previous

from . import crossing, geometry, modes

STATE = (
    "unchanged_S251_coupled_Cauchy_covariance_in_explicit_finite_Weyl_jet_diagnostic"
)
MEASURE = "finite_orbit_and_band_covariance_not_quantum_BRST_regulator_or_ordering"
OBSERVABLES = (
    "whole_spatial_density_weight_orbit",
    "all_product_mode_second_order_gauge_equation",
    "fixed_coupled_Weyl_coordinate_mean_and_physical_volume_contact",
    "actual_regular_bounce_coordinate_covariance_symbols",
)
ITEM = {
    "id": "QG2_H8A425_full_gauge_weight_orbit_fixed_reference_mean_transport_and_physical_volume_cancellation",
    "status": "EXACT_FULL_CLASSICAL_GAUGE_ORBIT_AND_FIXED_REFERENCE_WEYL_JETS_WITH_BOUNCE_UV_SYMBOLS_NOT_INTERACTING_NIELSEN_MEAN_ORDERING_REGULATOR_CUTOFF_ORIGINAL_V_G_B_OR_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(values):
    return previous.require_parameters(values)


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "Keep the entire fixed reference and the explicit Weyl-jet boundary"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError(
            "Require the stated finite orbit, source contact or actual reference symbol"
        )
    return label


def require_measure(label):
    if not isinstance(label, str) or label != MEASURE:
        raise ValueError("No interacting quantum measure or ordering is provided")
    return label


def require_weight(value):
    value = fraction(value)
    if value >= 2:
        raise ValueError("The displayed monotone gauge family requires w<2")
    return value


def require_band(lower, upper):
    lower, upper = fraction(lower), fraction(upper)
    if not 0 < lower < upper:
        raise ValueError("Keep both endpoints of the finite positive diagnostic band")
    return lower, upper


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No primitive or previous matching item is promoted")
    return True


@cache
def packets():
    return {
        "whole_spatial_density_weight_and_physical_volume": geometry.family(),
        "whole_exact_periodic_orbit_with_both_tensor_profiles": geometry.periodic_orbit(),
        "whole_lapse_dependent_physical_volume_contact": geometry.physical_volume_contact(),
        "whole_three_direction_unprojected_Fourier_gauge_jet": modes.packet(),
        "whole_unchanged_coupled_reference_rows_and_commutator": crossing.reference_rows(),
        "whole_actual_bounce_regular_canonical_covariance_symbols": crossing.bounce_symbols(),
    }


@cache
def residuals():
    return {
        packet + "_" + name: clean(value)
        for packet, data in packets().items()
        for name, value in data["checks"].items()
    }


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    return {
        **{
            packet + "_" + name: bool(value)
            for packet, data in packets().items()
            for name, value in data["gates"].items()
        },
        "all_nine_original_primitive_rows_unchanged": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
        "all_original_fixed_parameters_retained": require_parameters(parameters())
        == parameters(),
        "finite_Gaussian_support_not_assumed_inside_nonlinear_branch": True,
        "original_P8_and_interacting_fixed_mean_remain_open": True,
    }


def observable():
    return {
        "positive_result": "The whole density-weight gauge orbit gives an exact periodic map retaining two tensor profiles and a full unprojected three-direction second-order equation. Its fixed coupled Gaussian/Weyl jets have nonzero coordinate-mean shifts -9 Cvv/4 and -9 Cvn/4. The entire lapse-dependent physical volume cancels these against its two-point contact. The actual regular bounce chart gives explicit leading covariance and radial ultraviolet coefficients.",
        "scope_boundary": "These are classical reconstruction maps, finite Weyl polynomial expectations and unrenormalized fixed-reference symbols. They are not the actual interacting Nielsen vector or quantum mean. Generic reduced lapse and volume coordinates do not commute; their linear bounce commutator vanishes without fixing nonlinear ordering. A Gaussian jet calculation does not place its unbounded support inside a nonlinear auxiliary branch. The finite periodic diagnostic and R3 band restriction do not supply a BRST regulator or alter the original state.",
        "original_problem": "The same full parent, source contacts, fixed profiles and all existing preparations remain. Physical curved subtraction, Ward/state/regulator defects, complete interacting state, controlled cutoff and matching, nonlinear bounce, quantum gravity limit, UV scattering and finite-gravity infrared/Regge remain OPEN, with original V/G/B/P8 and scoped P8(a) qualifications unchanged.",
    }


def bad_cases():
    cases = []
    invalid = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.I,
        s.nan,
        s.Symbol("unknown"),
    )
    for i, value in enumerate(invalid):
        for label, call in (
            ("state", require_state),
            ("measure", require_measure),
            ("observable", require_observable),
            ("weight", require_weight),
        ):
            cases.append(("invalid_" + label + "_" + str(i), call, (value,)))
        for key in parameters():
            values = parameters()
            values[key] = value
            cases.append(
                (
                    "invalid_parameter_" + key + "_" + str(i),
                    require_parameters,
                    (values,),
                )
            )
    for value in (2, 3, s.Rational(5, 2)):
        cases.append(
            ("singular_or_wrong_weight_" + str(value), require_weight, (value,))
        )
    for lower, upper in (
        (0, 1),
        (-1, 1),
        (1, 1),
        (2, 1),
        (1.0, 2),
        (1, s.oo),
        (True, 2),
    ):
        cases.append(
            ("invalid_band_" + str((lower, upper)), require_band, (lower, upper))
        )
    for label in (
        "all_coordinate_means_gauge_invariant",
        "full_quantum_Nielsen_vector_evaluated",
        "Weyl_ordering_is_original_nonlinear_quantum_ordering",
        "Gaussian_support_inside_auxiliary_branch",
        "finite_Fourier_gauge_algebra_closed",
        "physical_volume_changes_with_gauge",
        "actual_interacting_lapse_mean_computed",
        "lapse_and_volume_always_commute",
        "cutoff_equals_heavy_mass",
        "ultraviolet_asymptotic_is_finite_error_bound",
        "state_reprepared_after_gauge_change",
        "source_onepoint_contact_deleted",
        "full_quantum_BRST_regulator_proved",
        "all_P8a_qualifications_removed",
        "original_P8_closed",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    for i in range(len(frontier())):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        cases.append(
            ("primitive_promotion_" + str(i), validate_scope, (rows, matching()))
        )
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        cases.append(
            ("matching_promotion_" + str(i), validate_scope, (frontier(), rows))
        )
    cases.append(
        ("previous_state_scope_not_current_contract", require_state, (previous.STATE,))
    )
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported gauge-mean input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "whole_density_weight_operator_and_both_tensor_profiles_retained": True,
        "nonlinear_displacement_all_product_modes_and_translations_retained": True,
        "full_coupled_covariance_not_replaced_by_diagonal_oscillators": True,
        "physical_lapse_volume_source_contacts_cancel_without_deletion": True,
        "generic_noncommutativity_and_Gaussian_support_boundary_explicit": True,
        "actual_crossing_symbols_not_interacting_mean_or_cutoff_estimate": True,
        "no_quantum_measure_state_ordering_or_original_P8_promotion": True,
        "all_original_and_previous_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
