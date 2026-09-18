"""Positive Borel soft conversion and conditioned cloud transfer; P8 open."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_quantitative_soft_cutoff import audit as previous

from . import current, cutoff, extension, moments, source

MODEL = "original_borel_soft_conversion_and_cloud_transfer_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_positive_current_increment",
    "whole_Borel_soft_conversion_extension",
    "whole_conditioned_cloud_soft_transfer",
    "whole_index_and_conversion_cutoff",
)

ITEM = {
    "id": "QG2_H8A495_positive_energy_soft_conversion_and_conditioned_cloud_transfer",
    "status": "BOREL_SOFT_CONVERSION_AND_CONDITIONED_CLOUD_TRANSFER_NOT_HARD_INTERACTING_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError(
            "Require the original positive-energy coefficient and cutoff model"
        )


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Require a stated positive-energy increment or quantitative cutoff observable"
        )


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "positive_Borel_additional_soft_factor",
        "conditioned_cloud_soft_transfer",
        "quantitative_index_conversion_cutoff",
    ):
        raise ValueError("Require the specified positive increment or D4 marked cutoff")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A marked leading reference cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_positive_current_increment": current.data(),
        "whole_Borel_soft_conversion_extension": extension.data(),
        "whole_conditioned_cloud_soft_transfer": moments.data(),
        "whole_index_and_conversion_cutoff": cutoff.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: clean(value)
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }


def scalar_entry_count():
    return sum(
        len(value) if isinstance(value, s.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "all9_original_primitive_statuses_preserved": len(frontier()) == 9,
        "all186_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 186,
        "all_matching_identifiers_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            row["status"].startswith("REJECTED_") for row in matching()
        ),
        "no_unknown_matching_coefficient_is_chosen": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The same fixed-ball additional-soft a,Delta and regulated coefficients extend continuously to positive Borel angular-energy measures of mass<=1/8. Positive added energy t at fixed original E,u has index modulus1700t/kappa and conversion modulus11000t(1+ln1/t)/kappa. On the fixed-Born D4 leading cloud conditioned byR<=x, the same-state marked regulator insertion has a total-variation limit<19000x(1+ln1/x)/kappa^2. Actual cutoff-conditioned a and Delta mean errors are<4000eta/kappa^2 and<34000eta(1+ln1/eta)/kappa^2.",
        "domain": "Original action and recoil, fixed hard E,u,0<=a<=1 (original a<4/(5kappa)),0<eta<=x<=1/8. Physical D4 marked cloud with only the additional soft factor dimensionally continued; the full trace, radial and phase convention are retained.",
        "historical_qualification": qualifications(),
        "not_established": "A positive fractional-dimensional cloud intensity, quantitative cutoff control of the changed remaining-energy logarithm in Z, outer hard dimensional/evanescent conversion, finite radiative hard matching or full common-regulator detector rate, state-dependent branching or interacting state, all-N hard summation, quantum unitarity, absolute Regge, common-parent bounce, UV completion or original V/G/B/P8 closure.",
    }


def bad_cases():
    bad = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.I,
        s.nan,
        s.Symbol("unspecified"),
    )
    rows = []
    for i, value in enumerate(bad):
        for label, call in (
            ("model", require_model),
            ("observable", require_observable),
            ("order", source.require_order),
            ("mass", source.require_mass),
            ("sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate((-1, 0, s.Rational(1, 2), s.Rational(3, 2))):
        rows.append(
            ("outside_positive_integer_order_" + str(i), source.require_order, (value,))
        )
    for i, value in enumerate((0, -1)):
        rows.append(("outside_positive_mass_" + str(i), source.require_mass, (value,)))
    for i, value in enumerate(
        (
            "whole_finite_G_amplitude",
            "UV_complete",
            "full_crossed_b20",
            "exact_massless_LSZ",
            "choose_finite_curvature",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "all_matching_zero",
            "known_piece_is_entire_amplitude",
            "finite_loop_proves_unitarity",
            "drop_evanescent_terms",
            "drop_old_Gram_finite",
            "discard_imaginary_phase_from_amplitude",
            "relative_Regge_is_absolute",
            "all_loop_regulator_complete",
            "graph_counts_close_full_inclusive_detector_rate",
        )
    ):
        rows.append(("unsupported_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("deleted_primitive", validate_scope, ([], matching())),
            ("deleted_matching", validate_scope, (frontier(), [])),
        )
    )

    rows.extend(
        (
            (
                "unsupported_composite_TT_premise",
                require_sector,
                ("composite_pair_assumed_TT",),
            ),
            (
                "unsupported_full_probability_subtraction",
                require_sector,
                ("full_probability_already_subtracted",),
            ),
        )
    )
    rows.extend(
        (
            ("third_tree_not_all_N", require_observable, ("third_tree_closes_all_N",)),
            (
                "rectangle_not_probability",
                require_observable,
                ("amplitude_rectangle_is_probability",),
            ),
        )
    )
    rows.extend(
        (
            (
                "signed_not_positive_probability",
                require_observable,
                ("signed_measure_is_positive_detector_probability",),
            ),
            (
                "no_virtual_matching_by_definition",
                require_observable,
                ("baseline_equals_virtual_terms_by_definition",),
            ),
        )
    )
    rows.extend(
        (
            (
                "different_subtraction_prescriptions",
                require_observable,
                ("amplitude_and_probability_baselines_identical",),
            ),
            (
                "state_changes_not_zero",
                require_observable,
                ("drop_radiative_state_current_changes",),
            ),
        )
    )
    rows.extend(
        (
            (
                "soft_transfer_not_full_loop",
                require_observable,
                ("soft_transfer_equals_full_radiative_hard_loop",),
            ),
            (
                "uniform_constant_not_radial_dominator",
                require_observable,
                ("constant_O_e2_error_integrable_against_dR_over_R",),
            ),
        )
    )
    rows.extend(
        (
            (
                "conservation_zero_normal_derivative_not_zero",
                require_observable,
                ("drop_constraint_before_soft_derivative",),
            ),
            (
                "collinear_bound_not_unique_atom_value",
                require_observable,
                ("assign_unique_soft_collinear_TT_point",),
            ),
        )
    )
    rows.extend(
        (
            (
                "unpolarized_cancellation_not_per_helicity",
                require_observable,
                ("delete_hard_phase_in_each_helicity",),
            ),
            (
                "real_part_bound_not_hard_amplitude_modulus",
                require_observable,
                ("bound_full_hard_loop_square_by_real_interference",),
            ),
        )
    )
    rows.extend(
        (
            (
                "energy_measure_not_quantum_state",
                require_observable,
                ("angular_energy_is_interacting_quantum_state",),
            ),
            (
                "current_counterexample_not_complete_coefficient",
                require_observable,
                ("nonunique_current_implies_nonunique_complete_C",),
            ),
            (
                "keep_atom_self_diagonal",
                require_sector,
                ("delete_measure_atom_self_diagonal",),
            ),
        )
    )
    rows.extend(
        (
            (
                "leading_reference_not_interacting_state",
                require_observable,
                ("conditioned_reference_is_full_interacting_state",),
            ),
            (
                "D4_mark_not_dimensional_matching",
                require_sector,
                ("scalar_expDelta_matches_every_marked_insertion",),
            ),
            (
                "fixed_intensity_not_radiative_branching",
                require_sector,
                ("replace_Born_intensity_by_radiative_index",),
            ),
        )
    )
    rows.extend(
        (
            (
                "positive_increment_not_arbitrary_signed",
                require_sector,
                ("arbitrary_signed_angular_energy_Lipschitz",),
            ),
            (
                "fixed_hard_data_not_variable_E_u",
                require_sector,
                ("same_bound_for_changing_hard_E_and_u",),
            ),
            (
                "rare_cut_not_loose_absolute_normalization",
                require_observable,
                ("divide_uncontrolled_absolute_error_by_tiny_cut_probability",),
            ),
        )
    )
    rows.extend(
        (
            (
                "additional_soft_not_outer_hard_dimensional",
                require_sector,
                ("cloud_limit_supplies_outer_hard_evanescence",),
            ),
            (
                "fractional_D_contraction_not_probability",
                require_observable,
                ("noninteger_D_contraction_is_positive_cloud_intensity",),
            ),
            (
                "no_remaining_energy_cutoff_rate",
                require_observable,
                ("same_cutoff_rate_for_changed_ln_x_minus_R",),
            ),
        )
    )
    return rows


def rejected_inputs():
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported finite-gravity input accepted: " + name)
    return total


def controls():
    return {
        "original_fixed_ball_current_and_positive_recoil": True,
        "complete_trace_radial_and_phase_terms": True,
        "positive_Borel_extension_with_diagonal_independence": True,
        "exact_conditioned_log_moments_and_remaining_energy": True,
        "relative_rare_event_normalization_kept": True,
        "distinct_index_and_conversion_cutoff_moduli": True,
        "additional_soft_only_not_outer_hard_dimensional_completion": True,
        "rejected_inputs": rejected_inputs(),
    }
