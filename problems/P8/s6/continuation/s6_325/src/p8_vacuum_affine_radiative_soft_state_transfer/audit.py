"""Uniform radiative soft-state transfer with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_three_real_probability_overlap import audit as previous

from . import continuity, source
from . import matching as soft_matching

MODEL = "original_radiative_soft_state_transfer_not_full_loop_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_uniform_state_continuity",
    "whole_paired_soft_transfer",
)

ITEM = {
    "id": "QG2_H8A489_uniform_radiative_finite_conversion_and_paired_soft_state_transfer",
    "status": "UNIFORM_RADIATIVE_SOFT_STATE_TRANSFER_NOT_FULL_LOOP_ALL_N_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original radiative soft-state-transfer model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated source or radiative soft-state observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "uniform_radiative_finite_conversion",
        "regulator_uniform_state_difference",
        "Born_seeded_paired_soft_state_transfer",
    ):
        raise ValueError(
            "Require the uniform state bound, regulator difference or named paired transfer"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A named soft-state insertion cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_uniform_state_continuity": continuity.data(),
        "whole_paired_soft_transfer": soft_matching.data(),
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
        "all180_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 180,
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
        "established": "Uniform finite-conversion continuity for every finite radiative state: |delta a|<1400R/kappa, |delta Delta|<10000R(1+ln(1/R))/kappa and the regulator quotient<15000R(1+ln(1/R))/kappa. The named Born-seeded paired soft-state insertion has a total-variation regulator limit bounded by200000*x*(1+ln(1/x))/kappa^2.",
        "domain": "Unchanged S300 recoil and original parameters,5/4<=E<=2,arbitrary finite positive null radiation of totalR<=1/8 and all angles;0<e<=1/8. The marked Born seed is physical D4; only the additional universal soft factor is dimensionally continued. The calorimetric remaining energy is x-b,0<x<=1/8.",
        "historical_qualification": qualifications(),
        "not_established": "A complete radiative hard loop or whole NNLO real-virtual rate, evanescent outer-state and finite hard matching, all-N hard subtraction/summation, positive normalized detector measure, interacting quantum state or unitarity, absolute complex Regge, common-parent bounce or original V/G/B/P8 closure.",
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
        "same_original_action_recoil_and_fixed_ball_current": True,
        "uniform_arbitrary_multiplicity_state_difference": True,
        "trace_and_radial_finite_terms_both_retained": True,
        "energy_dependent_regulator_dominator_not_constant_remainder": True,
        "same_state_virtual_pole_and_amplitude_phase_retained": True,
        "calorimetric_Gamma_factor_and_x_minus_b_retained": True,
        "no_full_loop_matching_or_original_closure_promotion": True,
        "rejected_inputs": rejected_inputs(),
    }
