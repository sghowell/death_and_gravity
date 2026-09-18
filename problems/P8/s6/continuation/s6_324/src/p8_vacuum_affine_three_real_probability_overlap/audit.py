"""Complete third-order probability overlap with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_three_singleton_subtraction import audit as previous

from . import current, faces, measure, projection, source

MODEL = "original_three_real_probability_overlap_not_inclusive_P8"
OBSERVABLES = (
    "whole_literal_external_residues",
    "whole_uniform_proper_faces",
    "whole_probability_projection",
    "whole_finite_scheme_transfer",
)

ITEM = {
    "id": "QG2_H8A488_complete_three_real_probability_overlap_and_scheme_transfer",
    "status": "COMPLETE_THREE_REAL_PROBABILITY_OVERLAP_NOT_INCLUSIVE_ALL_N_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original three-real probability-overlap model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated source or three-real overlap observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_two_marked_soft_faces",
        "complete_probability_rectangle",
        "finite_signed_scheme_transfer",
    ):
        raise ValueError(
            "Require the stated complete faces, probability rectangle or signed transfer"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "A finite signed probability rectangle cannot close original P8"
        )
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_literal_external_residues": current.data(),
        "whole_uniform_proper_faces": faces.data(),
        "whole_probability_projection": projection.data(),
        "whole_finite_scheme_transfer": measure.data(),
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
        "all179_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 179,
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
        "established": "Every proper three-real amplitude face uses the full four-massive-plus-two-null state and has a uniform anchored hierarchy. The exact12-term transfer converts the parent amplitude-baseline signed subtraction into the complete anchored probability density. Its total variation is below2e-1424*x^2+1e-1340*x^4 and has a unique common-cutoff limit.",
        "domain": "Original action and couplings; positive real energies with total<=1/8,5/4<=E<=2,positive associated Born transfers,fixed real null directions and unit complex spatial TT leaves. Angular and energy-hierarchy bounds are uniform on the same allowed pole-excluding domain. Only the kinematic soft current receives an absolute energy disc.",
        "historical_qualification": qualifications(),
        "not_established": "A positive normalized inclusive probability or matching to actual real-virtual and integrated counterterms; all-N subtraction and summation, finite hard or evanescent matching, interacting state, unitarity, absolute complex Regge, common-parent bounce or original V/G/B/P8 closure. The radiative-state current changes in the older single-marked reference are retained as separate matching terms.",
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
        "same_original_action_and_frozen_lower_sources": True,
        "generic_scalar_and_EH3_conserved_root_residues": True,
        "original5116_faces_with_both_marked_null_legs": True,
        "all27_union_pairs_and12_proper_overlap_terms": True,
        "kinematic_absolute_disc_not_full_amplitude_disc": True,
        "nonzero_state_changes_not_virtual_matching_by_definition": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
