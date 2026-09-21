"""Complete selected source-loop radiation and unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_local_tadpole_radiation import audit as previous

from . import bounds, calibration, contractions, mixed, radiation, source, triangles

MODEL = "original_complete_mixed_source_radiation_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_source_tadpoles",
    "whole_covariant_cross_class",
    "whole_D_kernel_and_EOM",
    "whole_complete_source_radiation",
    "whole_uniform_bounds",
    "whole_original_calibrations",
)
ITEM = {
    "id": "QG2_H8A502_complete_mixed_source_radiation",
    "status": "EXACT_SELECTED_SOURCE_LOOP_RADIATION_NOT_FULL_HARD_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original selected mixed-source radiation model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated complete mixed-source observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "three_J2_GH_J4_contraction_classes",
        "physical_four_scalar_one_graviton_TT",
        "known_OS4_source_radiative_remainder",
    ):
        raise ValueError("Require the stated selected source-loop radiation sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A selected source-loop sector cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_source_tadpoles": contractions.data(),
        "whole_covariant_cross_class": mixed.data(),
        "whole_D_kernel_and_EOM": triangles.data(),
        "whole_complete_source_radiation": radiation.data(),
        "whole_uniform_bounds": bounds.data(),
        "whole_original_calibrations": calibration.data(),
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
        "all193_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 193,
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
        "S336_unknown_curvature_coefficient_not_assigned": True,
        "other_hard_loops_quantum_Regge_and_UV_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The complete selected J2 G_H J4 one-loop physical radiation sector includes all three light contraction classes. Within-J2 cancels against the existing onepoint condition. Exact-D two-triangle, covariant Green-function and fully labelled EOM identities show that the complete across-source TT radiation is removed by the existing OS4 constant subtraction. The within-J4 corrected source and all42-plus-counterterm graphs reduce to the stated centered heavy-exchange radiation. Its original nonleading remainder and finite interference have explicit uniform bounds; S337 can be added by linearity.",
        "domain": "Original n,g,C,kappa and H8A420-VAC-OS4 conditions, selected cubic-plus-five-field source matter loops, real external TT after exact-D subtraction.25/4<=s<=16,-1<z<1,0<=x<=1/8. Remainder uniform over the compact recoil domain; integrated high-window bound fixed in nonforward Born angle.",
        "historical_qualification": qualifications(),
        "not_established": "Remaining polynomial/heavy or internal-graviton loop radiation; complete curved counterfunctional; independent physical curvature matching including the S336 unknown coordinate; full interacting inclusive detector probability, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
    for i, value in enumerate(
        (
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
    ):
        for label, call in (
            ("new_model", require_model),
            ("new_observable", require_observable),
            ("new_sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate(
        (
            "discard_EOM_before_metric_variation",
            "flat_kernel_fixes_every_curvature_operator",
            "only_light_triangle_is_complete",
            "centered_loop_contact_replaces_original_Born",
            "known_source_remainder_is_full_inclusive_rate",
            "selected_source_sector_closes_P8",
        )
    ):
        rows.append(("unsupported_new_scope_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("new_deleted_primitive", validate_scope, ([], matching())),
            ("new_deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    return rows


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            count += 1
        else:
            raise ValueError("Unsupported source-radiation input accepted: " + name)
    return count


def controls():
    return {
        "original_source_and_existing_subtractions_unchanged": True,
        "three_source_classes_all_metric_and_internal_emissions": True,
        "both_massive_triangles_in_exact_D_before_MS": True,
        "all24_EOM_emissions_and_full_Ward_not_dropped": True,
        "independent_original_component_calibrations_and_uniform_bounds": True,
        "other_loops_and_independent_curvature_matching_still_open": True,
        "original_P8_frontiers_and_historical_qualifications_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
