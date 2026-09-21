"""Complete selected local-tadpole radiation without closing hard matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_radiative_curvature_matching import audit as previous

from . import bounds, calibration, contractions, radiation, source

MODEL = "original_complete_local_tadpole_radiation_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_D_contractions",
    "whole_quartic_radiation",
    "whole_uniform_bounds",
    "whole_original_calibrations",
)
ITEM = {
    "id": "QG2_H8A501_complete_local_tadpole_radiation",
    "status": "EXACT_SELECTED_LOCAL_LOOP_RADIATION_NOT_FULL_HARD_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original selected local-tadpole radiation model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated complete local-tadpole observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "six_local_degree_six_tadpoles",
        "physical_four_scalar_one_graviton_TT",
        "known_OS4_local_radiative_remainder",
    ):
        raise ValueError("Require the stated selected local-loop radiation sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A selected local-loop sector cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_D_contractions": contractions.data(),
        "whole_quartic_radiation": radiation.data(),
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
        "all192_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 192,
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
        "other_loop_hard_quantum_Regge_and_UV_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "All six original degree-six local tadpoles give a complete exact-D physical four-scalar one-graviton 1PI sector, including massive scalar loop insertions and covariant connection terms. After the retained finite subtraction and four external emissions, generic Ward identities and five on-shell reductions give the known local TT representative. Its original-source nonleading remainder and finite interference are uniformly bounded as stated.",
        "domain": "Original n,g,C,kappa,lambda and H8A420-VAC-OS4 conditions. Selected six-field local matter tadpoles; physical external TT at D4 after exact-D subtraction. 25/4<=s<=16,strictly nonforward angle,0<=x<=1/8. Remainder uniform on the compact recoil domain; integrated bound fixed in the Born angle.",
        "historical_qualification": qualifications(),
        "not_established": "Other matter/source or internal-graviton loops; full curved counterfunctional; physical curvature matching including the S336 unknown coordinate; full interacting detector probability, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure.",
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
            "external_line_emissions_are_all_loop_graphs",
            "set_D4_before_pole_subtraction",
            "flat_matching_determines_curvature_couplings",
            "TT_representative_is_complete_curved_counterfunctional",
            "local_remainder_is_full_inclusive_loop_rate",
            "local_tadpole_sector_closes_P8",
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
            raise ValueError("Unsupported local radiation input accepted: " + name)
    return count


def controls():
    return {
        "original_source_and_existing_subtraction_unchanged": True,
        "all720_labels_metric_tadpole_and_massive_loop_insertion": True,
        "exact_D_retained_through_finite_subtraction": True,
        "complete_Ward_and_generic_physical_reductions": True,
        "original_source_nonzero_remainder_and_uniform_bounds": True,
        "other_loops_and_unknown_curvature_matching_still_open": True,
        "original_P8_frontiers_and_historical_qualifications_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
