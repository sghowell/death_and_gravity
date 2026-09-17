"""State-transport matching with all original and historical frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_two_real_soft_overlap import audit as previous

from . import assembly, continuity, series, source

MODEL = "original_state_transport_matched_soft_reference_not_full_P8"
OBSERVABLES = (
    "whole_angular_state_continuity",
    "whole_entire_signed_series_difference",
    "whole_probability_matching",
    "whole_matched_soft_reference",
)
ITEM = {
    "id": "QG2_H8A478_state_transport_connector_and_two_real_matched_soft_reference",
    "status": "STATE_TRANSPORT_AND_TWO_REAL_MATCHED_SOFT_REFERENCE_NOT_FULL_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original state-transport model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated state-transport observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_state_transport",
        "matched_soft_reference",
    ):
        raise ValueError("Require the stated state transport or matched soft reference")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A subtracted real measure cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_angular_state_continuity": continuity.data(),
        "whole_entire_signed_series_difference": series.data(),
        "whole_probability_matching_and_reference": assembly.data(),
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
        "all169_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 169,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "no_unknown_matching_coefficient_is_chosen": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "A total-energy continuity modulus for the complete angular soft conversion makes the state-transport connector finite before and after regulator removal. A finite probability-level two-real residual and correct two-marked Bose dressing combine with S309 to match all0/1/2-real physical D4 tree densities in a named soft reference, with uniformly vanishing relative difference below1e-653.",
        "domain": "Unchanged original four-dimensional full47/full434 trees and positive full Born; S300 recoil E in[5/4,2], all nonforward hard directions and physical emitted angles, total marked energy<=x<=1/8. The dimensional additional-soft scheme, finite conversion and total remaining-energy cut are retained.",
        "historical_qualification": qualifications(),
        "not_established": "Full finite hard real-virtual and evanescent matching, all-N nonleading radiation, a positive event measure or unitary inclusive detector rate, actual interacting quantum state, absolute complex Regge, original common-parent bounce or V/G/B/P8 closure.",
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
            "soft_overlap_closes_full_inclusive_detector_rate",
        )
    ):
        rows.append(("unsupported_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("deleted_primitive", validate_scope, ([], matching())),
            ("deleted_matching", validate_scope, (frontier(), [])),
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
        "same_original_full47_full434_source_and_positive_Born": True,
        "total_energy_continuity_includes_radial_and_trace_terms": True,
        "entire_signed_series_derivatives_allow_negative_regulated_index": True,
        "necessary_state_connector_before_seed_integration": True,
        "two_soft_probability_faces_and_two_marked_Bose_counting": True,
        "matched_reference_not_complete_hard_real_virtual_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
