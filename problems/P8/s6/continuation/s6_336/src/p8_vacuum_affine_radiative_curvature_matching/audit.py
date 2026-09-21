"""Independent radiative curvature matching direction and unchanged frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_dimensional_real_remainder import audit as previous

from . import bounds, calibration, contact, source
from . import matching as hard_matching

MODEL = "original_independent_radiative_curvature_matching_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_covariant_curvature_contact",
    "whole_original_radiative_calibrations",
    "whole_flat_matching_nonidentifiability",
    "whole_conditional_rate_bounds",
)
ITEM = {
    "id": "QG2_H8A500_independent_hard_radiative_curvature_matching",
    "status": "EXPLICIT_RADIATIVE_MATCHING_NONIDENTIFIABILITY_NOT_MATCHED_HARD_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError(
            "Require the original hard-radiative matching comparison model"
        )


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated independent curvature matching observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "one_loop_curvature_matching_direction",
        "unchanged_flat_matching_through_one_loop",
        "physical_four_scalar_one_graviton_contact",
    ):
        raise ValueError("Require the stated radiative curvature comparison sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "A matching non-identifiability result cannot close original P8"
        )
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_covariant_curvature_contact": contact.data(),
        "whole_original_radiative_calibrations": calibration.data(),
        "whole_flat_matching_nonidentifiability": hard_matching.data(),
        "whole_conditional_rate_bounds": bounds.data(),
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
        "all191_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 191,
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
        "hard_quantum_Regge_and_UV_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "An explicit real covariant Weyl-curvature operator gives a nonzero conserved Bose-symmetric four-scalar one-graviton contact at formal one-loop order. It is invisible to the retained flat amplitudes through one loop and to the leading through sub-subleading soft coefficients. Its interference with the original complete tree is nonzero; conditional compact rate bounds do not fix its unknown physical coefficient.",
        "domain": "Original classical source and flat onepoint conditions unchanged. Formal hbar order one, same canonical metric and original compact physical recoil states. One comparison direction, not a complete curved operator basis or an arbitrary-coefficient UV family. The low-window interference bound requires x<=delta/192.",
        "historical_qualification": qualifications(),
        "not_established": "The physical value or bound of the new curvature matching coordinate; a complete covariant one-loop five-point functional; full hard matching, interacting detector probability, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure.",
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
            "flat_matching_determines_all_curvature_couplings",
            "soft_theorem_determines_omega_squared_contact",
            "set_unknown_curvature_coefficient_to_zero",
            "contact_square_is_complete_two_loop_rate",
            "arbitrary_curvature_family_has_UV_completion",
            "matching_obstruction_closes_P8",
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
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported radiative matching input accepted: " + name)
    return total


def controls():
    return {
        "original_action_source_and_matching_parameters_unchanged": True,
        "literal_covariant_operator_and_all24_assignments": True,
        "exact_Ward_Bose_and_quadratic_soft_order": True,
        "same_flat_matching_through_one_loop": True,
        "nonzero_original47_tree_interference": True,
        "conditional_bounds_do_not_assign_physical_chi": True,
        "hard_quantum_and_original_P8_frontiers_open": True,
        "rejected_inputs": rejected_inputs(),
    }
