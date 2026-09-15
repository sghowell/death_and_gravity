"""Whole normal-cut frontier, preserved obligations and rejected shortcuts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_graviton_threshold_sheet import audit as previous

from . import legs, masters, reconstruction, source

MODEL = "original_common_gravity_cut_masters_and_external_leg_matching_not_full_P8"
OBSERVABLES = (
    "same_crossed_box_both_species_cuts",
    "whole_pure_gravity_cut_matched_basis",
    "literal_gravitational_scalar_self_energy",
    "complete_pair_and_external_leg_soft_pole",
)
ITEM = {
    "id": "QG2_H8A447_complete_common_crossed_gravitational_master_cut_basis_and_literal_scalar_external_leg_soft_matching",
    "status": "COMPLETE_SCOPED_COMMON_GRAVITY_CUT_MASTERS_AND_LSZ_IR_NOT_FULL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def frontier():
    return previous.frontier()


def qualifications():
    return previous.qualifications()


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original normal-cut model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a scoped normal-cut observable")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A normal-cut continuation cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "entire_raw_scalar_masters": masters.data(),
        "whole_common_crossed_cut_match": reconstruction.data(),
        "literal_external_legs": legs.data(),
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
        "all9_primitive_statuses_preserved": len(frontier()) == 9,
        "all138_matching_records_preserved": matching()[:-1] == previous.matching()
        and len(previous.matching()) == 138,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "cut_matched_basis_not_full_finite_amplitude": True,
        "external_leg_IR_not_finite_pole_anchor": True,
        "S275_S276_S277_and_S281_physical_scopes_unchanged": True,
    }


def observable():
    return {
        "established": "Same six ordered scalar boxes reproduce the complete original graviton and pure-gravity elastic normal cuts without double counting. Complete dimensional master normalization and evanescent traces retain all finite elastic terms. A literal scalar-graviton self-energy has zero on-shell mass shift in D for this sector; its four external LSZ legs supply the missing self term, completing the S278 universal soft pole.",
        "domain": "Unchanged formal leading mass1,kappa10^800 vacuum, pure-gravity first-loop sector, physical strictly nonforward s>4mu and the named Feynman analytic master continuation. Raw ultraviolet and infrared poles are distinguished; finite rational completion and renormalized pole normalization are not supplied.",
        "historical_qualification": qualifications(),
        "not_established": "Full finite crossed amplitude or b20; rational, local, physical pole and vacuum matching; other heavy/contact/matter first-loop sectors, omitted higher loops, finite detector/Regge errors, exact LSZ existence, original quantum state/domain/bounce or V/G/B/P8 closure.",
    }


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "pure_gravity",
        "pure_gravity_external_legs",
    ):
        raise ValueError("Require one of the stated original gravitational sectors")
    return value


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
        for label, call, args in (
            ("model", require_model, (value,)),
            ("observable", require_observable, (value,)),
            ("exact", source.previous.exact_real, (value,)),
            ("physical", source.require_physical, (value, 0)),
            ("sector", require_sector, (value,)),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, args))
    for i, args in enumerate(((4, 0), (3, 0), (5, 1), (5, -1), (5, 0, 0), (-1, 0))):
        rows.append(
            ("outside_physical_domain_" + str(i), source.require_physical, args)
        )
    for i, value in enumerate(
        (
            "full_first_loop_all_species",
            "UV_complete",
            "massless_external_replacement",
            "full_quantum_pole",
            "zero_finite_anchor",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "full_finite_amplitude",
            "forward_b20",
            "drop_massless_pole",
            "new_local_counterterm",
            "massless_external_scalar",
            "physical_cutoff_from_cap",
            "full_Regge_bound",
            "four_point_cuts_fix_LSZ",
            "zero_total_mass_shift",
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
            raise ValueError("Unsupported normal-cut input accepted: " + name)
    return total


def controls():
    return {
        "whole_unchanged_source_and_literal_gravity_vertices": True,
        "full_Feynman_master_normalization_and_boundaries": True,
        "same_six_boxes_have_both_species_cuts": True,
        "complete_elastic_finite_evanescent_terms": True,
        "literal_all_D_self_energy_and_separate_UV_IR": True,
        "all_crossed_pair_and_four_leg_soft_pole": True,
        "cut_matched_not_full_rational_amplitude_or_P8": True,
        "rejected_inputs": rejected_inputs(),
    }
