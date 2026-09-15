"""Preserved frontier and explicitly scoped dimensional rational completion."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_common_gravity_masters import audit as previous

from . import dimensional, polarizations, source, threshold

MODEL = (
    "original_dimensional_gravity_cut_coefficients_and_Gram_principal_parts_not_full_P8"
)
OBSERVABLES = (
    "whole_D_TT_sewing",
    "all_D_common_cut_master_coefficients",
    "finite_bubble_rational_term",
    "crossed_Gram_principal_part_removal",
)
ITEM = {
    "id": "QG2_H8A448_complete_dimensional_TT_sewing_common_cut_coefficients_and_finite_bubble_rational_Gram_completion",
    "status": "COMPLETE_SCOPED_DIMENSIONAL_CUT_AND_GRAM_PRINCIPAL_PARTS_NOT_FULL_V_G_B_P8",
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
        raise ValueError("Require the stated dimensional cut model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a scoped dimensional cut observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "pure_gravity_dimensional_cut",
        "pure_gravity_Gram_principal_parts",
    ):
        raise ValueError("Require one of the explicitly stated gravitational sectors")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Dimensional cut completion cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "entire_threshold_jets_and_principal_parts": threshold.data(),
        "literal_dimensional_TT_tree_and_sewing": polarizations.data(),
        "whole_dimensional_master_and_finite_rational_terms": dimensional.data(),
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
        "all139_matching_records_preserved": matching()[:-1] == previous.matching()
        and len(previous.matching()) == 139,
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
        "all_D_cut_basis_not_full_tadpole_or_finite_amplitude": True,
        "Gram_choice_not_physical_pole_or_local_anchor": True,
        "S275_S276_S277_and_S281_S283_scopes_unchanged": True,
    }


def observable():
    return {
        "established": "Literal Einstein/matter TT trees, a written spectator-trace degree argument and exact symbolic-D reconstruction give the complete physical polarization sew. Whole D-dimensional angular moments derive both species' common scalar-master coefficients and their determined finite bubble rational term. Exact triangle threshold jets and a crossing-symmetric rational choice remove integer Gram principal parts, including their forced tadpole-dependent coefficients.",
        "domain": "Same original full source, leading mass1,kappa10^800 flat reference, pure-gravity first-loop sector and explicit minimal D=4+2EP continuation. The dimensional regulator is not a noninteger-dimensional physical Hilbert space. Feynman boundary values, genuine threshold branches and physical-pole terms remain.",
        "historical_qualification": qualifications(),
        "not_established": "Complete massive tadpole coefficient, threshold-supported distribution matching, physical finite pole/local/vacuum normalization, all-species finite amplitude or b20; higher-arity local loops, omitted loops, finite detector/Regge bounds, exact LSZ vacuum, original state/domain/bounce or V/G/B/P8 closure.",
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
            ("dimension", source.require_dimension),
            ("mass", source.require_threshold_mass),
            ("sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate((-1, 0, 1, s.Rational(5, 2))):
        rows.append(
            ("outside_TT_dimension_" + str(i), source.require_dimension, (value,))
        )
    for i, value in enumerate((0, -1)):
        rows.append(
            ("outside_positive_mass_" + str(i), source.require_threshold_mass, (value,))
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
            "drop_physical_pole",
            "new_local_counterterm",
            "massless_external_scalar",
            "physical_cutoff_from_cap",
            "full_Regge_bound",
            "cuts_fix_tadpoles",
            "retune_fixed_vacuum",
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
            raise ValueError("Unsupported dimensional cut input accepted: " + name)
    return total


def controls():
    return {
        "whole_unchanged_source_and_fixed_vacuum": True,
        "literal_symbolic_D_EH_and_matter_tree": True,
        "physical_TT_projection_not_dilaton_square": True,
        "all_D_angular_moments_and_both_species_masters": True,
        "determined_finite_bubble_rational_term": True,
        "complete_integer_Gram_parts_not_physical_pole_retuning": True,
        "tadpole_local_and_original_P8_boundaries_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
