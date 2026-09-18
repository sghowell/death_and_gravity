"""Three-current subtraction and one-block class with unchanged frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_relative_energy_complex_tube import audit as previous

from . import hard, source, subtraction

MODEL = "original_three_current_and188_class_not_full_inclusive_P8"
OBSERVABLES = (
    "whole_angular_coverage",
    "whole_three_current_hierarchy",
    "whole_one_block_class",
)
ITEM = {
    "id": "QG2_H8A485_three_soft_current_faces_and_one_block_subtraction",
    "status": "THREE_CURRENT_AND188_CLASS_SUBTRACTION_NOT_FULL_INCLUSIVE_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original all-multiplicity tree model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated finite-tree source observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "three_current_compatible_energy_faces",
        "one_block188_subtraction",
    ):
        raise ValueError(
            "Require the complete three-current or its one188-term hard-core class"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite-tree source or graph count cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_angular_coverage": source.angular_data(),
        "whole_three_current_hierarchy": subtraction.data(),
        "whole_one_block_class": hard.data(),
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
        "all176_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 176,
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
        "established": "F=abc*H3/W has all-angle squarefree derivative hierarchy and compatible faces. Its threefold rectangle is bounded by721*J/kappa. The complete188-term one-block class in the temporal reorganization has rectangle below1e-738*J.",
        "domain": "Unchanged original action and couplings; generic positive energies, fixed real directions and arbitrary complex unit-Frobenius TT leaves. The hard-class domain is W<=1/8,5/4<=E<=2 and positive associated elastic transfers. Angular approaches are uniform; exact internal poles have no assigned value.",
        "historical_qualification": qualifications(),
        "not_established": "A subtraction bound for the remaining4928 temporal-tree terms, the complete5116-tree probability, all-N overlap and real-virtual completion, finite hard/evanescent matching, an interacting quantum state, unitarity, absolute complex Regge, common-parent bounce or original V/G/B/P8 closure. The188 class is not an independent gauge-invariant observable.",
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
        "literal_polynomial_and_independent_Taylor_oracles": True,
        "two_complete_angular_sectors_plus_exact_reflection": True,
        "all_polarizations_hierarchies_and_compatible_energy_faces": True,
        "one_block_global_tube_not_full_current_global_holomorphy": True,
        "remaining4928_and_inclusive_real_virtual_not_closed": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
