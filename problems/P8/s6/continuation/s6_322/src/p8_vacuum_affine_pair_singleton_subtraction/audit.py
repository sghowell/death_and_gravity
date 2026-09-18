"""Pair-plus-singleton subtraction with all inherited scientific boundaries."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_three_soft_current_subtraction import audit as previous

from . import algebra, hard, hierarchy, source

MODEL = "original_three_pair_classes_and1349_subtraction_not_inclusive_P8"
OBSERVABLES = (
    "whole_pair_hierarchy",
    "whole_offshell_core",
    "whole_pair_class_subtraction",
)
ITEM = {
    "id": "QG2_H8A486_pair_singleton_subtraction_and_all_nonsingleton_faces",
    "status": "ALL1349_NON_SINGLETON_SUBTRACTION_NOT_FULL5116_INCLUSIVE_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original pair-plus-singleton tree model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated partial tree-subtraction observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "pair387_subtraction",
        "all1349_nonsingleton_subtraction",
    ):
        raise ValueError("Require one387 class or all1349 nonsingleton terms")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A partial tree subtraction cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_pair_hierarchy": hierarchy.data(),
        "whole_offshell_core": algebra.data(),
        "whole_pair_class_subtraction": hard.data(),
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
        "all177_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 177,
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
        "established": "The complete temporal pair has uniform lower energy derivatives and compatible faces. Each387 pair-plus-singleton class has threefold rectangle below1e-723*c*I(a,b), with permutations for the other pairs. Together with S321, all1349 nonsingleton terms have rectangle below2e-723*J.",
        "domain": "Unchanged original action and couplings; W<=1/8,5/4<=E<=2, positive associated Born transfers, real null directions and unit complex TT external leaves. Composite pair momentum may be timelike and its spatial tensor traceful and nontransverse. Exact original internal poles are excluded; angular approaches are uniform.",
        "historical_qualification": qualifications(),
        "not_established": "The remaining3767 three-singleton terms or the complete5116-tree probability subtraction; all-N overlap, real-virtual completion, finite hard and evanescent matching; an interacting quantum state, unitarity, absolute complex Regge, common-parent bounce or original V/G/B/P8 closure. Partial gauge-fixed classes are not independent observables.",
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
                "unsupported_full3767_subtraction",
                require_sector,
                ("full3767_already_subtracted",),
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
        "all24_complete_pair_coefficients_and96_hierarchy_certificates": True,
        "off_shell_trace_and_nontransverse_corrections_retained": True,
        "independent_original387_140_and_temporal_class_calibrations": True,
        "anisotropic_complex_points_outside_old_relative_c_tube": True,
        "remaining3767_and_full5116_probability_not_closed": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
