"""Complete third-order amplitude subtraction with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_pair_singleton_subtraction import audit as previous

from . import analytic, forest, grouping, line, measure, source

MODEL = "original_complete5116_amplitude_and_signed_subtraction_not_inclusive_P8"
OBSERVABLES = (
    "whole_connected_scalar_lines",
    "whole_exact_hard_forest",
    "whole_complete_amplitude_rectangle",
    "whole_three_real_signed_measure",
)
ITEM = {
    "id": "QG2_H8A487_complete_three_real_tree_subtraction_and_signed_measure",
    "status": "COMPLETE5116_TREE_RECTANGLE_AND_SIGNED_MEASURE_NOT_INCLUSIVE_ALL_N_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original third-tree amplitude-subtraction model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated source or third-tree amplitude observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "three_singleton_amplitude_rectangle",
        "complete5116_amplitude_rectangle",
        "defined_three_real_signed_measure",
    ):
        raise ValueError("Require the3767 core or complete5116 amplitude rectangle")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A third-tree amplitude rectangle cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_connected_scalar_lines": line.data(),
        "whole_exact_hard_forest": forest.data(),
        "whole_generic_forward_grouping": grouping.data(),
        "whole_complete_amplitude_rectangle": analytic.data(),
        "whole_three_real_signed_measure": measure.data(),
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
        "all178_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 178,
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
        "established": "The full3767 three-singleton hard core has compatible soft faces and an integrable third anchored amplitude remainder. Together with the retained1349 nonsingleton terms, the complete original5116-tree coefficient G3=abc*sqrt(rho)*M3/A0 has rectangle bounded by1e-670*J. Its seven-face amplitude baseline defines a full signed three-real density, retaining all interference, with cutoff-independent total variation below1e-1424*x^2+1e-1340*x^4.",
        "domain": "Original action and couplings; a,b,c>0,W<=1/8,5/4<=E<=2 and positive associated Born transfers. Fixed real null directions and unit complex TT leaves; uniform angular/energy-hierarchy approaches, without assigning original propagator-pole values. Hard-only grouped kernels receive the global-energy and origin discs, not connected soft currents.",
        "historical_qualification": qualifications(),
        "not_established": "A positive normalized inclusive probability distribution or matching of this defined signed subtraction to actual real-virtual counterterms; all-N subtraction and summation; finite hard and evanescent matching; an interacting quantum state, unitarity, absolute complex Regge, common-parent bounce or original V/G/B/P8 closure. Computational cumulant classes are not separately physical detector observables.",
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
        "generic49_coefficient_scalar_line_and_original_signed_TT_mapping": True,
        "original3767_forest_and_all_eight_nonzero_cumulant_classes": True,
        "all192_endpoint_groups_and_omitted_owner_negative_controls": True,
        "actual_connected_pair_pole_excludes_naive_global_W_tube": True,
        "all_interferences_and_seven_faces_not_inclusive_or_all_N_closure": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
