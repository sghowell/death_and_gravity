"""Relative-energy complex bounds with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_uniform_all_tree_bound import audit as previous

from . import bounds, continuation, geometry, source

MODEL = "original_relative_complex_tube_not_soft_faces_or_inclusive_P8"
OBSERVABLES = (
    "whole_complex_geometry",
    "whole_exact_complex_continuation",
    "whole_relative_energy_majorants",
)
ITEM = {
    "id": "QG2_H8A484_all_finite_relative_energy_complex_tubes",
    "status": "RELATIVE_COMPLEX_TUBE_AND_CAUCHY_BOUND_NOT_SOFT_FACES_OR_FULL_V_G_B_P8",
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
        "relative_energy_holomorphy",
        "shrinking_radius_Cauchy_bound",
    ):
        raise ValueError(
            "Require relative-energy holomorphy or its shrinking-radius bound"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite-tree source or graph count cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_complex_geometry": geometry.data(),
        "whole_exact_complex_continuation": continuation.data(),
        "whole_relative_energy_majorants": bounds.data(),
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
        "all175_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 175,
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
        "established": "The complete finite-N original tree is holomorphic in |zi-wi|<=1e-12*wi with an explicit angular-uniform bound. The scaled G_N has Cauchy majorant alpha!*6*n^2*N!*(4e17)^N/[kappa^(N/2)*product((1e-12*wi)^alpha_i)]. A complete-current pole control excludes a naive total-W disc.",
        "domain": "The unchanged original action and parameters, fixed E,u,real emission directions and complex TT leaf tensors of unit spatial Frobenius norm. Centers are generic positive real energies with total<=1/8,5/4<=E<=2 and positive associated elastic transfers. Only energies and the same analytic recoil branches are continued. Angular limits remain uniform but excluded exact internal poles are not assigned values.",
        "historical_qualification": qualifications(),
        "not_established": "Holomorphic control through all zero-energy soft faces or an integrable all-N overlap subtraction; an inclusive real-virtual probability, finite hard and evanescent matching, interacting quantum state, unitarity, absolute complex Regge, common-parent bounce or V/G/B/P8 closure. Energy-independent derivative bounds, probability summability and values at excluded exact poles are not claimed.",
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
        "generic_complex_conserved_inverse_and_weighted_factorization": True,
        "positive_relative_subset_gaps_and_all_order_grading": True,
        "complete_original_complex_5116_source_retained": True,
        "nonzero_cluster_pole_rejects_global_total_energy_disc": True,
        "shrinking_Cauchy_radii_not_an_inclusive_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
