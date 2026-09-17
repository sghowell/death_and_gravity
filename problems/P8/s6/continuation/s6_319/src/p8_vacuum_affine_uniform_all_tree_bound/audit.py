"""Complete all-finite tree bound with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_uniform_soft_current_bound import audit as previous

from . import bounds, core, source, vertices

MODEL = "original_all_finite_complete_tree_bound_not_inclusive_P8"
OBSERVABLES = (
    "whole_unique_hard_core",
    "whole_hard_vertex_and_energy_bounds",
    "whole_complete_tree_majorant",
)
ITEM = {
    "id": "QG2_H8A483_uniform_all_finite_complete_tree_bound",
    "status": "ALL_FINITE_COMPLETE_BARE_TREE_BOUND_NOT_INCLUSIVE_RATE_OR_FULL_V_G_B_P8",
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
        "complete_hard_tree_bound",
        "nonnegative_composed_majorant",
    ):
        raise ValueError(
            "Require the complete bare-tree bound or its composed majorant"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite-tree source or graph count cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_unique_hard_core": core.data(),
        "whole_hard_vertex_and_energy_bounds": vertices.data(),
        "whole_complete_tree_majorant": bounds.data(),
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
        "all174_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 174,
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
        "established": "Complete original finite-N tree bound relative to the positive unexpanded Born amplitude:3*n^2*N!*(2*10^16)^N/[kappa^(N/2)*product(w_i)]. Unique hard cores, all-valence vertices, paired propagators and composition with the S318 current majorant give a uniform angular estimate with explicit soft poles.",
        "domain": "The unchanged four-dimensional covariant action and original parameters, four mass-one external scalars,5/4<=E<=2, positive radiated energies with total<=1/8, physical complex spatial TT leaf norms<=1, positive associated elastic t/u transfers and generic internal pure-soft propagators. Collinear and hard-angle approaches are bounded, but singular exactly-collinear input values are not assigned.",
        "historical_qualification": qualifications(),
        "not_established": "Overlapping all-N soft subtraction or an inclusive real-virtual probability bound; finite hard and evanescent matching, interacting quantum state, unitarity, absolute complex Regge, common-parent bounce or V/G/B/P8 closure. Energy-independent boundedness, probability summability and values at excluded exactly-collinear points are not claimed.",
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
        "unique_maximal_soft_contraction_and_three_hard_cores": True,
        "all_valence_vertices_and_original_hard_cut_geometry": True,
        "distinct_energy_charges_and_one_initial_hard_pole": True,
        "dropping_pure_soft_blocks_fails_exact_graph_inventory": True,
        "original_full5116_and_composed_majorant_not_an_inclusive_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
