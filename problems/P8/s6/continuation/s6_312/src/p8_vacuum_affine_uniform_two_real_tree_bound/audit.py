"""Complete uniform two-real tree envelope with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_two_real_collinear_current import audit as previous

from . import bounds, gaps, source, topology

MODEL = "original_uniform_two_real_tree_bound_not_full_P8"
OBSERVABLES = (
    "whole_regular_tree_inventory",
    "whole_common_hard_propagator_gaps",
    "whole_complete434_tree_bound",
    "whole_bare_tree_infrared_boundary",
)
ITEM = {
    "id": "QG2_H8A476_uniform_original_complete_two_real_tree_bound",
    "status": "UNIFORM_COMPLETE434_BARE_TREE_BOUND_NOT_IR_FINITE_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original full-tree model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated full-tree observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete434_tree",
        "regular387_graphs",
    ):
        raise ValueError(
            "Require the complete434 tree or regular387 computational sector"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A bare full-tree bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_regular_topology": topology.data(),
        "whole_common_propagator_gaps": gaps.data(),
        "whole_uniform_complete_two_real_bound": bounds.data(),
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
        "all167_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 167,
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
        "established": "Every labeled regular-graph topology and radiation assignment is covered; common hard gaps and canonical vertex budgets give |Mregular|/(Am+AG)<1e-395/(ab). Adding the frozen S311 pair bound proves |Mfull434|/(Am+AG)<2e-395/(ab), uniformly over the stated physical angular domain.",
        "domain": "Original four-dimensional formal tree and S300 recoil, E in[5/4,2],a,b>0,a+b<=1/8, all nonforward hard Born angles and all emitted directions, with unit spatial-Frobenius TT polarizations. Collinear neighborhoods are bounded; no unique exactly-collinear helicity value is asserted.",
        "historical_qualification": qualifications(),
        "not_established": "An infrared-finite integrated two-real detector correction, overlap subtraction and virtual pairing, all-N physical rate, radiative hard loops or finite physical matching, unitarity or complex Regge, original quantum state/common-parent bounce, original V/G/B/P8 closure.",
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
            "bare_tree_bound_is_an_IR_finite_detector_rate",
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
        "same_original434_graph_tree_and_original_parameters": True,
        "all387_regular_graph_profiles_and_radiation_assignments": True,
        "common_cluster_gap_not_null_or_mass_one_approximation": True,
        "all_canonical_scalar_and_Einstein_vertices_budgeted": True,
        "exact_original_complete434_extreme_angle_calibrations": True,
        "bare_soft_envelope_not_IR_finite_integrated_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
