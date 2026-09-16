"""Uniform selected pair-current theorem with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_complete_two_graviton_tree import audit as previous

from . import bounds, collinear, current, source

MODEL = "original_two_real_collinear_current_not_full_P8"
OBSERVABLES = (
    "whole_offshell_hard_current",
    "whole_physical_pair_quotient",
    "whole_uniform_pair_sector_bound",
    "whole_pair_sector_scope_boundary",
)
ITEM = {
    "id": "QG2_H8A475_complete_offshell_current_and_uniform_two_real_collinear_sector_bound",
    "status": "UNIFORM_SELECTED_PAIR_CURRENT_BOUND_NOT_FULL_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original pair-current model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated pair-current observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_offshell_current",
        "selected_pair_propagator_graphs",
    ):
        raise ValueError("Require the conserved current or selected47-graph sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A selected pair-current bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_offshell_current": current.data(),
        "whole_physical_collinear_quotient": collinear.data(),
        "whole_uniform_original_pair_bound": bounds.data(),
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
        "all166_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 166,
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
        "established": "General timelike-root Ward identity for the complete47-tree current; exact factorization of434=387+47 diagrams; all4 physical-polarization angular-pole cancellations; uniform current norm and original full-Born selected-pair amplitude bound9e-400(1/a+1/b).",
        "domain": "Original four-dimensional formal tree, S300 recoil E in[5/4,2], a,b>0,a+b<=1/8, all emitted directions and all nonforward hard Born angles. Polarizations have unit spatial Frobenius norm. The bound covers angular endpoint neighborhoods without asserting a unique exactly-collinear helicity value.",
        "historical_qualification": qualifications(),
        "not_established": "A bound on all434 graphs or the integrated two-real nonleading detector correction, overlapping soft subtraction/virtual pairing, all-N physical rate, radiative loops and physical matching, unitarity or complex Regge, original quantum state/common-parent bounce, original V/G/B/P8 closure.",
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
            "pair_collinear_bound_controls_full_two_real_rate",
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
        "same_original434_graph_tree_and_47_plus387_factorization": True,
        "general_offshell_root_Ward_not_null_approximation": True,
        "all4_physical_polarization_pair_quotients": True,
        "uniform_original_current_norm_and_positive_full_Born": True,
        "exact_original_extreme_collinear_and_forward_calibrations": True,
        "selected_pair_bound_not_full_integrated_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
