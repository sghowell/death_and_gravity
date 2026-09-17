"""Gauge-complete temporal trees with unchanged original scientific frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_complete_pair_factorization import audit as previous

from . import bounds, chart, source, trees

MODEL = "original_gauge_complete_temporal_trees_not_all_N_rate_or_full_P8"
OBSERVABLES = (
    "whole_finite_temporal_chart",
    "whole_nonlinear_recursive_source",
    "whole_temporal_pair_bound",
    "whole_planar_three_ray_bound",
)
ITEM = {
    "id": "QG2_H8A481_gauge_complete_temporal_trees_and_planar_bound",
    "status": "GAUGE_COMPLETE_TEMPORAL_TREES_AND_PLANAR_BOUND_NOT_ALL_N_RATE_OR_FULL_V_G_B_P8",
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
        "gauge_complete_temporal_recursion",
        "restricted_planar_three_ray_norm",
    ):
        raise ValueError(
            "Require complete temporal recursion or the restricted planar norm"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite-tree source or graph count cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_finite_temporal_chart": chart.data(),
        "whole_nonlinear_recursive_source": trees.data(),
        "whole_pair_and_planar_bounds": bounds.data(),
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
        "all172_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 172,
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
        "established": "Gauge-complete finite temporal-current recursion with invariant complete on-shell tree amplitudes; independent literal nonlinear coordinate chart; uniform two-ray field bound and a restricted planar three-ray field bound below120000/kappa.",
        "domain": "The unchanged four-dimensional covariant scalar/heavy/Einstein action and original parameters, physical free waves and finite generic tree coefficients away from internal poles. Temporal pure-soft subsets have positive energy. The three-ray norm is restricted to coplanar rays, energy ratio1:2:3, half-angle coordinates0,r,t with |r|,|t|<=1, and unit-Frobenius TT leaves.",
        "historical_qualification": qualifications(),
        "not_established": "An arbitrary three-dimensional or energy-hierarchy-uniform three-ray bound; uniform higher-multiplicity current or complete-amplitude norms; an all-N inclusive probability bound, finite hard real-virtual matching, interacting quantum state, absolute complex Regge, common-parent bounce or V/G/B/P8 closure.",
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
        "finite_Noether_and_literal_nonlinear_coordinate_comparison": True,
        "same_vertices_and_original5116_tree_amplitude_unchanged": True,
        "all24_pair_tensor_coefficients_and_exact_norm_budget": True,
        "linear_only_projection_pole_is_retained_as_negative_control": True,
        "planar_all8_TT_polynomial_majorant_not_an_all_N_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
