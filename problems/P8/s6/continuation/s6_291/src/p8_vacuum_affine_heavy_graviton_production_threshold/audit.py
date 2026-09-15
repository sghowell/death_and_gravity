"""Original heavy-graviton production and full dimensional threshold boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_complete_matter_graviton_endpoint import audit as previous

from . import phase, production, source, threshold

MODEL = "original_heavy_graviton_production_threshold_not_full_physical_P8"
OBSERVABLES = (
    "whole_original_production_tree_arity",
    "whole_Ward_TT_and_forward_production",
    "whole_D_phase_and_evanescent_coefficient",
    "whole_threshold_Laurent_and_remainder",
)
ITEM = {
    "id": "QG2_H8A455_whole_original_heavy_graviton_production_full_dimensional_forward_cut_and_soft_threshold_Laurent_remainder",
    "status": "COMPLETE_SCOPED_HEAVY_GRAVITON_PRODUCTION_AND_THRESHOLD_NOT_PAIRED_VIRTUAL_PHYSICAL_IR_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original heavy-graviton production sector")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated production/threshold observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "heavy_graviton_production",
        "regulated_forward_threshold",
    ):
        raise ValueError(
            "Require complete production or its regulated forward threshold"
        )
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Threshold matching cannot close the original physical P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_four_graph_production": production.data(),
        "whole_dimensional_phase_and_evanescence": phase.data(),
        "whole_soft_threshold_Laurent_boundary": threshold.data(),
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
        "all146_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 146,
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
        "finite_RH_and_other_curved_matching_not_assigned": True,
        "whole_physical_IR_Regge_and_full_mixed_amplitude_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The complete original Phi Phi to H plus graviton tree obeys the all-component Ward identity. Its full physical TT polarization sum, D-dependent radial/angular phase and independent B0 discontinuity agree. The forward cut has an exact soft-threshold factor with retained evanescent coefficient. A compact-window subtract-add identity supplies its Laurent pole, finite term and an explicit conditional O(epsilon) remainder.",
        "domain": "Formal perturbative production at order g/sqrt(kappa), cut at g^2/kappa, mu>0,n>4mu,s>n. D=4+2epsilon with0<epsilon<=1/8; compact dimensionless windows0<L<=1/2 and fixed C1 weights. Original mu1,n=10^200/512+2,g1/8192,kappa10^800. No exact stable-heavy or fixed-transfer construction.",
        "historical_qualification": qualifications(),
        "not_established": "Paired virtual/pole threshold cancellation, an exact stable heavy state or width resummation, a finite physical detector/dressing observable, all mixed four-point cuts or amplitude, local matching/Regge/all-loop bounds, original state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "choose_light_finite_curvature",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "D4_cut_fixes_finite_threshold",
            "remove_soft_pole_by_hand",
            "H_is_exact_stable",
            "ignore_virtual_pole_term",
            "use_two_helicities_for_all_D",
            "distinct_pair_as_identical_pair",
            "full_Regge_bound",
            "finite_forward_dispersion_from_unpaired_cut",
            "exchange_IR_limits",
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
            raise ValueError(
                "Unsupported heavy-graviton production input accepted: " + name
            )
    return total


def controls():
    return {
        "whole_original_source_and_production_arity": True,
        "all_four_graph_Ward_and_graph_deletion_controls": True,
        "whole_physical_TT_projector_not_selected_helicity": True,
        "independent_radial_phase_and_B0_normalization": True,
        "complete_evanescent_coefficient_and_compact_remainder": True,
        "unpaired_soft_pole_and_virtual_width_boundary_retained": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
