"""Complete known mixed heavy-gravity coefficient, not original physical closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_whole_quartic_gravity_sector import audit as previous

from . import cuts, forward, graphs, source

MODEL = "original_whole_mixed_heavy_gravity_sector_not_full_physical_P8"
OBSERVABLES = (
    "whole_original_mixed_heavy_gravity_graph_ownership",
    "complete_offshell_vertex_selfenergy_contact_and_box_sum",
    "entire_mixed_two_body_cuts_and_covariant_resonance_matching",
    "whole_finite_forward_coefficient_and_explicit_known_bound",
)
ITEM = {
    "id": "QG2_H8A458_complete_mixed_heavy_gravity_whole_offshell_amplitude_all_two_body_cuts_finite_forward_and_explicit_known_bound",
    "status": "COMPLETE_KNOWN_MINIMAL_MIXED_GRAVITY_COEFFICIENT_AND_EXPLICIT_BOUND_NOT_FULL_MATCHING_PHYSICAL_IR_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original mixed heavy-gravity sector")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated complete mixed-sector observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "whole_mixed_heavy_gravity",
        "formal_soft_divided_forward",
    ):
        raise ValueError(
            "Require the whole mixed sector or its formal forward coefficient"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A known mixed coefficient cannot close original physical P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_offshell_proper_selfenergy_contact_and_box_graphs": graphs.data(),
        "whole_mixed_scalar_heavy_and_graviton_cuts": cuts.data(),
        "whole_finite_forward_representative_and_known_bound": forward.data(),
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
        "all149_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 149,
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
        "known_representative_bound_not_unmatched_RH_or_higher_matching": True,
        "whole_physical_IR_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The complete minimal one-loop g^2/kappa graph sum has full-D offshell tensor and scalar-master reductions, all twelve mixed boxes, the two metric-cubic contact bubble, heavy selfenergy and four external residues, and the whole inherited matter-metric endpoint. All light/heavy/Hh two-body cuts and the finite H-metric resonance term are retained. Exact infrared subtraction gives a finite subthreshold representative and an explicit known forward bound below10^-1001 at original parameters.",
        "domain": "Original source unchanged. General formulas use mu>0,n>4mu and the retained Feynman boundary. Positive subthreshold masters support the forward coefficient. The numeric majorant requires mu=1,n>=8,|v|<=1/2 and1/4<=E^2<=1; the final exponent uses the exact original n,g,kappa. The S278 analytic soft-division premise remains conditional.",
        "historical_qualification": qualifications(),
        "not_established": "Independent finite R H/heavy-residue and higher-EFT matching, other coupling and all-loop errors, a stable heavy resonance or physical detector/dressed infrared measure, fixed-transfer Regge remainder, original state/domain/measure/bounce or V/G/B/P8 closure. Passing the known loop bound does not prove positivity or a UV completion.",
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
            raise ValueError("Unsupported mixed-sector input accepted: " + name)
    return total


def controls():
    return {
        "whole_original_mixed_graphs_OS_and_onepoint_terms_retained": True,
        "complete_D_offshell_tensor_and_twelve_box_reductions": True,
        "all_three_two_body_channels_and_full_RH_covariance": True,
        "whole_Gamma_evanescent_and_negative_finite_box_term": True,
        "explicit_analytic_original_parameter_magnitude_bound": True,
        "independent_RH_and_higher_matching_not_bounded_by_representative": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
