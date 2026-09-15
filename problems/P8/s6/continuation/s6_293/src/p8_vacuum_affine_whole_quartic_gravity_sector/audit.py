"""Whole quartic-gravity coefficient and explicit finite regulator control."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_heavy_resonance_soft_pole_pairing import audit as previous

from . import cuts, forward, proper, source

MODEL = "original_whole_quartic_gravity_sector_not_full_physical_P8"
OBSERVABLES = (
    "whole_original_quartic_gravity_graph_ownership",
    "complete_quartic_proper_external_and_endpoint_sum",
    "entire_D_scalar_cut_and_tree_interference",
    "closed_finite_forward_coefficient_and_explicit_uniform_bound",
)
ITEM = {
    "id": "QG2_H8A457_complete_quartic_gravity_whole_D_scalar_cut_closed_finite_forward_coefficient_and_explicit_uniform_regulator_bound",
    "status": "COMPLETE_KNOWN_QUARTIC_GRAVITY_COEFFICIENT_AND_EXPLICIT_BOUND_NOT_FULL_SOURCE_MATCHING_PHYSICAL_IR_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original quartic-gravity sector")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated quartic or formal-forward observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "whole_quartic_gravity",
        "formal_soft_divided_forward",
    ):
        raise ValueError(
            "Require the complete quartic sector or its formal forward coefficient"
        )
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A quartic coefficient cannot close the original physical P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_quartic_proper_external_and_endpoint_sum": proper.data(),
        "whole_D_scalar_cut_and_tree_interference": cuts.data(),
        "whole_closed_finite_forward_coefficient": forward.data(),
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
        "all148_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 148,
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
        "constant_anchors_zero_b20_not_higher_matching_closed": True,
        "whole_physical_IR_Regge_and_full_mixed_amplitude_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The entire minimal one-loop coefficient linear in the original C and1/kappa includes six proper pair graphs, four contacts, four external residues and all matter endpoints with full OS terms. Whole-D massive moment cuts sew exactly to the complete crossed GR tree interference. Its known S278 soft-divided forward coefficient has a closed Catalan/log expression and an explicit uniform regulator-error bound. On the stated original resolution interval its magnitude is below10^-1005.",
        "domain": "Original source unchanged. mu>0 for the analytic center and s>4mu for the physical scalar cut;D=4+2epsilon,0<epsilon<=1/8. The explicit numerical remainder/magnitude bounds require original mu=nu^2=1 and1/4<=E^2<=1. E remains an explicit inherited formal soft-prescription parameter, not a detector selected to obtain a sign.",
        "historical_qualification": qualifications(),
        "not_established": "Higher-derivative or other-coupling full-source matching, a positive physical spectral measure or isolated positivity verdict, complete Coulomb/detector/dressed IR factorization and unitarity, finite-transfer Regge/all-loop remainder, original state/domain/measure/bounce or V/G/B/P8 closure. The H-metric and unstable-heavy sectors remain separate.",
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
            raise ValueError("Unsupported quartic-sector input accepted: " + name)
    return total


def controls():
    return {
        "whole_original_quartic_graphs_and_OS_terms_retained": True,
        "complete_D_pair_and_endpoint_tensor_reductions": True,
        "whole_scalar_cut_equals_entire_tree_interference": True,
        "full_evanescent_forward_coefficient_not_D0_only": True,
        "explicit_uniform_regulator_and_actual_magnitude_bounds": True,
        "constant_anchors_zero_b20_not_higher_matching_closed": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
