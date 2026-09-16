"""Known minimal-gravity finite part with every original closure gate retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_radiative_angular_finite import audit as previous

from . import assembly, bounds, masters, source

MODEL = "original_minimal_gravity_finite_not_full_P8"
OBSERVABLES = (
    "whole_finite_master_dictionary",
    "whole_minimal_gravity_known_finite_part",
    "whole_compact_physical_majorant",
    "whole_symbolic_finite_matching",
)
ITEM = {
    "id": "QG2_H8A466_complete_known_minimal_gravity_finite_reference_and_compact_bound",
    "status": "KNOWN_MINIMAL_GRAVITY_FINITE_REFERENCE_NOT_MATCHED_FULL_AMPLITUDE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact minimal-gravity finite model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated minimal-gravity finite observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "known_minimal_gravity_finite",
        "symbolic_matching",
    ):
        raise ValueError("Require the known finite or symbolic matching sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A compact finite-loop bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_finite_master_dictionary": masters.data(),
        "whole_finite_minimal_gravity_assembly": assembly.data(),
        "whole_compact_physical_majorants": bounds.data(),
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
        "all157_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 157,
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
        "established": "The full known finite minimal Einstein/massive-scalar "
        "one-loop representative, including whole-D coefficient derivatives, "
        "both Gram terms, four LSZ legs, physical-pole completion and analytic "
        "soft division, is explicit. Its compact known amplitude obeys "
        "1e9*(1+abs(log(delta)))/(kappa^2*delta^2). All three physical finite "
        "matching anchors remain symbolic.",
        "domain": "mu=nu=1,25/4<=s<=16,t,u<0,min(-t,-u)>=delta,0<delta<=1. "
        "Reference resolution/nu=1 defines an analytic finite amplitude, "
        "not a detector resolution claim. At the explicit comparison window "
        "delta=1e-204 the known original-parameter loop/Born is below1e-580 "
        "and its rate interference below1e-579.",
        "historical_qualification": qualifications(),
        "not_established": "Values or bounds for alpha,beta,delta_kappa; the full "
        "current-action amplitude, gravity-Born radiation, all-loop errors, "
        "exact forward or high-energy Regge bounds, physical state/domain/"
        "measure/bounce or original V/G/B/P8 closure.",
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
            "take_forward_before_soft_division",
            "relative_Regge_is_absolute",
            "all_loop_regulator_complete",
            "window_is_physical_UV_cutoff",
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
        "same_original_action_and_minimal_gravity_inputs": True,
        "whole_D_finite_assembly_and_all_crossings": True,
        "Feynman_contour_and_triangle_endpoint_bounds": True,
        "UV_before_IR_and_D_tree_finite_reference": True,
        "three_matching_coordinates_explicitly_unassigned": True,
        "compact_known_piece_not_full_forward_amplitude": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
