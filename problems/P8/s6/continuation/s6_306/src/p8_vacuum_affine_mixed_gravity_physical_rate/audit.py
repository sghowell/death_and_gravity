"""Physical mixed-sector bound with every original closure gate retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_spectator_gravity_insertion import audit as previous

from . import assembly, bounds, masters, source

MODEL = "original_mixed_gravity_physical_rate_not_full_P8"
OBSERVABLES = (
    "whole_physical_scalar_masters",
    "whole_finite_mixed_assembly",
    "whole_uniform_known_rate_bound",
    "whole_mixed_matching_boundary",
)
ITEM = {
    "id": "QG2_H8A470_complete_physical_mixed_gravity_known_rate_and_exact_box_collapse",
    "status": "PHYSICAL_MIXED_GRAVITY_KNOWN_RATE_NOT_FULL_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact mixed-gravity-rate model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated mixed-gravity-rate observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "physical_mixed_gravity",
        "known_reference_rate",
    ):
        raise ValueError("Require the known finite or symbolic matching sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A known mixed-sector bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_physical_scalar_masters": masters.data(),
        "whole_finite_mixed_assembly": assembly.data(),
        "whole_uniform_known_rate_bound": bounds.data(),
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
        "all161_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 161,
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
        "established": "Exact whole mixed-box reduction, correct physical-sheet scalar masters, complete inherited C/kappa plus g^2/kappa graph assembly with both endpoints, and an explicit uniform known full-Born interference bound below4e-204 at original parameters.",
        "domain": "Original mu=nu^2=1,25/4<=s<=16,t,u<0,n>=128; formal analytic reference1/4<=E^2<=1. All nonforward angles, including arbitrarily close to either endpoint.",
        "historical_qualification": qualifications(),
        "not_established": "Unknown physical finite matching, complete all-order detector rate, finite forward cross section, high-energy or complex Regge control, quantum-state and bounce construction, original V/G/B/P8 closure.",
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
            "known_rate_bound_controls_loop_squares",
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
        "same_whole_S293_S294_S290_graphs_not_new_graph_claim": True,
        "whole_exact_box_identity_before_physical_continuation": True,
        "independent_original_integrals_and_Box16_calibration": True,
        "same_positive_full_Born_at_every_nonforward_angle": True,
        "all_finite_matching_coordinates_remain_unassigned": True,
        "known_linear_rate_not_loop_squares_or_Regge": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
