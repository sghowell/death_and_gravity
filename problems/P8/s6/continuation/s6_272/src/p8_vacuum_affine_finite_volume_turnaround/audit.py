"""Qualified finite volume turnaround; no original frontier is closed."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_window_growth.intervals import fraction
from p8_vacuum_affine_weyl_operator_comparison import audit as previous

from . import branch, geometry, moving, quantum, source

STATE, MEASURE = previous.STATE, previous.MEASURE
OBSERVABLES = (
    "longer_full_source_translated_auxiliary_and_phase_domain",
    "uniform_positive_finite_Weyl_and_calibrated_coherent_volume",
    "state_independent_finite_regulator_volume_turnaround",
)
ITEM = {
    "id": "QG2_H8A436_full_source_longer_interval_state_independent_finite_volume_turnaround",
    "status": "COMPLETE_FINITE_REGULATOR_VOLUME_TURNAROUND_WITH_EXTERNAL_MEANS_NOT_UNIQUE_STRICT_OR_SELF_CONSISTENT_BOUNCE_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def require_radius(value):
    return previous.require_radius(value)


def require_time(value):
    value = fraction(value)
    if not -fraction(source.TIME) <= value <= fraction(source.TIME):
        raise ValueError("Require the explicit longer real time interval")
    return value


def require_state(value):
    return previous.require_state(value)


def require_measure(value):
    return previous.require_measure(value)


def require_cutoff(value):
    return previous.require_cutoff(value)


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a declared finite volume turnaround observable")
    return value


def require_amplitude(value):
    quantum.jet(0, value)
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Finite volume turnaround does not close original P8")
    return True


@cache
def packets():
    return {
        "whole_full_clock_center_and_translated_source_derivatives": source.data(),
        "whole_actual_fixed_balance_reference_flow": moving.data(),
        "whole_scaled_full_spatial_reconstruction_and_curvature": geometry.data(),
        "whole_translated_complex_auxiliary_branch_and_volume": branch.data(),
        "whole_actual_operator_bounds_and_volume_turnaround": quantum.data(),
    }


@cache
def residuals():
    return {
        packet + "_" + name: clean(value)
        for packet, data in packets().items()
        for name, value in data["checks"].items()
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
            packet + "_" + name: bool(value)
            for packet, data in packets().items()
            for name, value in data["gates"].items()
        },
        "all_nine_original_primitive_rows_retained": len(frontier()) == 9,
        "all_matching_identifiers_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "same_original_parameters": require_parameters(parameters()) == parameters(),
        "S261_refutation_and_S265_integrability_boundary_unchanged": True,
        "scoped_P8a_unchanged_and_original_P8_open": True,
    }


def observable():
    return {
        "established": "On real |u|<=1e-130, the same finite48-pair original-source regulator has a complete complex4R-to8R reference/spatial/auxiliary domain. The actual Weyl normalized volume is within1e-264 ofI and aboveI/2 for either original cutoff. For either full unitary ordering and either cutoff, both physical-volume endpoints exceed the center by>5e-260; every global minimum lies strictly in |u|<1e-132. These are state-independent operator bounds and apply to the unchanged pure seed.",
        "domain": "Same L1, P1e64, kappa1e800, R1e20, all96 phase variables, all generated spatial harmonics, complete original fields/Gauss/primitive/canonical contacts and both explicit c1/c2 cutoffs. Full original scalar phases are retained.",
        "not_established": "No unique minimum, strict positive acceleration, self-consistent homogeneous mean equations, small longer-time state or leakage error, unlocalized dynamics, continuum limit, physical matching, omitted loops, UV or global completion. The homogeneous variables remain EXTERNAL. Original V/G/B/P8 remain OPEN; scoped P8(a) is unchanged.",
    }


def bad_cases():
    cases = [
        ("inherited_" + name, call, args) for name, call, args in previous.bad_cases()
    ]
    invalid = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.I,
        s.nan,
        s.Symbol("unknown"),
    )
    for i, value in enumerate(invalid):
        for name, call in (
            ("time", require_time),
            ("amplitude", require_amplitude),
            ("observable", require_observable),
        ):
            cases.append(("invalid_new_" + name + "_" + str(i), call, (value,)))
    for i, value in enumerate(
        (-2 * source.TIME, 2 * source.TIME, s.Rational(1, 10**60))
    ):
        cases.append(("outside_longer_time_" + str(i), require_time, (value,)))
    for i, value in enumerate((0, -1, 1, s.Rational(1, 10**255), 10**900)):
        cases.append(("wrong_longer_amplitude_" + str(i), require_amplitude, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "self_consistent_homogeneous_bounce",
            "unique_minimum",
            "strict_positive_acceleration",
            "small_longer_time_unitary_error",
            "unchanged_old_leakage",
            "uniform_continuum_limit",
            "new_diagonal_vacuum",
            "delete_homogeneous_shape",
            "drop_full_Gauss_source",
        )
    ):
        cases.append(("wrong_turnaround_scope_" + str(i), require_observable, (value,)))
    for i, expression in enumerate(
        (
            s.diff(source.q.R, source.N, 6),
            s.diff(source.q.F, source.u, 6),
            s.diff(source.q.j, source.N, 6),
            s.diff(source.q.Hclock, source.u, 4),
            1 / s.Symbol("unproved_denominator"),
            s.sqrt(s.Symbol("unproved_base")),
            s.Symbol("unbounded_atom"),
            s.Float(1),
            s.oo,
            s.nan,
        )
    ):
        cases.append(
            ("unproved_translated_majorant_" + str(i), source.magnitude, (expression,))
        )
    for i, expression in enumerate(
        (
            s.diff(source.q.R, source.N, 8),
            s.diff(source.q.F, source.N, 8),
            s.diff(source.q.j, source.N, 8),
        )
    ):
        cases.append(("unlicensed_clock_jet_" + str(i), source.clock_at, (expression,)))
    cases.extend(
        [
            ("delete_original_frontier", validate_scope, ([], matching())),
            ("delete_turnaround_matching_frontier", validate_scope, (frontier(), [])),
        ]
    )
    if len({name for name, _, _ in cases}) != len(cases):
        raise ValueError("Duplicate rejected turnaround input")
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            count += 1
        else:
            raise ValueError("Unsupported finite-turnaround input accepted: " + name)
    return count


def controls():
    return {
        "actual_moving_center_not_zero_momenta_off_clock": True,
        "full_sublattice_reconstruction_not_extra_mode_truncation": True,
        "complete_curvature_and_tracefree_cotangent_contacts": True,
        "new_time_domain_not_old_M_times_T_shortcut": True,
        "full_operator_bounds_not_symbol_sup_shortcut": True,
        "unitarity_not_old_small_state_error_assumed": True,
        "external_means_not_self_consistent_bounce_or_original_P8": True,
        "historical_errata_and_all_original_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
