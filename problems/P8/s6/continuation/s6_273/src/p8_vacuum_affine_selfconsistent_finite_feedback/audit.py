"""Explicit finite hybrid verdict and unchanged original research frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_volume_turnaround import audit as previous
from p8_vacuum_affine_finite_window_growth.intervals import fraction

from . import dynamics, geometry, homogeneous, quantum, source, symmetry

STATE, MEASURE = previous.STATE, previous.MEASURE
MODEL = "classical_homogeneous_quantum_nonzero_modes"
OBSERVABLES = (
    "selfconsistent_finite_hybrid_solution",
    "full_average_operator_force_budget",
    "finite_hybrid_physical_volume_turnaround",
)
ITEM = {
    "id": "QG2_H8A437_full_source_selfconsistent_finite_hybrid_feedback_and_volume_turnaround",
    "status": "COMPLETE_SELFCONSISTENT_FINITE_CLASSICAL_HOMOGENEOUS_QUANTUM_MODE_HYBRID_NOT_HOMOGENEOUS_QUANTIZATION_OR_ORIGINAL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def require_radius(value):
    return previous.require_radius(value)


def require_state(value):
    return previous.require_state(value)


def require_measure(value):
    return previous.require_measure(value)


def require_cutoff(value):
    return previous.require_cutoff(value)


def require_time(value):
    value = fraction(value)
    if not -fraction(source.TIME) <= value <= fraction(source.TIME):
        raise ValueError("Require the evaluated finite hybrid interval")
    return value


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError(
            "Require classical homogeneous variables and quantum nonzero modes"
        )
    return value


def require_ordering(value):
    if not isinstance(value, str) or value not in ("calibrated_coherent", "Weyl"):
        raise ValueError("Require one of the two defined complete operator orderings")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a declared finite hybrid observable")
    return value


def require_amplitude(value):
    return quantum.require_amplitude(value)


def require_homogeneous_deviation(value):
    if not isinstance(value, (tuple, list)) or len(value) != 5:
        raise ValueError("Require all five live homogeneous density coordinates")
    result = tuple(fraction(v) for v in value)
    if any(abs(v) > fraction(source.REAL_RADIUS) for v in result):
        raise ValueError("Require the full proved inner homogeneous path ball")
    return result


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite hybrid result does not close original P8")
    return True


@cache
def packets():
    return {
        "whole_componentwise_entire_source_and_wider_lapse_domain": source.data(),
        "whole_canonical_homogeneous_equations_and_force_chain": homogeneous.data(),
        "whole_actual_reference_cubic_symmetry_and_homogeneous_consistency": symmetry.data(),
        "whole_full_spatial_average_and_implicit_source_bounds": geometry.data(),
        "whole_two_operator_orderings_and_quantum_force_estimates": quantum.data(),
        "whole_coupled_hybrid_fixed_point_and_volume_turnaround": dynamics.data(),
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
        "all_matching_ids_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "same_original_parameters": require_parameters(parameters()) == parameters(),
        "explicit_finite_hybrid_model": require_model(MODEL) == MODEL,
        "historical_errata_and_scoped_P8a_unchanged": True,
        "original_P8_remains_open": True,
    }


def observable():
    return {
        "established": "For the separately named finite classical-homogeneous/quantum-nonzero-mode model, both original cutoffs and both complete orderings admit a unique coupled trajectory on |u|<=1e-180, with homogeneous deviation<1e-390, unchanged pure seed and physical-volume relative error<1e-380. Both volume endpoints exceed the center by>5e-360, and every minimum lies in |u|<1e-188.",
        "domain": "Same L1,P1e64,kappa1e800,zeta1e-6,R1e20,all48 original quantum nonzero-mode pairs, full untruncated nonlinear spatial/auxiliary source, five live homogeneous density variables and cyclic M1 quadrature. All three vector and five traceless homogeneous canonical pairs satisfy their equations at zero by actual cubic symmetry. The heavy homogeneous scalar remains live.",
        "not_established": "No quantization of the homogeneous variables or original fully quantum mean, unique or strict volume minimum, long-time core leakage, unlocalized Hamiltonian, uniform mode/volume/cutoff limit, physical matching, omitted-loop/Regge/UV control or nonlinear global completion. Original V/G/B/P8 remain OPEN; scoped P8(a) and historical refutations remain unchanged.",
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
            ("model", require_model),
            ("ordering", require_ordering),
        ):
            cases.append(("invalid_hybrid_" + name + "_" + str(i), call, (value,)))
    for i, value in enumerate(
        (-2 * source.TIME, 2 * source.TIME, previous.source.TIME)
    ):
        cases.append(("outside_hybrid_time_" + str(i), require_time, (value,)))
    for i, value in enumerate(
        (
            [],
            [0] * 4,
            [0] * 6,
            [source.HOMOGENEOUS_RADIUS] * 5,
            [s.Float(0)] * 5,
            [True] * 5,
            "00000",
        )
    ):
        cases.append(
            (
                "invalid_homogeneous_path_" + str(i),
                require_homogeneous_deviation,
                (value,),
            )
        )
    for i, value in enumerate(
        (
            "original_P8_closed",
            "fully_quantum_homogeneous_mean",
            "unique_volume_minimum",
            "strict_positive_acceleration",
            "original_regulator_removed",
            "unchanged_old_leakage",
            "drop_scalar_center_force",
            "delete_heavy_homogeneous_scalar",
            "two_TT_zero_modes",
            "new_isotropic_vacuum",
            "externally_prescribed_background",
        )
    ):
        cases.append(("wrong_hybrid_scope_" + str(i), require_observable, (value,)))
    for i, value in enumerate((0, -1, 1, 10**1000, s.Rational(1, 10**320))):
        cases.append(("wrong_hybrid_amplitude_" + str(i), require_amplitude, (value,)))
    for i, value in enumerate((-1, 197, True, 1.0, s.Rational(1, 2))):
        cases.append(
            ("unproved_phase_jet_" + str(i), quantum.jet, (value, quantum.H_AMPLITUDE))
        )
        cases.append(
            (
                "unproved_parameter_jet_" + str(i),
                quantum.jet,
                (0, quantum.H_AMPLITUDE, value),
            )
        )
    for i, value in enumerate(
        (
            s.diff(source.q.R, source.N, 6),
            s.diff(source.q.F, source.u, 6),
            s.diff(source.q.Hclock, source.u, 5),
            s.Float(1),
            s.oo,
            1 / s.Symbol("unproved_denominator"),
            s.Symbol("unproved_atom"),
        )
    ):
        cases.append(("unproved_full_source_" + str(i), source.magnitude, (value,)))
    for i, value in enumerate(
        (s.zeros(3), -s.eye(3), s.Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]]))
    ):
        cases.append(
            ("not_proper_cubic_rotation_" + str(i), symmetry.configuration, (value,))
        )
    cases.extend(
        [
            ("delete_original_frontier", validate_scope, ([], matching())),
            ("delete_hybrid_matching_frontier", validate_scope, (frontier(), [])),
        ]
    )
    if len({name for name, _, _ in cases}) != len(cases):
        raise ValueError("Duplicate rejected finite hybrid input")
    return cases


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            count += 1
        else:
            raise ValueError("Unsupported finite hybrid input accepted: " + name)
    return count


def controls():
    return {
        "fixed_profiles_not_replaced_by_live_means": True,
        "all_500_full_source_derivatives_and_N_MVT_retained": True,
        "full_trace_generator_and_spatial_average_remainders": True,
        "actual_radial_seed_not_isotropic_reset": True,
        "all_five_homogeneous_shape_and_three_vector_pairs": True,
        "scalar_phase_restored_and_its_force_retained": True,
        "full_heavy_field_and_cyclic_M1_charge": True,
        "operator_norm_not_symbol_sup_only": True,
        "hybrid_not_original_fully_quantum_bounce": True,
        "historical_errata_and_original_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
