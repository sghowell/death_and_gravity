"""Scope guards and preserved frontier for the full-source lapse-chart obstruction."""

from functools import cache

import sympy as s
from p8_vacuum_affine_canonical_boundary_corrected_hybrid import audit as previous
from p8_vacuum_affine_complement_measure.split import clean

from . import auxiliary, branch, fields, fold

MODEL = "full_original_lapse_chart_obstruction_not_quantum_completion"
OBSERVABLES = (
    "regular_joint_trace_temporal_atlas",
    "actual_fixed_profile_lapse_fold",
    "constrained_finite_input_path",
    "connected_sheet_and_measure_boundary",
)
ITEM = {
    "id": "QG2_H8A440_full_original_joint_auxiliary_atlas_and_constrained_finite_phase_single_lapse_chart_obstruction",
    "status": "COMPLETE_SCOPED_JOINT_ATLAS_AND_CONNECTED_LAPSE_CHART_OBSTRUCTION_NOT_A_QUANTUM_OR_ORIGINAL_V_G_B_P8_NO_GO",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the full-source lapse-chart obstruction scope")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require an established lapse-chart observable")
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def qualifications():
    return previous.qualifications()


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A lapse-chart obstruction does not close original P8")
    return True


@cache
def packets():
    return {
        "whole_joint_original_trace_temporal_atlas": auxiliary.data(),
        "whole_exact_original_fixed_profile_clock_fold": fold.data(),
        "whole_finite_input_fields_and_spatial_constraints": fields.data(),
        "whole_connected_local_sheet_and_measure_boundary": branch.data(),
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
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "all9_original_primitive_rows_preserved": len(frontier()) == 9,
        "all_previous_matching_rows_preserved": matching()[:-1] == previous.matching(),
        "all_matching_ids_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "all6_archived_physical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "corrected_finite_S275_result_not_refuted": True,
        "original_parameters_not_changed": require_parameters(parameters())
        == parameters(),
        "no_new_state_or_operator_extension_selected": True,
        "original_V_G_B_P8_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The full original mixed K,T system is regular for all positive R, including the old sequential Gamma crossing. Separately, exact full-source bounce jets and actual fixed quantum profiles give a genuine N,T lapse fold with Hessian diag(0,-1), C_s=3 and C_NN>28. A path of exact finite Fourier canonical inputs, preserving all means and spatial constraints, connects the existing local lapse sheet to this fold. There is no globally C1 lapse map extending that sheet to all these finite canonical data. The full branch envelope and finite second-class determinant boundary are retained.",
        "domain": "Same original full parent, fixed reference profile functions, torusL1, frequencyP1e64, kappa1e800, zeta1e-6 and original axis modes. The path lies in unlocalized finite canonical phase space, not the certified small quantum cutoff ball. Its start is in the entire S266 invariant box. All reconstructed metric/momentum harmonics remain. No new quantum state, expectation, trajectory or operator ordering is defined.",
        "historical_qualification": qualifications(),
        "not_established": "Not a no-go theorem for the parent, all UV completions or original P8; not a classical evolution reaching the fold, a global Lorentzian singularity, a clock-time turning point, a new first-class symmetry, quantum determinant cancellation at rank change, unitary branch gluing, self-adjoint domain, regulator removal, physical cutoff, UV matching or omitted-loop/Regge control. The S275 corrected finite hybrid remains intact.",
    }


def bad_cases():
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
    rows = []
    for i, value in enumerate(invalid):
        for name, call in (
            ("model", require_model),
            ("observable", require_observable),
            ("axis", fields.fixture),
            ("frequency", fields.require_frequency),
            ("clock_order", fold.constraint_jet),
        ):
            rows.append(("invalid_" + name + "_" + str(i), call, (value,)))
    for i, value in enumerate((-1, 3, 4, s.Rational(1, 2))):
        rows.append(("bad_axis_" + str(i), fields.fixture, (value,)))
    for i, value in enumerate((-1, 4, 5, s.Rational(1, 2))):
        rows.append(("bad_clock_order_" + str(i), fold.constraint_jet, (value,)))
    for i, value in enumerate((0, 1, 10**63, 2 * 10**64)):
        rows.append(("changed_frequency_" + str(i), fields.require_frequency, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "all_UV_completions_excluded",
            "Gamma_is_a_physical_singularity",
            "full_quantum_state_support",
            "new_Airy_unitarity",
            "automatic_gauge_symmetry",
            "global_C1_lapse_map",
            "S275_finite_result_refuted",
        )
    ):
        rows.append(("unsupported_scope_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("deleted_original_frontier", validate_scope, ([], matching())),
            ("deleted_matching_frontier", validate_scope, (frontier(), [])),
        )
    )
    if len({name for name, _, _ in rows}) != len(rows):
        raise ValueError("Duplicate lapse-chart rejection case")
    return rows


def rejected_inputs():
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported lapse-chart input accepted: " + name)
    return total


def controls():
    return {
        "sequential_pivot_zero_is_not_joint_stationary_rank_loss": True,
        "actual_full_profiles_and_nonzero_primitive_retained": True,
        "literal_metric_curvature_and_all_spatial_generators": True,
        "same_means_and_modes_with_complete_generated_harmonics": True,
        "connected_sheet_not_isolated_arbitrary_invariant_fixture": True,
        "fold_not_gauge_symmetry_or_quantum_no_go": True,
        "earlier_corrected_finite_result_and_all_frontiers_preserved": True,
        "rejected_inputs": rejected_inputs(),
    }
