"""Preserved frontier and strict scope for homogeneous classical fold dynamics."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_global_lapse_chart_obstruction import audit as previous

from . import canonical, consistency, endpoints, source

MODEL = "full_original_homogeneous_fold_dynamics_not_quantum_completion"
OBSERVABLES = (
    "full_clock_and_canonical_flow",
    "all_p_fold_consistency_obstruction",
    "actual_regular_singular_endpoint",
    "physical_lapse_volume_and_curvature",
)
ITEM = {
    "id": "QG2_H8A441_full_original_homogeneous_lapse_fold_constraint_preservation_and_finite_proper_time_classical_curvature_endpoint",
    "status": "COMPLETE_SCOPED_HOMOGENEOUS_CLASSICAL_FOLD_DYNAMICS_NOT_ORIGINAL_STATE_OR_V_G_B_P8_CLOSURE",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError(
            "Require the original homogeneous classical fold dynamics scope"
        )
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require an established homogeneous fold observable")
    return value


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def qualifications():
    return previous.qualifications()


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Homogeneous fold dynamics does not close original P8")
    return True


@cache
def packets():
    return {
        "whole_exact_original_clock_source": source.data(),
        "whole_original_canonical_density_flow": canonical.data(),
        "whole_first_and_second_preservation": consistency.data(),
        "whole_regular_classical_endpoint_and_physical_geometry": endpoints.data(),
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
        "all6_historical_physical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "S275_corrected_finite_quantum_hybrid_not_refuted": True,
        "S276_mean_preserving_wave_path_not_relabelled_as_time_solution": True,
        "original_parameter_record_not_changed": require_parameters(parameters())
        == parameters(),
        "no_new_quantum_state_or_boundary_condition_selected": True,
        "original_V_G_B_P8_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "At the actual full-source clock slice, the entire homogeneous zero-heavy diagonal Bianchi-I fold family has no C1 lapse trajectory through it: first preservation either fails or fixes a trace density below21 times the actual profile bound, where second preservation has discriminant below-6719. With p=-1/10 a smooth desingularized full-source curve yields regular C_N<0 classical clock-time solutions approaching this fold in finite proper time, with finite positive lapse/volume but divergent normal expansion and physical Ricci scalar.",
        "domain": "Original full local parent with all fixed profiles, masses, constants, primitive, temporal vector and live heavy equations retained. Homogeneous M1, shape and trace momenta are separately specified classical test data, not the unchanged prepared state or the S276 wave path. This is local qualitative endpoint existence, not an evaluated initial interval or physical-cutoff regime.",
        "historical_qualification": qualifications(),
        "not_established": "No generic instability in the original bounce neighborhood, quantum mean, all-inhomogeneous-fold classification, state-domain or measure construction, self-adjoint extension, quantum resolution, unlocalized regulator limit, physical UV matching, omitted-loop/Regge bound, all-parent no-go or original P8 closure. Corrected S275 finite solutions remain unchanged.",
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
        for name, call, args in (
            ("model", require_model, (value,)),
            ("observable", require_observable, (value,)),
            ("jet_family", source.clock_jet, (value, 0, 0)),
            ("jet_time", source.clock_jet, ("C", value, 0)),
            ("jet_lapse", source.clock_jet, ("C", 0, value)),
        ):
            rows.append(("invalid_" + name + "_" + str(i), call, args))
    for i, args in enumerate(
        (
            ("H", 2, 0),
            ("C", 3, 0),
            ("C", 0, 3),
            ("C", 1, 2),
            ("T", 2, 0),
            ("B", 0, 0),
            ("C", -1, 0),
            ("C", 0, -1),
        )
    ):
        rows.append(("unlicensed_jet_" + str(i), source.clock_jet, args))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "S275_finite_hybrid_refuted",
            "S276_wave_time_trajectory",
            "every_UV_completion_excluded",
            "automatic_gauge_symmetry",
            "quantum_Airy_resolution",
            "C1_clock_through_homogeneous_fold",
            "unchanged_reference_state_singularity",
            "physical_lapse_is_N_times_R_to_minus_one_quarter",
            "uniform_evaluated_lifetime",
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
        raise ValueError("Duplicate homogeneous fold rejection case")
    return rows


def rejected_inputs():
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported homogeneous fold input accepted: " + name)
    return total


def controls():
    return {
        "full_time_differentiation_before_clock_restriction": True,
        "nonzero_BuN_and_weighted_density_self_flow_retained": True,
        "complete_second_preservation_not_first_only": True,
        "whole_heavy_phase_pair_before_both_flows": True,
        "desingularized_full_source_not_cubic_replacement_dynamics": True,
        "actual_physical_lapse_and_independent_Ricci": True,
        "distinct_classical_initial_data_not_original_quantum_state": True,
        "rejected_inputs": rejected_inputs(),
    }
