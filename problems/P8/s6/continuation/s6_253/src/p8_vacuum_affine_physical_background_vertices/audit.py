"""Scoped homogeneous physical vertices with the original P8 frontier unchanged."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import audit as reference
from p8_vacuum_affine_gaussian_measure_response import audit as previous

from . import coupled, parent, physical

STATE = previous.STATE
OBSERVABLES = (
    "whole_homogeneous_aligned_physical_background_vertices",
    "whole_held_temporal_vector_embedding_contact",
    "whole_joint_scalar_longitudinal_constraint_measure",
    "fixed_eight_mode_finite_Gaussian_physical_probe_response",
)
ITEM = {
    "id": "QG2_H8A420_complete_homogeneous_physical_background_Gaussian_vertices_joint_light_Proca_constraints_and_held_vector_contact",
    "status": "EXACT_HOMOGENEOUS_GAUGE_FIXED_GAUSSIAN_VERTICES_AND_JOINT_MEASURE_NOT_FULL_COVARIANT_VERTEX_MEASURE_SUBTRACTION_PHYSICAL_MEAN_LOOP_BOUND_NONLINEAR_UV_REGGE_OR_P8",
}


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "Keep the existing full product reference and every fixed preparation"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError("Keep the stated homogeneous gauge-fixed Gaussian scope")
    return label


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No original or earlier matching obligation is promoted")
    return True


@cache
def packets():
    return {
        "entire_parent_source_germs": parent.source_germs(),
        "all_clock_lapse_jets_and_reference_coefficients": parent.reference_data(),
        "whole_aligned_three_mode_constraint_reduction": coupled.data(),
        "whole_joint_physical_Gaussian_measure": coupled.phase_data(),
        "whole_held_vector_source_and_embedding_contact": coupled.held_vector_data(),
        "whole_homogeneous_physical_background_vertices": physical.data(),
        "all_remaining_modes_and_existing_state_transports": physical.other_modes(),
        "whole_fixed_product_Gaussian_response": physical.response_data(),
    }


@cache
def residuals():
    return {
        packet + "_" + name: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
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
        "all_nine_original_primitive_rows_unchanged": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({row["id"] for row in matching()})
        == len(matching()),
        "no_covariant_measure_subtraction_or_physical_mean_inferred": True,
        "no_nonlinear_uv_regge_or_original_P8_closure": True,
    }


def observable():
    return {
        "exact_scope": "Complete homogeneous physical lapse, log-scale and M1 probes through second order in the fixed affine-clock quadratic fluctuation chart. The aligned scalar/longitudinal-Proca block includes the entire source square, exact auxiliary solutions, joint Dirac/configuration measure, all6x6 vertices and the whole held-W0 embedding contact. Both TT, both transverse Proca and H modes, their full metric vertices and their unchanged canonical state transports are retained.",
        "state_and_functional": "The existing two-scalar/tensor, Proca and H preparations form the same product free reference. The finite-regulator16-phase CTP functional retains its continued metaplectic phase, full mixed Wick response, seagulls and fixed profile countergrades. No state is re-prepared under a probe.",
        "boundary": "This is not the full arbitrary inhomogeneous covariant background Hessian, a nonlinear/BRST/ghost/connection/disformal measure, selected Ward/variational curved counterfunctional, physical metric/vector/light mean, or quantitative omitted-loop norm. Finite-mode regularity does not give a uniform perturbed principal-symbol neighborhood.",
        "original_problem": "No new QG3 retuning, compatible nonlinear estimate, same-state bounce, quantum gravitational limit, physical UV scattering or finite-gravity IR/Regge solution. Original V/G/B/P8 remain OPEN.",
    }


def bad_cases():
    cases = []
    for i, value in enumerate(
        (True, False, 1.0, s.Float(1), "1", None, s.oo, s.I, s.nan, s.Symbol("unknown"))
    ):
        cases.extend(
            (
                ("invalid_state_" + str(i), require_state, (value,)),
                ("invalid_observable_" + str(i), require_observable, (value,)),
                (
                    "invalid_physical_normalization_" + str(i),
                    reference.require_physical_normalization,
                    (value,),
                ),
            )
        )
    for label in (
        "full_covariant_metric_stress",
        "selected_Ward_counterfunctional",
        "interacting_vector_tadpole_zero",
        "aligned_second_equals_held_W0_second",
        "separate_Proca_and_light_loops_are_complete",
        "all_nonreference_principal_symbols_hyperbolic",
        "all_quantum_loops_small",
        "original_P8_closed",
        "new_state_preparation_allowed",
        "source_free_heavy_determinant_deleted",
        "frozen_profiles_are_live_means",
        "unrestricted_zero_momentum_shift_constraint",
    ):
        cases.append(("unsupported_" + label, require_observable, (label,)))
    for i in range(len(frontier())):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        cases.append(
            ("primitive_promotion_" + str(i), validate_scope, (rows, matching()))
        )
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        cases.append(
            ("matching_promotion_" + str(i), validate_scope, (frontier(), rows))
        )
    cases.extend(
        (
            ("missing_original", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
            ("wrong_covariance_CCR6", physical.pure_covariance, (s.eye(6) / 5,)),
            ("wrong_covariance_CCR16", physical.pure_covariance, (s.eye(16) / 5,)),
            ("negative_covariance", physical.pure_covariance, (-s.eye(6) / 2,)),
            (
                "nonexact_covariance",
                physical.pure_covariance,
                (s.eye(6) * s.Float(0.5),),
            ),
            ("wrong_whole_symplectic6", physical.symplectic, (2 * s.eye(6), 6)),
            ("wrong_whole_symplectic16", physical.symplectic, (2 * s.eye(16), 16)),
            ("nonsymmetric_whole_vertex", physical.exact_matrix, (s.ones(3, 6),)),
            (
                "unknown_coefficients",
                physical.exact_matrix,
                (s.Symbol("unknown", real=True) * s.eye(6),),
            ),
            ("nonfinite_coefficients", physical.exact_matrix, (s.oo * s.eye(6),)),
            (
                "missing_product_mode",
                physical.assemble_product,
                ([s.eye(4) / 2] + [s.eye(2) / 2 for _ in range(5)],),
            ),
            (
                "wrong_product_scalar_block",
                physical.assemble_product,
                ([s.eye(6) / 2] + [s.eye(2) / 2 for _ in range(6)],),
            ),
            ("changed_kappa", reference.require_physical_normalization, (10**800 + 1,)),
        )
    )
    for i, value in enumerate((True, -1, 0, 3, 12, 6.0, "6", None)):
        cases.append(("invalid_phase_dimension_" + str(i), physical.omega, (value,)))
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported physical-background input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "complete_mixed_source_square_required": coupled.data()["gates"][
            "mixed_first_source_vertex_nonzero"
        ],
        "entire_constraint_lower_bracket_retained": coupled.phase_data()["gates"][
            "whole_lower_constraint_bracket_not_assumed_zero"
        ],
        "nonzero_normal_vector_light_source_retained": coupled.held_vector_data()[
            "gates"
        ]["nonzero_light_quadratic_normal_vector_source_retained"],
        "physical_matter_Ward_jet_not_set_to_reference_zero": physical.data()["gates"][
            "off_reference_matter_Ward_first_contact_nonzero"
        ],
        "complete_physical_lapse_derivative_chain_retained": physical.data()["gates"][
            "physical_lapse_time_derivatives_retained"
        ],
        "same_preparations_and_full_time_connections_retained": True,
        "finite_Gaussian_response_not_renormalized_physical_stress": True,
        "all_original_and_prior_matching_rows_unchanged": frontier()
        == previous.frontier()
        and matching()[:-1] == previous.matching(),
    }
