"""Scoped full Gaussian functional and unchanged original P8 frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import audit as previous

from . import functional, measure, response

STATE = previous.STATE
OBSERVABLES = (
    "complete_quadratic_auxiliary_Gaussian_measure",
    "fixed_state_whole_coefficient_CTP_kernels",
    "continuum_reference_Wick_noise",
    "time_retarded_distribution_extension_family",
)
ITEM = {
    "id": "QG2_H8A420_complete_reference_Gaussian_measure_CTP_vertices_and_continuum_Wick_noise_with_retarded_extension_family",
    "status": "EXACT_QUADRATIC_MEASURE_AND_FINITE_CTP_WITH_WRITTEN_CONTINUUM_WICK_NOISE_AND_EXTENSION_FAMILY_NOT_COVARIANT_COUNTERFUNCTIONAL_PHYSICAL_LOOP_MEAN_BOUND_NONLINEAR_UV_REGGE_OR_P8",
}


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError("Keep the same fixed S251 Cauchy state")
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError(
            "Keep the explicit quadratic coefficient and Wick-extension scope"
        )
    return label


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No original or prior matching obligation is promoted")
    return True


@cache
def packets():
    return {
        "whole_physical_auxiliary_Gaussian_measure": measure.data(),
        "whole_fixed_state_quadratic_response_and_noise": response.data(),
        "complete_CTP_coefficient_vertices_and_continuum_Wick_boundary": functional.data(),
    }


@cache
def residuals():
    return {
        packet + "_" + key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for packet, data in packets().items()
        for key, value in data["checks"].items()
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
            packet + "_" + key: bool(value)
            for packet, data in packets().items()
            for key, value in data["gates"].items()
        },
        "all_nine_original_primitive_rows_preserved": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({r["id"] for r in matching()})
        == len(matching()),
        "no_selected_covariant_counterfunctional_or_physical_mean_inferred": True,
        "no_loop_size_nonlinear_or_original_P8_closure": True,
    }


def observable():
    return {
        "exact_finite": "Entire current quadratic lapse/shift Dirac and configuration measure, full fixed-state Gaussian CTP trace with continued metaplectic phase, complete independent coefficient first/second vertices, mean, time-retarded susceptibility, seagulls and connected noise.",
        "continuum": "Written two-cone contraction-graph and smoothing proof defines the specified reference Wick bilinears and positive connected noise without a momentum cutoff. Full chart symbol orders supply a conservative finite diagonal scaling bound and a family of time-retarded extensions.",
        "physical_boundary": "Independent coefficient probes and reference normal ordering are not a covariant physical metric/light stress. Full parent vertex maps, nonlinear off-shell contacts, covariant/gauge measure, Ward/variational-compatible curved subtraction and finite conditions remain separate.",
        "not_selected": "The local coefficients of the retarded extension family are not fixed. Its existence is not the construction of a physical counterfunctional, a numerical loop bound, spacelike microcausality or an interacting state.",
        "remaining": "Existing H/Proca/M1 data remain fixed. No QG3 retuning, compatible nonlinear estimate, same-state nonlinear bounce, quantum gravitational limit, physical UV scattering, finite-gravity IR/Regge or original V/G/B/P8 closure.",
    }


def bad_cases():
    out = []
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
        out.extend(
            (
                (f"invalid_state_{i}", require_state, (value,)),
                (f"invalid_observable_{i}", require_observable, (value,)),
                (
                    f"invalid_normalization_{i}",
                    previous.require_physical_normalization,
                    (value,),
                ),
            )
        )
    for label in (
        "full_covariant_BRST_measure",
        "selected_curved_counterfunctional",
        "renormalized_physical_stress_zero",
        "all_quantum_loops_small",
        "complete_interacting_state",
        "same_state_nonlinear_bounce",
        "original_P8_closed",
        "normal_ordering_is_QG3_retuning",
        "coefficient_derivative_is_metric_stress",
        "principal_endpoint_square_root_is_full_CTP_phase",
        "all_continuum_functional_determinants_multiply",
        "time_retarded_means_spacelike_microcausal",
    ):
        out.append(("unsupported_" + label, require_observable, (label,)))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_promotion_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_promotion_{i}", validate_scope, (frontier(), rows)))
    out.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
            ("live_state_reset", require_state, ("reminimize_on_varied_history",)),
            ("wrong_pure_CCR", response.pure_covariance, (s.eye(4) / 5,)),
            ("negative_covariance", response.pure_covariance, (-s.eye(4) / 2,)),
            ("float_covariance", response.pure_covariance, (s.eye(4) * s.Float(0.5),)),
            ("odd_phase", response.omega, (3,)),
            ("boolean_phase", response.omega, (True,)),
            ("noncanonical_propagator", response.symplectic, (2 * s.eye(4), 4)),
            ("nonsymmetric_vertex", response.exact_matrix, (s.ones(3, 2),)),
            (
                "singular_auxiliary",
                measure.finite_dirac_block,
                (s.zeros(2), s.zeros(2)),
            ),
            (
                "nonskew_constraint_block",
                measure.finite_dirac_block,
                (s.eye(2), s.eye(2)),
            ),
            (
                "nonexact_constraint_block",
                measure.finite_dirac_block,
                (s.eye(2) * s.Float(1), s.zeros(2)),
            ),
            (
                "formal_unknown_constraint_block",
                measure.finite_dirac_block,
                (s.eye(2) * s.Symbol("r", real=True), s.zeros(2)),
            ),
            ("changed_kappa", previous.require_physical_normalization, (10**800 + 1,)),
        )
    )
    for i, value in enumerate((True, -1, 1.5, "2", None, s.oo)):
        out.append(
            (
                f"unsupported_bilinear_order_{i}",
                functional.bilinear_scaling_bound,
                (value, 4),
            )
        )
    for i, value in enumerate((True, 0, 7, 2.0, None)):
        out.append(
            (
                f"unsupported_formal_depth_{i}",
                functional.logdet_coefficients,
                ([s.eye(2)], value),
            )
        )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported quantum-functional claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_auxiliary_lower_bracket_retained": measure.data()["gates"][
            "constraint_lower_commutator_not_deleted"
        ],
        "missing_configuration_Jc_factor_changes_result": measure.data()["gates"][
            "missing_Jc_measure_changes_configuration_functional"
        ],
        "Wick_transpose_not_conjugate_transpose": response.data()["gates"][
            "transpose_not_adjoint_negative_control"
        ],
        "CTP_phase_not_only_endpoint_principal_root": functional.data()["gates"][
            "endpoint_principal_root_would_lose_metaplectic_sign"
        ],
        "complete_second_vertices_and_off_shell_contacts_required": True,
        "fixed_state_and_formal_countergrades_preserved": True,
        "extension_family_not_selected_covariant_subtraction": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
