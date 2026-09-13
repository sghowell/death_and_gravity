"""Read-only physical free-state boundary and unchanged original P8 frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_source_filtration import audit as previous

from . import gaussian, phase, uv

ITEM = {
    "id": "QG2_H8A420_complete_coupled_scalar_and_tensor_free_reference_CCR_state_with_fixed_matrix_sampling_and_all_order_two_cone_spectrum",
    "status": "EXACT_PHYSICAL_QUADRATIC_REFERENCE_GAUSSIAN_STATE_AND_WRITTEN_REDUCED_TWO_CONE_SHORT_DISTANCE_PROOF_NOT_COVARIANT_COUNTERFUNCTIONAL_INTERACTING_STATE_LOOP_BOUND_NONLINEAR_UV_REGGE_OR_P8",
}
STATE = "QG2_H8A420_FREE_REF_S251_FIXED_SAMPLING_AND_CAUCHY_DATA"
OBSERVABLES = (
    "complete_coupled_scalar_CCR_reference",
    "two_physical_tensor_CCR_references",
    "reduced_oriented_two_cone_short_distance_condition",
)


def require_state(label):
    if not isinstance(label, str) or label != STATE:
        raise ValueError(
            "Keep the whole fixed sampled reference, not a live state reset"
        )
    return label


def require_observable(label):
    if not isinstance(label, str) or label not in OBSERVABLES:
        raise ValueError(
            "Keep the explicit physical free-reference and reduced-spectrum scope"
        )
    return label


def require_physical_normalization(value):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, s.Integer, s.Rational))
        or value != 10**800
    ):
        raise ValueError("Keep the exact current kappa at both physical endpoints")
    return s.Integer(value)


def require_sampling(left, right):
    if any(
        isinstance(x, bool) or not isinstance(x, (int, s.Integer, s.Rational))
        for x in (left, right)
    ):
        raise ValueError("Keep the exact full preparation interval")
    if (left, right) != (-s.Rational(15, 32), -s.Rational(13, 32)):
        raise ValueError("Changing the bump support changes this specified state")
    return left, right


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No original or previous physical obligation is promoted")
    return True


@cache
def packets():
    return {
        "complete_current_canonical_phase_and_preparation": phase.data(),
        "whole_matrix_Gaussian_covariance_and_CCR": gaussian.data(),
        "full_double_cone_symbol_and_Riccati_input": uv.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for name, data in packets().items()
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
            name + "_" + key: bool(value)
            for name, data in packets().items()
            for key, value in data["gates"].items()
        },
        "all_nine_original_primitive_rows_retained": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({r["id"] for r in matching()})
        == len(matching()),
        "no_covariant_curved_counterfunctional_or_interacting_state_inferred": True,
        "no_omitted_loop_size_or_nonlinear_closure_from_free_state": True,
    }


def observable():
    return {
        "exact": "The complete current physical two-scalar local Hamiltonian and two classical tensor modes admit the fixed sampled positive pure Gaussian reference, retaining full mixing, both normalization factors, all time boundaries and original coherent M1 mean.",
        "continuum": "The written all-order homogeneous matrix proof supplies the oriented two-cone short-distance condition across the complete chart cover, with full noncommuting lower coefficients and smooth anomalous-block sampling. Finite diagnostics are not its formalization.",
        "not_a_live_preparation": "Reference Cauchy covariances are selected once and held fixed; no re-minimization, operational measurement protocol, instantaneous vacuum or future response boundary condition is supplied.",
        "not_covariant_subtraction": "This reduced physical free state is not a covariant gauge-fixed propagator, full constraint/ghost determinant, curved counterfunctional or renormalized physical stress calculation.",
        "remaining": "Complete nonlocal/interacting state, all omitted-loop and compatible nonlinear estimates, same-state bounce, physical quantum gravitational limit, UV scattering, finite-gravity IR/Regge and original P8.",
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
    templates = (
        ("state", require_state, (STATE,)),
        ("observable", require_observable, (OBSERVABLES[0],)),
        ("normalization", require_physical_normalization, (s.Integer(10**800),)),
        ("sampling", require_sampling, (-s.Rational(15, 32), -s.Rational(13, 32))),
    )
    out = []
    for i, value in enumerate(invalid):
        for name, call, args in templates:
            for j in range(len(args)):
                changed = list(args)
                changed[j] = value
                out.append((f"{name}_invalid_{i}_{j}", call, tuple(changed)))
    for label in (
        "complete_interacting_state",
        "full_curved_covariant_counterfunctional",
        "all_quantum_loops_small",
        "nonlinear_bounce_solved",
        "original_P8_closed",
        "single_scalar_Hadamard_theorem_applies_automatically",
        "instantaneous_vacuum",
        "full_nonlocal_quantum_response_state",
        "all_times_one_Fock_implementation",
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
            ("changed_kappa", require_physical_normalization, (10**800 + 1,)),
            (
                "changed_sampling",
                require_sampling,
                (-s.Rational(1, 2), -s.Rational(13, 32)),
            ),
            ("live_state_reset", require_state, ("reminimize_on_varied_history",)),
            ("wrong_matrix_shape", gaussian.ground_covariance, (s.eye(2),)),
            (
                "nonpositive_preparation",
                gaussian.ground_covariance,
                (s.diag(1, 1, 1, -1),),
            ),
            (
                "nonsymmetric_preparation",
                gaussian.ground_covariance,
                (s.eye(4) + s.eye(4).elementary_row_op("n->n+km", row=0, k=1, row2=1),),
            ),
            (
                "uncontrolled_float_preparation",
                gaussian.ground_covariance,
                (s.eye(4) * s.Float(1),),
            ),
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
        raise ValueError("Unsupported free-reference claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "missing_volume_connection_breaks_CCR": phase.data()["gates"][
            "omitting_density_time_connection_breaks_CCR"
        ],
        "real_positivity_not_complex_uncertainty": gaussian.data()["gates"][
            "strictly_positive_real_covariance_is_not_enough"
        ],
        "anomalous_blocks_cannot_be_deleted": gaussian.data()["gates"][
            "deleting_anomalous_blocks_changes_generic_state"
        ],
        "full_two_cones_not_single_normal_hyperbolicity": uv.data()["gates"][
            "two_cones_not_scalar_normally_hyperbolic_symbol"
        ],
        "complete_lower_kinetic_terms_retained": uv.data()["gates"][
            "finite_q_kinetic_corrections_not_deleted"
        ],
        "initial_state_not_reselected_on_varied_history": True,
        "curved_counterfunctional_and_loop_sizes_still_required": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
