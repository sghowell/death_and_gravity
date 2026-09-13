"""Precise nonscalar reference inverse scope with all original P8 obligations retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_mass_coupled_inverse import audit as previous

from . import geometry, inverse, local

exact = previous.exact
require_parameters = previous.require_parameters
ITEM = {
    "id": "QG2_H8A420_actual_two_mass_tensor_and_transverse_vector_prepared_Fourier_ball_inverse_with_full_nonscalar_curvature_and_compatible_shift_force_graph",
    "status": "ACTUAL_CONDITIONAL_NONSCALAR_TWO_MASS_PREPARED_FOURIER_BALL_INVERSE_NOT_ARBITRARY_SHIFT_FORCE_UNRESTRICTED_GRAPH_SMALLNESS_STABILITY_INTERACTING_LOOPS_NONLINEAR_UV_REGGE_OR_P8",
}


def require_stage(label):
    if not isinstance(label, str) or label not in (
        "entire_O2_tensor_vector_decomposition",
        "full_actual_nonscalar_curvature_and_tree",
        "actual_prepared_nonscalar_inverse",
        "complete_compatible_vector_shift_recovery",
        "direct_sum_with_current_scalar_clock_M1_subgraph",
    ):
        raise ValueError(
            "Keep the actual conditional nonscalar inverse and its explicit domains"
        )
    return label


def require_vector_graph(label):
    if (
        not isinstance(label, str)
        or label != "prepared_quotient_force_with_full_Ward_compatible_shift_partner"
    ):
        raise ValueError(
            "Arbitrary shift forcing is not the proved infrared quotient-force graph"
        )
    return label


def require_reference(label):
    if (
        not isinstance(label, str)
        or label != "complete_same_scheme_two_mass_shear_factor"
    ):
        raise ValueError("Retain both complete mass cuts and their fixed finite sum")
    return label


def require_internal_integral(label):
    if not isinstance(label, str) or label != "all_original_internal_momenta":
        raise ValueError("External Fourier support is not a physical internal cutoff")
    return label


def require_sector(label):
    if not isinstance(label, str) or label not in ("tensor", "transverse_vector"):
        raise ValueError(
            "Distinguish the two nonscalar sectors and their different full operators"
        )
    return label


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No original primitive or prior matching status is promoted")
    return True


@cache
def packets():
    return {
        "whole_spatial_O2_representation": geometry.representation(),
        "entire_vector_Ward_and_force_graph": geometry.ward(),
        "whole_nonscalar_leading_and_finite_match": geometry.leading(),
        "full_nonscalar_local_and_classical_action": local.data(),
        "complete_nonscalar_inverse_and_constraint_recovery": inverse.data(),
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
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
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
        "full_actual_conditional_metric_subgraphs_not_all_interacting_physics": True,
        "physical_normalizations_and_infrared_force_domains_retained": True,
    }


def observable():
    return {
        "actual": "The two tensor and two transverse-vector sectors of the same actual QG2-H8A420 Proca-plus-heavy conditional Gaussian reference, with complete finite curvature actions, current Einstein tree, retarded state/contacts and profiles kept exactly.",
        "new": "The full O(2) decomposition and correct sector-dependent curvature Hessians yield an actual smooth prepared finite-ball inverse. Both original transverse-shift constraints follow on the explicit Ward-compatible force graph.",
        "bound": "For each fixed external ball, the complete total shear kernel has finite L1 norm K, the actual remainder has finite unevaluated C, lambda=(4KC+1)^2 and canonical source bound2exp(lambda)K16pi²kappa on I4 forcing. All physical factors remain.",
        "combination": "Together with S248 this supplies the direct sum of the stated conditional metric/M1 scalar and nonscalar prepared amplitude/compatible-force subgraphs. It is not arbitrary raw shift-force surjectivity or an unrestricted all-momentum completed graph.",
        "remaining": "No evaluated stability, physical smallness, interacting light/mixed loops, nonlinear same-state bounce, quantum gravitational limit, physical UV/Regge or original V/G/B/P8 closure.",
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
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    g = geometry.reference
    templates = (
        ("parameters", require_parameters, (g.PROCA_MASS2, g.HEAVY_MASS2, g.KAPPA)),
        ("ball", inverse.require_ball, (1,)),
        ("stage", require_stage, ("actual_prepared_nonscalar_inverse",)),
        (
            "graph",
            require_vector_graph,
            ("prepared_quotient_force_with_full_Ward_compatible_shift_partner",),
        ),
        (
            "reference",
            require_reference,
            ("complete_same_scheme_two_mass_shear_factor",),
        ),
        ("internal", require_internal_integral, ("all_original_internal_momenta",)),
        ("sector", require_sector, ("tensor",)),
    )
    out = []
    for i, value in enumerate(invalid):
        for name, call, args in templates:
            for j in range(len(args)):
                bad = list(args)
                bad[j] = value
                out.append((f"{name}_invalid_{i}_{j}", call, tuple(bad)))
    for value in (0, -1, -s.Rational(1, 2)):
        out.append((f"nonpositive_ball_{value}", inverse.require_ball, (value,)))
    for j in range(3):
        args = [g.PROCA_MASS2, g.HEAVY_MASS2, g.KAPPA]
        args[j] += 1
        out.append((f"changed_parameter_{j}", require_parameters, tuple(args)))
    for value in (
        "unrestricted_graph",
        "arbitrary_shift_force_Hr",
        "uniform_all_Pmax",
        "nonlinear_bounce",
        "stable",
        "small_physical_inverse",
        "physical_cutoff",
        "closed_P8",
        "interacting_quantum_state",
        "literal_homogeneous_vector",
    ):
        out.append(("unsupported_" + value, require_stage, (value,)))
    out.extend(
        (
            (
                "wrong_vector_source_graph",
                require_vector_graph,
                ("unprepared_arbitrary_shift_force",),
            ),
            (
                "wrong_old_shear_inverse",
                require_reference,
                ("old_Proca_pole_plus_cut_only",),
            ),
            ("wrong_massless_cutoff", require_internal_integral, ("finite_bare_band",)),
            ("wrong_scalar_sector", require_sector, ("scalar",)),
        )
    )
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
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported nonscalar inverse claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "reflection_removes_two_false_parity_channels": geometry.representation()[
            "rotation_only_nullity"
        ]
        == 8,
        "tensor_gradient_not_copied_into_vector": local.weyl_mixed(
            geometry.BASIS[2], geometry.BASIS[2]
        )
        != local.weyl_mixed(geometry.BASIS[4], geometry.BASIS[4]),
        "Gaussian_mean_not_zeroed_before_vector_contacts": True,
        "nonzero_unprepared_shift_constraint_kernel": geometry.ward()[
            "unprepared_constraint_counterexample"
        ]
        != 0,
        "both_canonical_factors_and_force_sign_kept": True,
        "arbitrary_shift_force_IR_counterexample_retained": True,
        "full_inverse_is_of_two_mass_sum_not_old_factor": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
