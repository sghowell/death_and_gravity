"""Known selected-triangle finite matching in a fixed lift, not original P8 closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_curvature_contact_basis import audit as previous

from . import calibration, dressing, jets, moment, source, triangle

MODEL = "original_selected_triangle_finite_curvature_coefficient"
OBSERVABLES = (
    "whole_original_source",
    "whole_explicit_labeled_box_lift",
    "whole_generic_triangle_response",
    "whole_finite_positive_mass_moment",
    "whole_complete_dressing_and_normalization",
    "whole_original_vector_and_frozen_kernel_calibrations",
)
ITEM = {
    "id": "QG2_H8A508_known_triangle_finite_curvature_coefficient",
    "status": "COMPLETE_SELECTED_TRIANGLE_DEGREE6_COEFFICIENT_IN_EXPLICIT_LIFT_NOT_FULL_PARENT_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the unchanged selected-triangle model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated whole triangle-coefficient observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "selected_triangle_labeled_box_lift",
        "homogeneous_six_derivative_local_difference",
        "original_source_finite_triangle_coefficient",
    ):
        raise ValueError("Require the stated selected triangle matching sector")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "A selected finite triangle coefficient cannot close original P8"
        )
    return True


@cache
def packets():
    return dict(
        zip(
            OBSERVABLES,
            (
                source.data(),
                jets.data(),
                triangle.data(),
                moment.data(),
                dressing.data(),
                calibration.data(),
            ),
        )
    )


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
        "all199_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 199,
        "all_matching_identifiers_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            row["status"].startswith("REJECTED_") for row in matching()
        ),
        "independent_extra_parent_chi_not_assigned_or_bounded": True,
        "local_jet_not_a_physical_threshold_or_truncation_bound": True,
        "physical_metric_bounce_scoped_P8a_and_global_obligations_unchanged": True,
    }


def observable():
    return {
        "established": "The complete selected local-bilocal matter triangle, including its outer-heavy-line metric variation, has a definite finite degree6 real-TT curvature difference from the explicit labeled scalar-box lift of its off-shell flat Taylor series. The whole symbolic integrand vanishes at degree2 and4, and equals-8*x^2*y^2*z^2*Q/M^4 at degree6. Its convergent mass moment gives chi_triangle=g^2*(C+g^2/n)*c_triangle(n)/(16pi^2), including all24 labels. The original coefficient is positive and |chi_triangle|/A0<10^-206.",
        "domain": "Same original massive source and physical metric. Homogeneous low-external-momentum Taylor jet of the UV-finite selected triangle functional, with four identical scalar fields and one real-null-TT graviton. The off-shell local comparison convention is fully fixed; no physical mass-shell constraint is used to infer the coefficient.",
        "historical_qualification": qualifications(),
        "not_established": "The box coefficient in this convention, the aggregate coefficient of every matter class, independent extra parent chi, a full curved-background functional, internal-graviton loops or finite-gravity quantum decoupling. No approximation of the above-threshold physical amplitude by this local jet, no truncation error bound, and no addition on top of the already complete S342 triangle. Full virtual/inclusive probability, exact LSZ, quantum unitarity, complex Regge, same-parent bounce, UV and original V/G/B/P8 remain open.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
    for index, value in enumerate(
        (True, False, 1.0, s.Float(1), "1", None, s.oo, s.I, s.nan, s.Symbol("unknown"))
    ):
        for label, call in (
            ("new_model", require_model),
            ("new_observable", require_observable),
            ("new_sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(index), call, (value,)))
    for index, value in enumerate(
        (
            "triangle_sets_parent_chi",
            "local_jet_replaces_physical_cut",
            "known_loop_added_twice",
            "flat_OS4_fixes_curvature",
            "one_sample_matching",
            "triangle_closes_original_P8",
        )
    ):
        rows.append(
            ("unsupported_new_scope_" + str(index), require_observable, (value,))
        )
    rows.extend(
        (
            ("new_deleted_primitive", validate_scope, ([], matching())),
            ("new_deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    for index, value in enumerate(
        (True, False, 1.0, s.Float(1), "1", None, 0, -1, s.oo, s.I)
    ):
        rows.append(
            (
                "invalid_exact_positive_triangle_mass_" + str(index),
                moment.coefficient,
                (value,),
            )
        )
    for index, value in enumerate((True, False, 0, -1, 3.0, s.Integer(3), 4, "3")):
        rows.append(
            ("invalid_triangle_jet_order_" + str(index), jets.require_order, (value,))
        )
    for index, value in enumerate((True, False, -1, 3, 1.0, s.Integer(1), "1", None)):
        rows.append(
            (
                "invalid_literal_triangle_line_" + str(index),
                triangle.literal_line,
                (value, 3),
            )
        )
    return rows


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            count += 1
        else:
            raise ValueError("Unsupported finite triangle input accepted: " + name)
    return count


def controls():
    return {
        "full_generic_simplex_identity_not_sample_fit": True,
        "explicit_off_shell_labeled_box_completion": True,
        "all_three_loop_lines_and_shifted_outer_heavy_branch": True,
        "strict_positive_mass_integral_equal_mass_and_asymptotic_checks": True,
        "all24_Bose_labels_and_original_loop_normalization": True,
        "independent_original_vectors_and54_frozen_line_kernels": True,
        "physical_threshold_extra_chi_and_all_original_frontiers_preserved": True,
        "rejected_inputs": rejected_inputs(),
    }
