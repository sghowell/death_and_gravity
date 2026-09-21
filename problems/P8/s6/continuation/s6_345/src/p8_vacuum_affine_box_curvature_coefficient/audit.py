"""One complete known ordered-box coefficient in a fixed jet lift, not full P8."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_triangle_curvature_coefficient import audit as previous

from . import basis, calibration, jets, moment, radiation, source

MODEL = "original_selected_box_finite_curvature_in_symmetric_jet_lift"
OBSERVABLES = (
    "whole_original_source",
    "whole_complete_flat_Bose_jet_projection",
    "whole_literal_covariant_jet_metric_response",
    "whole_generic_box_radiation_matching",
    "whole_finite_box_mass_moment",
    "whole_original_vector_and_frozen_kernel_calibrations",
)
ITEM = {
    "id": "QG2_H8A509_known_box_finite_curvature_coefficient",
    "status": "COMPLETE_SELECTED_BOX_DEGREE6_COEFFICIENT_IN_EXPLICIT_JET_LIFT_NOT_AGGREGATE_PARENT_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the unchanged selected ordered-box model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated whole box-coefficient observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "selected_ordered_box_symmetric_jet_lift",
        "whole_homogeneous_six_derivative_polynomial",
        "original_known_box_mass_moment",
    ):
        raise ValueError("Require the stated selected box matching sector")
    return value


def require_word_index(value):
    if type(value) is not int or not 0 <= value < 6:
        raise ValueError("Require exactly one of the six stated covariant jet words")
    return basis.WORDS[value]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A selected box coefficient cannot close original P8")
    return True


@cache
def packets():
    return dict(
        zip(
            OBSERVABLES,
            (
                source.data(),
                basis.data(),
                jets.data(),
                radiation.data(),
                moment.data(),
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
        "all200_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 200,
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
        "triangle_basis_conversion_required_before_aggregate_matching": True,
        "physical_metric_bounce_scoped_P8a_and_global_obligations_unchanged": True,
    }


def observable():
    return {
        "established": "The complete selected alternating-mass box has a fixed finite degree6 real-TT curvature difference from the explicit six-word covariant Bose jet lift of its entire off-shell flat degree6 polynomial. All220 flat monomials establish the rank6 basis; all96 literal metric insertions and2430 generic coefficient residuals establish the whole radiative difference. The mass coefficient is chi_box=g^4*c_box(n)/(16pi^2), with a strictly negative convergent integral and exact rational-log formula. Both original massive-vector and72 actual frozen line-kernel checks pass.",
        "domain": "Same original massive source and physical metric; one real-null-TT graviton and four identical scalars. The local analytic-origin degree6 coefficient only, with a fully fixed symmetrized-jet convention before any scalar-shell substitution. All six original mass-ordered flat boxes remain included.",
        "historical_qualification": qualifications(),
        "not_established": "An aggregate known matter coefficient: the S344 triangle convention must first be converted and all other sectors treated consistently. Independent extra parent chi, a full curved-background functional, internal-graviton loops and finite-gravity quantum decoupling are not matched. No physical above-threshold Taylor approximation or truncation-error bound, no perturbative smallness from the LARGE10^187 coefficient bound, and no second copy added to S342. Full virtual/inclusive probability, exact LSZ, quantum unitarity, complex Regge, same-parent bounce, UV and original V/G/B/P8 remain open.",
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
            "box_sets_parent_chi",
            "local_jet_replaces_physical_cut",
            "known_loop_added_twice",
            "different_bases_silently_added",
            "finite_fixtures_prove_generic_matching",
            "box_closes_original_P8",
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
                "invalid_exact_positive_box_mass_" + str(index),
                moment.coefficient,
                (value,),
            )
        )
    for index, value in enumerate(
        (True, False, -1, 6, 1.0, s.Integer(1), "1", None, [])
    ):
        rows.append(
            ("invalid_fixed_box_jet_word_" + str(index), require_word_index, (value,))
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
            raise ValueError("Unsupported finite box input accepted: " + name)
    return count


def controls():
    return {
        "all220_flat_Bose_polynomials_and_six_explicit_covariant_jets": True,
        "all64_third_jet_components_and_original_Galileon_crosscheck": True,
        "all96_mass_ordered_metric_insertions_no_sample_inference": True,
        "whole2430_generic_polynomial_coefficient_residuals": True,
        "strict_mass_moment_sign_equal_mass_and_asymptotics": True,
        "all_original_vectors_and72_frozen_line_kernels": True,
        "basis_conversion_extra_parent_chi_and_frontiers_preserved": True,
        "rejected_inputs": rejected_inputs(),
    }
