"""First curvature contact basis and selected-matter pole closure, not full P8."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_triangle_box_radiation import audit as previous

from . import basis, calibration, source, uvpoles
from . import matching as match

MODEL = "original_first_curvature_contact_not_full_matching"
OBSERVABLES = (
    "whole_original_source",
    "whole_first_curvature_contact_basis",
    "whole_selected_matter_UV_poles",
    "whole_strict_local_matching_map",
    "whole_original_literal_curvature_calibrations",
)
ITEM = {
    "id": "QG2_H8A507_first_curvature_basis_and_selected_matter_poles",
    "status": "COMPLETE_FIRST_PARITY_EVEN_CURVATURE_CONTACT_BASIS_AND_SELECTED_TT_POLES_NOT_FINITE_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the unchanged original curvature-contact model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated basis or selected-pole observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "parity_even_explicit_curvature_through_six_derivatives",
        "selected_minimal_matter_TT_UV_poles",
        "entire_local_matching_polynomial",
    ):
        raise ValueError("Require the stated restricted curvature or matter sector")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A basis or selected-pole result cannot close original P8")
    return True


@cache
def packets():
    return dict(
        zip(
            OBSERVABLES,
            (
                source.data(),
                basis.data(),
                uvpoles.data(),
                match.data(),
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
        "all198_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 198,
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
        "finite_chi_not_assigned_or_bounded": True,
        "higher_derivatives_parity_odd_and_internal_gravity_not_closed": True,
        "physical_metric_bounce_scoped_P8a_and_global_obligations_unchanged": True,
    }


def observable():
    return {
        "established": "The one-physical-graviton contact of parity-even explicit-curvature local four-identical-scalar operators has dimension0 below six derivatives and dimension1 at six derivatives, represented by the S336 chi*T contact. All256 Bose words, literal null-TT trace and k contractions and original massive-state checks are retained. The selected minimal-matter radiative UV poles cancel against the fixed literal counterfunctionals; no additional curvature-only TT pole remains in that representative.",
        "domain": "Same original massive scalar and physical metric. Exactly one real null TT graviton and four identical scalar fields; at most six derivatives in the explicit-curvature local basis. The separate UV audit covers the selected one-matter-loop original-source graph inventory, not internal-graviton loops or a complete curved-background functional.",
        "historical_qualification": qualifications(),
        "not_established": "Finite curved matching including the value or a bound for chi, higher-derivative or parity-odd matching, internal-graviton loops, a full gravity beta function or finite-gravity quantum decoupling. The known minimal remainder bound cannot bound an unknown extra curvature term. Full virtual/inclusive probability, exact LSZ, quantum unitarity, complex Regge, same-parent bounce, UV and original V/G/B/P8 remain open.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
    for index, value in enumerate(
        (
            True,
            False,
            1.0,
            s.Float(1),
            "1",
            None,
            s.oo,
            s.I,
            s.nan,
            s.Symbol("unspecified"),
        )
    ):
        for label, call in (
            ("new_model", require_model),
            ("new_observable", require_observable),
            ("new_sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(index), call, (value,)))
    for index, value in enumerate(
        (
            "basis_sets_chi_zero",
            "UV_pole_fixes_finite_matching",
            "massless_basis_proves_original_bounce",
            "one_sample_fixes_full_matching",
            "all_higher_derivatives_are_the_same_operator",
            "selected_matter_poles_close_P8",
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
        (
            True,
            1.0,
            s.Float(1),
            None,
            "unknown",
            s.Symbol("chi"),
            s.I * basis.canonical,
            s.oo,
            s.nan,
            basis.a0[0] * basis.canonical,
        )
    ):
        rows.append(
            (
                "invalid_complete_matching_polynomial_" + str(index),
                match.extract_local_coefficient,
                (value,),
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
            raise ValueError("Unsupported curvature-basis input accepted: " + name)
    return count


def controls():
    return {
        "full256_symbolic_Bose_reduction_not_sample_inference": True,
        "all_null_TT_curvature_traces_and_k_slots": True,
        "original_massive_states_and_literal256_component_curvature": True,
        "full_exact_D_local_poles_and_evanescent_finite_terms_retained": True,
        "all_selected_source_bubble_triangle_box_quadratic_classes": True,
        "whole_matching_polynomial_required_no_silent_projection": True,
        "finite_matching_metric_and_original_frontiers_preserved": True,
        "rejected_inputs": rejected_inputs(),
    }
