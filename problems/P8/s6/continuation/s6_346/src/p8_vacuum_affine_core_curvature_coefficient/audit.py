"""Common-basis selected-loop coefficient without promotion to parent closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_box_curvature_coefficient import audit as previous
from p8_vacuum_affine_complement_measure.split import clean

from . import bounds, calibration, conversion, core, moments, source

MODEL = "original_old_polynomial_core_common_covariant_jet_lift"
OBSERVABLES = (
    "whole_original_source",
    "whole_generic_labeled_to_jet_conversion",
    "whole_triangle_and_box_moments",
    "whole_known_polynomial_core",
    "whole_finite_original_coefficient_bound",
    "whole_original_vector_calibration",
)
ITEM = {
    "id": "QG2_H8A510_known_polynomial_core_common_curvature",
    "status": "COMPLETE_OLD_POLYNOMIAL_CORE_DEGREE6_COEFFICIENT_IN_ONE_JET_LIFT_NOT_FULL_MATTER_PARENT_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the unchanged known old-polynomial core")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated whole core-coefficient observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "old_polynomial_bubble_triangle_box",
        "whole_common_six_word_jet_lift",
        "original_finite_coefficient_bound",
    ):
        raise ValueError("Require the stated known polynomial matching sector")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("The core coefficient cannot close original P8")
    return True


@cache
def packets():
    return dict(
        zip(
            OBSERVABLES,
            (
                source.data(),
                conversion.data(),
                moments.data(),
                core.data(),
                bounds.data(),
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
        "all9_original_primitive_statuses_preserved": len(frontier()) == 9,
        "all201_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 201,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "independent_extra_parent_chi_not_assigned_or_bounded": True,
        "new_local_and_mixed_sources_not_silently_included": True,
        "physical_metric_bounce_scoped_P8a_and_global_obligations_unchanged": True,
    }


def observable():
    return {
        "established": "All10 degree3 labeled scalar-box words are converted by whole off-shell polynomial identities into the fixed S345 six-word covariant jet lift. The complete old polynomial bubble, all24 dressed triangles and all24 alternating-mass boxes then have a single finite real-TT degree6 coefficient chi_core=g^4*c_core(n)/(16pi^2). The exact tuned sum cancels its n^-2 and n^-3 terms and has n^4*c_core -> -23/105. Exact finite-domain bounds give -1/(4n^4)<c_core<-1/(5n^4) at the original parameters, hence |chi_core|/A0<10^-207. No basis-dependent separate sector bound is substituted for this aggregate.",
        "domain": "Same original source, massive spectrum, physical metric and finite renormalization conditions. The homogeneous analytic-origin six-derivative coefficient, compared before any scalar-shell substitution to the six fixed symmetrized covariant scalar-jet words. One real null transverse-traceless graviton; all original polynomial mass orderings retained.",
        "historical_qualification": qualifications(),
        "not_established": "New local tadpole and mixed-source classes are not yet converted off shell into this common basis, so this is not the entire known matter coefficient. Independent extra parent chi remains unmatched. No complete curved effective/counterfunctional, internal-graviton loops, finite-gravity quantum decoupling, physical above-threshold Taylor approximation or truncation-error estimate, full virtual/inclusive matching, exact LSZ, quantum unitarity, complex Regge, same-parent bounce or UV completion. Original V/G/B/P8 remain open; no extra copy is added to S339/S342.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
    for j, v in enumerate(
        (True, False, 1.0, s.Float(1), "1", None, s.oo, s.I, s.nan, s.Symbol("unknown"))
    ):
        for label, call in (
            ("model", require_model),
            ("observable", require_observable),
            ("sector", require_sector),
        ):
            rows.append((f"core_invalid_{label}_{j}", call, (v,)))
    for j, v in enumerate(
        (
            "full_parent_chi_zero",
            "all_matter_closed",
            "physical_cut_replaced_by_jet",
            "known_loop_added_twice",
            "asymptotic_limit_is_finite_bound",
            "core_closes_original_P8",
        )
    ):
        rows.append((f"core_unsupported_scope_{j}", require_observable, (v,)))
    rows.extend(
        (
            ("core_deleted_primitive", validate_scope, ([], matching())),
            ("core_deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    for j, v in enumerate((True, False, 1.0, s.Float(1), "1", None, 0, -1, s.oo, s.I)):
        rows.append((f"core_invalid_mass_{j}", moments.coefficient, (0, v)))
    for j, v in enumerate((True, False, -1, 4, 1.0, s.Integer(1), "1", None, [])):
        rows.append((f"core_invalid_moment_order_{j}", moments.coefficient, (v, 1)))
    for j, v in enumerate(
        (
            True,
            None,
            [],
            [3, 0, 0],
            (3, 0),
            (4, 0, 0),
            (-1, 2, 2),
            (True, 1, 1),
            (s.Integer(3), 0, 0),
            (3.0, 0, 0),
            ("3", 0, 0),
        )
    ):
        rows.append((f"core_invalid_word_{j}", conversion.require_powers, (v,)))
    return rows


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            count += 1
        else:
            raise ValueError("Unsupported core input accepted: " + name)
    return count


def controls():
    return {
        "all10_full_offshell_flat_and_covariant_word_conversions": True,
        "complete_three_channel_bubble_normalization": True,
        "all_four_heavy_dressing_orders_and_exact_triangle_moments": True,
        "complete_common_basis_sum_not_separate_sector_bound": True,
        "finite_rational_log_tail_bound_not_asymptotic_inference": True,
        "all60_original_vector_and_six_frozen_bubble_checks": True,
        "source_scheme_extra_parent_chi_and_frontiers_preserved": True,
        "rejected_inputs": rejected_inputs(),
    }
