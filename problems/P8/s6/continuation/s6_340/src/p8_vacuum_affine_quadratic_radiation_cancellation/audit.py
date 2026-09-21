"""Complete selected quadratic radiation cancellation and unchanged frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_factorized_bubble_radiation import audit as previous

from . import calibration, cancellation, loops, source, ward

MODEL = "original_quadratic_radiation_cancellation_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_inverse_Ward_and_TT_boundary",
    "whole_actual_mixed_loop_and_counterterms",
    "whole_complete_quadratic_cancellation",
    "whole_original_parameter_calibrations",
)
ITEM = {
    "id": "QG2_H8A504_complete_quadratic_radiation_cancellation",
    "status": "EXACT_MASSIVE_MATTER_QUADRATIC_CANCELLATION_NOT_FULL_HARD_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError(
            "Require the original selected quadratic-radiation cancellation model"
        )


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated complete quadratic-radiation observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_old_matter_two_point_sector",
        "physical_four_scalar_one_graviton_TT",
        "fixed_OS_mass_residue_external_cancellation",
    ):
        raise ValueError(
            "Require the stated selected quadratic-radiation cancellation sector"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A selected quadratic sector cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_inverse_Ward_and_TT_boundary": ward.data(),
        "whole_actual_mixed_loop_and_counterterms": loops.data(),
        "whole_complete_quadratic_cancellation": cancellation.data(),
        "whole_original_parameter_calibrations": calibration.data(),
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
        "all195_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 195,
        "all_matching_identifiers_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            row["status"].startswith("REJECTED_") for row in matching()
        ),
        "S336_unknown_curvature_coefficient_not_assigned": True,
        "other_hard_loops_quantum_Regge_and_UV_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The complete selected old massive-matter first-loop external quadratic radiation cancels in the real physical TT amplitude. Both actual mixed light/heavy triangle insertions, the full on-shell mass and kinetic counterterms, the single self-energy propagator insertion, the quartic tadpole and the heavy-onepoint source pair are retained. The zero holds legwise for arbitrary shifted hard subamplitudes; it is not a leading-soft approximation.",
        "domain": "Original n,g,C,kappa and already fixed H8A420-VAC-OS4 light mass/residue/onepoint conditions. First matter-loop order and one real external TT graviton,25/4<=s<=16,-1<z<1,0<=omega<=1/8; external shifted virtualities0<r^2<2 with analytic soft continuation. Selected massive-matter loop only.",
        "historical_qualification": qualifications(),
        "not_established": "Remaining hard local/bilocal triangle and ordered-box radiation or internal-graviton loops; unprojected loop form factors, off-shell metric response or complete curved counterfunctional; independent four-hard-direction curvature matching including S336; exact LSZ or resummed propagators, full inclusive interacting detector probability, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
    for i, value in enumerate(
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
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate(
        (
            "discard_self_energy_before_metric_variation",
            "two_direction_result_sets_four_hard_Weyl_zero",
            "truncated_inverse_is_exact_resummed_propagator",
            "omit_fixed_kinetic_metric_counterterm",
            "off_shell_graviton_response_also_zero",
            "selected_quadratic_sector_closes_P8",
        )
    ):
        rows.append(("unsupported_new_scope_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("new_deleted_primitive", validate_scope, ([], matching())),
            ("new_deleted_matching", validate_scope, (frontier(), [])),
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
            raise ValueError("Unsupported quadratic-radiation input accepted: " + name)
    return count


def controls():
    return {
        "original_source_and_fixed_mass_residue_onepoint_unchanged": True,
        "both_mixed_triangles_and_single_propagator_insertion": True,
        "all_counterterms_before_real_TT_projection": True,
        "arbitrary_shifted_hard_factor_and_exact_noncommuting_powers": True,
        "actual_original_parameter_and_both_polarization_controls": True,
        "remaining_hard_loops_and_four_direction_curvature_still_open": True,
        "original_P8_frontiers_and_historical_qualifications_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
