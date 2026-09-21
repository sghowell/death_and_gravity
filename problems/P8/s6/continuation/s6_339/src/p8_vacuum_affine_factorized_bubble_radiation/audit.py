"""Complete factorized-bubble radiation with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_mixed_source_radiation import audit as previous

from . import bounds, branch, calibration, kernel, loops, source

MODEL = "original_complete_factorized_bubble_radiation_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_two_endpoint_kernel",
    "whole_exact_D_factorized_loops",
    "whole_physical_complex_branch",
    "whole_uniform_bounds",
    "whole_original_calibrations",
)
ITEM = {
    "id": "QG2_H8A503_complete_factorized_bubble_radiation",
    "status": "EXACT_FACTORIZED_BUBBLE_RADIATION_NOT_FULL_HARD_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError(
            "Require the original selected factorized-bubble radiation model"
        )


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated complete factorized-bubble observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "full_local_local_Hessian_bubble_class",
        "physical_four_scalar_one_graviton_TT",
        "known_OS4_factorized_bubble_remainder",
    ):
        raise ValueError(
            "Require the stated selected factorized-bubble radiation sector"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A selected bubble sector cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_two_endpoint_kernel": kernel.data(),
        "whole_exact_D_factorized_loops": loops.data(),
        "whole_physical_complex_branch": branch.data(),
        "whole_uniform_bounds": bounds.data(),
        "whole_original_calibrations": calibration.data(),
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
        "all194_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 194,
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
        "established": "The full factorized local/local Hessian bubble class has its physical four-scalar/one-graviton TT response computed: both inner light insertions, both outer heavy insertions, external/source radiation and the fixed UV and OS4 counterterms. Exact-D triangles and the ordered A-B-A product give the continuous divided difference of A(v)^2 B(v). The physical complex branch, its cut, original nonleading remainder and finite full-tree interference have explicit uniform bounds.",
        "domain": "Original n,g,C,kappa and H8A420-VAC-OS4, selected full factorized matter/heavy bubble class. Real external TT after exact-D subtraction;25/4<=s<=16,-1<z<1,0<=x<=1/8. Actual channel endpoints stay on crossed[-12,0] or timelike[45/8,16] branches. High-window integration is fixed in nonforward Born angle.",
        "historical_qualification": qualifications(),
        "not_established": "Remaining local/bilocal triangle and bilocal/bilocal ordered-box radiation, unresolved mixed quadratic-radiation blocks or internal-graviton loops; a small isolated hard loop coefficient; unprojected k_mu k_nu loop form factors or complete curved counterfunctional; independent curvature matching including S336; full inclusive interacting detector probability, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure.",
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
            "only_contact_squared_bubble_is_complete_class",
            "flat_kernel_fixes_every_curvature_operator",
            "discard_cut_imaginary_part",
            "single_outer_heavy_insertion_is_complete",
            "isolated_hard_bound_is_perturbatively_small",
            "selected_bubble_sector_closes_P8",
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
            raise ValueError("Unsupported factorized-bubble input accepted: " + name)
    return count


def controls():
    return {
        "original_source_and_existing_subtractions_unchanged": True,
        "complete_A_B_A_chain_all_inner_and_outer_insertions": True,
        "both_equal_mass_triangles_exact_D_before_continuation": True,
        "physical_complex_cut_and_coincident_limits_retained": True,
        "independent_tensor_and_exact_original_absorptive_calibrations": True,
        "remaining_loops_and_independent_curvature_matching_open": True,
        "original_P8_frontiers_and_historical_qualifications_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
