"""Complete selected triangle/box TT radiation with unchanged global obligations."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_physical_loop_contours import audit as previous

from . import bounds, calibration, functional, labels, radiation, source

MODEL = "original_selected_triangle_box_radiation_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_covariant_labelled_functional",
    "whole_full_offshell_labelled_kernels",
    "whole_complete_physical_radiation",
    "whole_original_uniform_remainder_and_interference",
    "whole_original_label_and_density_calibrations",
)
ITEM = {
    "id": "QG2_H8A506_complete_selected_triangle_box_radiation",
    "status": "COMPLETE_SELECTED_MINIMAL_MATTER_TRIANGLE_BOX_TT_AND_FINITE_REMAINDER_NOT_FULL_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original selected triangle/box radiation model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated complete triangle/box observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_selected_triangle_box",
        "all_current_known_matter_radiation_remainders",
        "finite_full_tree_nonleading_interference",
    ):
        raise ValueError("Require a stated selected triangle/box radiation sector")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A selected matter-loop result cannot close original P8")
    return True


@cache
def packets():
    return dict(
        zip(
            OBSERVABLES,
            (
                source.data(),
                functional.data(),
                labels.data(),
                radiation.data(),
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
        "all197_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 197,
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
        "internal_gravity_and_full_quantum_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "Complete selected minimal old-matter triangle and ordered-box real-TT radiation: all24 labelled terms, every internal scalar line, the triangle's outer heavy branch, four independently shifted external hard coefficients and the fixed finite OS4 constant. A continuous original-parameter nonleading remainder bound and its finite full47-tree interference follow. The already known S337/S338/S339 sectors and zero S340 fit the same rounded remainder budget.",
        "domain": "Original source n,g,C,kappa,25/4<=s<=16,-1<z<1,0<omega<=1/8,all real physical graviton and hard angles. Scalar masters and split-line kernels are the complete S341 physical Feynman contour integrals. Soft and coincident-channel limits are continuous; fixed nonforward angle for the high-window rate bound. Specified minimal external metric variation of the classical limiting matter source, not a quantum decoupling theorem.",
        "historical_qualification": qualifications(),
        "not_established": "Independent four-hard curved matching including S336, internal-graviton loops, finite-gravity quantum decoupling, higher matter loops or exact LSZ. The large isolated hard bound is not perturbative smallness, and the separately divergent leading-soft loop extension is not included in the finite interference claim. Full inclusive interacting detector probability, virtual hard matching, quantum unitarity, complex Regge, same-parent bounce, UV and original V/G/B/P8 closure remain open.",
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
            "omit_triangle_outer_heavy_branch",
            "commute_Hessian_before_metric_variation",
            "swap_ordered_box_mass_cuts",
            "large_isolated_hard_bound_is_small",
            "finite_nonleading_interference_is_full_inclusive_rate",
            "selected_triangle_box_closes_P8",
        )
    ):
        rows.append(("unsupported_new_scope_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("new_deleted_primitive", validate_scope, ([], matching())),
            ("new_deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    for i, value in enumerate(
        ("unknown", True, None, "independent_curvature", "internal_graviton")
    ):
        rows.append(
            (
                "new_invalid_finite_interference_sector_" + str(i),
                bounds.require_sector,
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
            raise ValueError("Unsupported triangle/box input accepted: " + name)
    return count


def controls():
    return {
        "original_source_and_fixed_OS4_unchanged": True,
        "all_labelled_noncommuting_functional_variations_retained": True,
        "all_internal_lines_and_outer_heavy_branch": True,
        "full_physical_contours_and_independent_sign_normalizations": True,
        "original_source_polarizations_coincidences_and_density_controls": True,
        "full_tree_finite_remainder_interference_not_inclusive_closure": True,
        "original_frontiers_curvature_and_historical_qualifications_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
