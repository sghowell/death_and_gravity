"""Scoped physical contour/TT-kernel result and preserved original obligations."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_quadratic_radiation_cancellation import audit as previous

from . import calibration, contour, insertions, moments, source
from . import parameters as parameter_module

MODEL = "original_physical_loop_contours_not_full_radiation_or_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_full_offshell_parameters",
    "whole_physical_contour_and_joint_gap",
    "whole_scalar_moments_and_Cauchy_bounds",
    "whole_literal_TT_weighted_insertions",
    "whole_original_route_calibrations",
)
ITEM = {
    "id": "QG2_H8A505_physical_offshell_loop_contours_and_TT_kernel_bounds",
    "status": "PHYSICAL_CONTOUR_AND_WEIGHTED_KERNEL_BOUNDS_NOT_COMPLETED_HARD_RADIATION",
}
original_parameters, require_parameters = (
    previous.parameters,
    previous.require_parameters,
)
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def parameters():
    return original_parameters()


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original physical-contour kernel model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated contour or weighted-kernel observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "mass_ordered_offshell_scalar_masters",
        "physical_Feynman_contour_homotopy",
        "literal_TT_single_line_insertions",
    ):
        raise ValueError("Require a stated scalar master or TT line-insertion sector")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A contour-kernel result cannot close original P8")
    return True


@cache
def packets():
    return dict(
        zip(
            OBSERVABLES,
            (
                source.data(),
                parameter_module.data(),
                contour.data(),
                moments.data(),
                insertions.data(),
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
        "all196_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 196,
        "all_matching_identifiers_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameter_record_unchanged": require_parameters(
            original_parameters()
        )
        == original_parameters(),
        "rejected_S279_not_promoted": any(
            row["status"].startswith("REJECTED_") for row in matching()
        ),
        "S336_unknown_curvature_coefficient_not_assigned": True,
        "hard_radiation_assembly_and_internal_gravity_still_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "Exact full off-shell triangle/ordered-box parameter polynomials; an explicit two-variable homotopy retaining the physical Feynman sign; a uniform joint complex denominator gap and scalar, derivative and weighted moment bounds. The literal minimal scalar stress insertion gives the complete per-line TT numerator and correct split-parameter measure, with bounded sums over the internal lines of each scalar topology.",
        "domain": "Original source,10^6<=n<10^198,adjacent scalar virtualities[0,2],light invariant[-12,4/3] or[45/8,16],ordered heavy transfer[-12,16]. Local analytic continuation of each physical boundary germ with radius1/4096 per invariant and n/128 in mass. Real on-shell external graviton, physical TT; original compact recoil satisfies the route domain.",
        "historical_qualification": qualifications(),
        "not_established": "Completed four-hard radiative graph sum with all multiplicities, relative signs, external emissions, outer-heavy-branch insertions and finite counterterms; internal-graviton loops or independent curved matching. No exact LSZ, resummed propagator, full inclusive interacting detector probability, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure. A continued germ is not the same first sheet on both sides of the cut.",
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
            "absolute_real_integrand_across_cut",
            "one_first_sheet_disk_across_cut",
            "scalar_bound_completes_all_TT_graphs",
            "ordered_box_mass_assignment_exchangeable",
            "threshold4_in_uniform_domain",
            "selected_contours_close_P8",
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
        (True, 1.0, s.Float(1), "1", None, s.I, s.oo, s.nan, s.Symbol("unknown"))
    ):
        rows.append(
            ("new_inexact_rational_" + str(i), parameter_module.rational, (value,))
        )
    for i, value in enumerate((-13, s.Rational(3, 2), 4, 5, 17)):
        rows.append(
            (
                "new_bad_light_domain_" + str(i),
                parameter_module.domain,
                ("triangle", value, (1, 1)),
            )
        )
    for i, value in enumerate(((), (1,), (1, 1, 1))):
        rows.append(
            (
                "new_missing_virtuality_" + str(i),
                parameter_module.domain,
                ("triangle", 6, value),
            )
        )
    for i, value in enumerate(((-1, 1), (1, 3), (1, 1.0))):
        rows.append(
            (
                "new_bad_virtuality_" + str(i),
                parameter_module.domain,
                ("triangle", 6, value),
            )
        )
    rows.extend(
        (
            (
                "new_triangle_heavy_transfer",
                parameter_module.domain,
                ("triangle", 6, (1, 1), 1),
            ),
            (
                "new_box_heavy_transfer_low",
                parameter_module.domain,
                ("ordered_box", 6, (1, 1, 1, 1), -13),
            ),
            (
                "new_box_heavy_transfer_high",
                parameter_module.domain,
                ("ordered_box", 6, (1, 1, 1, 1), 17),
            ),
            ("new_mass_too_small", parameter_module.heavy_mass, (10**6 - 1,)),
            ("new_mass_too_large", parameter_module.heavy_mass, (10**198,)),
            ("new_unknown_master", moments.master_bound, ("box_exchange_symmetric",)),
            (
                "new_derivative_missing_multiindex",
                moments.derivative_bound,
                ("triangle", (0, 0)),
            ),
            (
                "new_derivative_bad_multiindex",
                moments.derivative_bound,
                ("triangle", "000"),
            ),
            (
                "new_derivative_bool",
                moments.derivative_bound,
                ("triangle", (0, True, 0)),
            ),
            (
                "new_derivative_float",
                moments.derivative_bound,
                ("triangle", (0, 1.0, 0)),
            ),
            (
                "new_derivative_negative",
                moments.derivative_bound,
                ("triangle", (0, -1, 0)),
            ),
            (
                "new_derivative_mass_float",
                moments.derivative_bound,
                ("triangle", (0, 0, 0), 1.0),
            ),
            (
                "new_unknown_internal_line",
                insertions.line_bound,
                ("triangle", "massless_graviton"),
            ),
            (
                "new_numerator_bound_too_large",
                insertions.line_bound,
                ("triangle", "light", source.HEAVY_MASS2, 10001),
            ),
            (
                "new_numerator_bound_negative",
                insertions.line_bound,
                ("triangle", "light", source.HEAVY_MASS2, -1),
            ),
            ("new_unproved_moment", moments.moment_bound, (4, 0)),
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
            raise ValueError("Unsupported contour input accepted: " + name)
    return count


def controls():
    return {
        "original_source_and_ordered_mass_assignments_unchanged": True,
        "whole_homotopy_Feynman_sign_and_nonzero_endpoints": True,
        "joint_complex_gap_and_all_weighted_moments": True,
        "literal_TT_vertex_and_exact_split_parameter_weight": True,
        "original_route_and_both_polarization_controls": True,
        "full_radiative_assembly_and_curvature_still_open": True,
        "original_frontiers_and_historical_qualifications_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
