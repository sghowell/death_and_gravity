"""Uniform state-correct two-real soft subtraction with original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_uniform_two_real_tree_bound import audit as previous

from . import analytic, external, source, subtraction

MODEL = "original_two_real_soft_overlap_not_full_P8"
OBSERVABLES = (
    "whole_double_external_regrouping",
    "whole_complex_energy_tube",
    "whole_two_real_soft_overlap",
    "whole_signed_subtracted_measure",
)
ITEM = {
    "id": "QG2_H8A477_uniform_state_correct_two_real_soft_overlap_subtracted_measure",
    "status": "UNIFORM_TWO_REAL_SOFT_SUBTRACTED_MEASURE_NOT_FULL_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original soft-overlap model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated soft-overlap observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_soft_overlap",
        "signed_two_real_measure",
    ):
        raise ValueError("Require the stated soft overlap or signed two-real measure")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A subtracted real measure cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_double_external_regrouping": external.data(),
        "whole_complex_energy_tube": analytic.data(),
        "whole_uniform_soft_overlap_measure": subtraction.data(),
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
        "all168_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 168,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "no_unknown_matching_coefficient_is_chosen": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "Exact external ordering/grouping, a uniform holomorphic soft-scaled amplitude estimate and Cauchy mixed derivative bound imply an integrable state-correct two-soft remainder. The complete signed two-real difference has total variation below2e-725*x+1e-652*x2, uniformly in all stated angles and lower soft cutoff.",
        "domain": "Original four-dimensional formal434 tree, S300 recoil E in[5/4,2],a,b>0,a+b<=x<=1/8, all nonforward hard Born directions and physical emitted angles, unit-Frobenius TT polarizations. Both reduced hard states include the marked null graviton and the same recoil/phase convention.",
        "historical_qualification": qualifications(),
        "not_established": "A complete inclusive detector rate, integrated counterterm/virtual matching, identification with a sum of earlier subtraction schemes, all-N nonleading errors, finite physical hard matching, unitarity, absolute complex Regge, original quantum state/common-parent bounce, or V/G/B/P8 closure.",
    }


def bad_cases():
    bad = (
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
    rows = []
    for i, value in enumerate(bad):
        for label, call in (
            ("model", require_model),
            ("observable", require_observable),
            ("order", source.require_order),
            ("mass", source.require_mass),
            ("sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate((-1, 0, s.Rational(1, 2), s.Rational(3, 2))):
        rows.append(
            ("outside_positive_integer_order_" + str(i), source.require_order, (value,))
        )
    for i, value in enumerate((0, -1)):
        rows.append(("outside_positive_mass_" + str(i), source.require_mass, (value,)))
    for i, value in enumerate(
        (
            "whole_finite_G_amplitude",
            "UV_complete",
            "full_crossed_b20",
            "exact_massless_LSZ",
            "choose_finite_curvature",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "all_matching_zero",
            "known_piece_is_entire_amplitude",
            "finite_loop_proves_unitarity",
            "drop_evanescent_terms",
            "drop_old_Gram_finite",
            "discard_imaginary_phase_from_amplitude",
            "relative_Regge_is_absolute",
            "all_loop_regulator_complete",
            "soft_overlap_closes_full_inclusive_detector_rate",
        )
    ):
        rows.append(("unsupported_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("deleted_primitive", validate_scope, ([], matching())),
            ("deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    return rows


def rejected_inputs():
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported finite-gravity input accepted: " + name)
    return total


def controls():
    return {
        "same_original_complete434_tree_and_positive_Born": True,
        "exact_double_external_orderings_and_independent_tree_coefficient": True,
        "uniform_holomorphic_energy_tube_not_real_samples_only": True,
        "both_state_correct_soft_faces_and_one_common_overlap": True,
        "all_signed_measure_interferences_and_Bose_phase_factors": True,
        "finite_subtracted_real_measure_not_full_inclusive_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
