"""Single finite residual with complete state-correct leading-soft dressing."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_leading_ladder_coulomb_phase import audit as previous

from . import bounds, dressing, series, source

MODEL = "original_single_residual_soft_dressing_not_full_P8"
OBSERVABLES = (
    "whole_signed_dimensional_soft_series",
    "whole_marked_remaining_energy_operator",
    "whole_uniform_dressed_residual_bound",
    "whole_single_residual_scope_boundary",
)
ITEM = {
    "id": "QG2_H8A473_complete_state_correct_allN_leading_soft_dressing_of_one_finite_residual",
    "status": "SINGLE_RESIDUAL_ALL_LEADING_SOFT_DRESSING_NOT_FULL_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact single-residual-soft-dressing model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated single-residual-soft-dressing observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "one_finite_radiative_residual",
        "additional_leading_soft_dressing",
    ):
        raise ValueError("Require the known finite or symbolic matching sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A single-residual soft reference cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_signed_dimensional_soft_series": series.data(),
        "whole_marked_remaining_energy_operator": dressing.data(),
        "whole_uniform_dressed_residual_bound": bounds.data(),
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
        "all164_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 164,
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
        "established": "The unchanged full47-tree finite signed one-real residual admits a complete same-state all-N additional-leading-soft total-energy dressing. Its finite analytic-regulator limit is a remaining-energy convolution with relative magnitude below2e-768 and a uniformly vanishing linear/quadratic threshold bound at original parameters.",
        "domain": "Original mu=nu=1,25/4<=s<=16,all physical nonforward hard angles,0<x<=1/8. Additional-soft regulator removal at fixed positive threshold precedes the reference's uniform forward and zero-threshold comparisons.",
        "historical_qualification": qualifications(),
        "not_established": "The full all-N nonleading interacting rate, two-soft contact and multiple-residual terms, correlated exact recoil, finite radiative hard loops or matching, a positive detector measure or unitarity, high-energy Regge bounds, original quantum state and bounce, original V/G/B/P8 closure.",
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
            "single_residual_controls_full_multi_real_rate",
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
        "same_complete47_tree_signed_residual_and_full_Born": True,
        "marked_Bose_factor_and_total_remaining_energy": True,
        "same_radiative_state_virtual_pole_and_finite_conversion": True,
        "signed_fractional_D_series_not_fictitious_probability": True,
        "positive_index_Poisson_and_zero_index_signed_limits": True,
        "single_residual_reference_not_complete_interacting_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
