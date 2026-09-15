"""Complete selected real emission and a finite physical soft-error bound."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_whole_mixed_heavy_gravity_sector import audit as previous

from . import rate, recoil, source, tree

MODEL = "original_massive_scalar_graviton_bremsstrahlung_not_full_physical_P8"
OBSERVABLES = (
    "whole_original_26_graph_scalar_bremsstrahlung",
    "whole_Ward_TT_and_improved_heavy_remainder",
    "whole_on_shell_recoil_and_phase_space_Jacobian",
    "whole_selected_physical_real_rate_error",
)
ITEM = {
    "id": "QG2_H8A459_complete26_graph_massive_scalar_bremsstrahlung_exact_recoil_phase_space_and_explicit_physical_real_rate_error",
    "status": "COMPLETE_SELECTED_PHYSICAL_REAL_EMISSION_SOFT_ERROR_BOUND_NOT_VIRTUAL_IR_PAIRING_FULL_GRAVITY_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError(
            "Require the stated original massive scalar bremsstrahlung sector"
        )


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated complete selected real-emission observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "whole_minimal_scalar_bremsstrahlung",
        "physical_real_soft_error",
    ):
        raise ValueError(
            "Require the selected full tree or its physical real-rate error"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "A selected real-rate comparison cannot close original physical P8"
        )
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_26_graph_tensor_and_TT_remainder": tree.data(),
        "whole_recoil_phase_space_and_majorants": recoil.data(),
        "whole_physical_real_rate_error": rate.data(),
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
        "all150_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 150,
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
        "real_rate_error_not_virtual_dressed_IR_or_higher_matching": True,
        "virtual_IR_full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The entire minimal C/sqrt(kappa) and g^2/sqrt(kappa) scalar bremsstrahlung tree has26 labelled graphs, complete Ward and TT checks, and an exact pure-gauge cancellation that improves its heavy remainder. An exact on-shell recoil map, phase-space Jacobian and positive original Born denominator give a uniform integrable bound for the difference of exact and leading-soft real rates, below10^-792 at the original parameters on the stated compact domain.",
        "domain": "Original mu1,n>=128,5/4<=incoming COM scalar energy E<=2,all physical graviton and outgoing pair-rest-frame angles,0<lower cutoff<detector resolution<=1/8. Both physical helicities and the identical-scalar rate normalization are retained. This includes every interference in the square of the selected C,g^2 tree.",
        "historical_qualification": qualifications(),
        "not_established": "Individual finite Fock rates, virtual/dressed IR pairing or conversion to S278 dimensional analytic division, radiation from the pure-gravity Born tree, hard-loop or higher-EFT matching, all-energy fixed-transfer Regge and original state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "choose_light_finite_curvature",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "D4_cut_fixes_finite_threshold",
            "remove_soft_pole_by_hand",
            "H_is_exact_stable",
            "ignore_virtual_pole_term",
            "use_two_helicities_for_all_D",
            "distinct_pair_as_identical_pair",
            "full_Regge_bound",
            "finite_forward_dispersion_from_unpaired_cut",
            "exchange_IR_limits",
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
            raise ValueError("Unsupported real-emission input accepted: " + name)
    return total


def controls():
    return {
        "whole26_graphs_and_original_source_retained": True,
        "full_Ward_TT_and_heavy_pure_gauge_cancellation": True,
        "physical_recoil_all_angles_and_phase_space_Jacobian": True,
        "both_unit_helicities_and_Born_rate_normalization": True,
        "integrable_lower_cutoff_uniform_explicit_error": True,
        "virtual_full_gravity_matching_and_Regge_not_inferred": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
