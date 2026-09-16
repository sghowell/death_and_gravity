"""Selected physical inclusive reference conversion, not full P8 closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import audit as previous

from . import continuity, inclusive, soft, source

MODEL = "original_physical_virtual_soft_pairing_not_full_P8"
OBSERVABLES = (
    "whole_original_selected_real_virtual_coupling_inventory",
    "whole_dimensional_soft_kernel_and_finite_conversion",
    "whole_fixed_domain_dimensional_recoil_continuity",
    "whole_selected_inclusive_reference_and_explicit_error",
)
ITEM = {
    "id": "QG2_H8A460_selected_gravitational_dressing_real_virtual_pairing_finite_analytic_reference_conversion_and_explicit_error",
    "status": "COMPLETE_SELECTED_ONE_NEWTON_PHYSICAL_RATE_CONVERSION_NOT_FULL_MATCHING_GRAVITY_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated selected physical real-virtual sector")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError(
            "Require a stated selected physical inclusive-reference observable"
        )


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "selected_matter_inclusive_rate",
        "finite_physical_analytic_conversion",
    ):
        raise ValueError(
            "Require the selected inclusive rate or physical-analytic conversion"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "A selected inclusive-reference comparison cannot close original P8"
        )
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_dimensional_soft_kernel_and_conversion": soft.data(),
        "whole_fixed_domain_continuity_and_majorants": continuity.data(),
        "whole_selected_inclusive_reference": inclusive.data(),
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
        "all151_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 151,
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
        "selected_rate_conversion_not_full_IR_or_higher_matching": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The selected gravitational-dressing contribution to the matter-Born rate has an explicit finite conversion to the exact inherited analytic soft-divided hard reference. All-D angular and phase factors, moving diagonal terms and the physical Coulomb sheet are retained. A fixed-domain dominated-convergence proof gives the physical S295 nonsoft remainder. The total conversion/recoil error is below10^-792 at original parameters on the stated compact domain.",
        "domain": "Original mu1,n>=128,5/4<=incoming COM scalar energy<=2,all physical hard and graviton angles,0<detector resolution<=1/8. The dimensional regulator is removed at fixed resolution before any further limit. The reference scale for the explicit regulator majorant is nu1.",
        "historical_qualification": qualifications(),
        "not_established": "The independent finite hard matching coefficient or its sign, pure-gravity-exchange radiation, hard-loop and higher-operator corrections, nonperturbative/dressed all-channel unitarity, a uniform forward gravitational or all-energy Regge limit, or original state/domain/measure/bounce and V/G/B/P8 closure.",
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
            "full_inclusive_rate_at_same_Newton_order",
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
        "whole_selected_virtual_and26_graph_real_sectors": True,
        "full_D_soft_phase_and_moving_diagonal_terms": True,
        "physical_Coulomb_sheet_and_local_UV_subtraction": True,
        "fixed_domain_positive_continuation_and_recoil_Jacobian": True,
        "explicit_finite_conversion_and_regulator_majorants": True,
        "hard_matching_and_full_gravity_Regge_not_inferred": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
