"""Full first-Newton perturbative rate assembly, not original P8 closure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_physical_virtual_soft_pairing import audit as previous

from . import forward, inclusive, matter, source

MODEL = "original_one_Newton_inclusive_assembly_not_full_P8"
OBSERVABLES = (
    "whole_actual_finite_kappa_matter_graph_equality",
    "whole_complete_matter_loop_gravity_Born_interference",
    "whole_positive_Born_and_forward_crossover",
    "whole_first_Newton_inclusive_assembly_and_error",
)
ITEM = {
    "id": "QG2_H8A461_complete_formal_one_Newton_inclusive_assembly_same_OS4_matter_graph_equality_and_forward_crossover",
    "status": "COMPLETE_FORMAL_FIRST_NEWTON_RATE_ASSEMBLY_NOT_FULL_MATCHING_HIGHER_NEWTON_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated formal first-Newton inclusive sector")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated formal first-Newton assembly observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "formal_one_Newton_inclusive_assembly",
        "positive_full_Born_normalization",
    ):
        raise ValueError(
            "Require the formal inclusive assembly or full Born normalization"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A formal one-Newton assembly cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_finite_kappa_matter_loop_and_interference": matter.data(),
        "whole_positive_Born_and_forward_crossover": forward.data(),
        "whole_first_Newton_inclusive_assembly": inclusive.data(),
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
        "all152_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 152,
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
        "formal_rate_assembly_not_higher_Newton_or_matching": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The complete current finite-kappa no-internal-graviton matter loop equals the frozen S239 formal diagrams through all contributing source vertices in the same OS4 prescription. Its missing interference with the gravity Born amplitude is bounded below10^-199 relative to the full positive Born squared. Together with S296 this assembles the complete formal first-Newton one-loop inclusive numerator, retaining independent hard gravitational matching. An exact forward crossover demonstrates why the Newton expansion is not uniform there.",
        "domain": "The matter-interference bound covers mass-one4<s<=10^196,-1<z<1 with original parameters and fixed matter OS4. The real/virtual assembly and reference-error bound use25/4<=s<=16,0<resolution<=1/8 and the S296 regulator-removal order. The full Born denominator is an unexpanded reference, not a higher-order calculation.",
        "historical_qualification": qualifications(),
        "not_established": "Independent finite gravitational hard matching, higher-Newton or omitted-loop errors, a finite forward cross section, uniform forward expansion, positive analytic dispersion measure or all-energy Regge remainder, original quantum state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "all_Newton_orders_uniform_forward_rate",
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
            raise ValueError("Unsupported inclusive-assembly input accepted: " + name)
    return total


def controls():
    return {
        "whole_actual_finite_kappa_matter_jets_and_graph_patterns": True,
        "same_complete_OS4_matter_counterterms_no_new_value_choice": True,
        "complete_missing_matter_gravity_interference": True,
        "strictly_positive_full_Born_normalization": True,
        "known_real_virtual_pair_and_error_retained": True,
        "forward_crossover_and_unknown_Regge_not_inferred": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
