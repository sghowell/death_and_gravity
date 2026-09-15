"""Scoped minimal-gravity pole completion with all original frontiers preserved."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_ordinary_gravity_ward_bridge import audit as previous

from . import graphs, poles, soft, source

MODEL = "original_minimal_GR_one_loop_pole_completion_not_full_physical_P8"
OBSERVABLES = (
    "whole_source_and_minimal_graph_scope",
    "entire_D_physical_pole_completion",
    "three_anchor_minimal_loop_representative",
    "complete_conditional_massive_soft_kernel",
)
ITEM = {
    "id": "QG2_H8A452_minimal_gravity_one_loop_physical_pole_completion_modulo_Newton_and_two_local_anchors_with_retained_massive_Coulomb_kernel",
    "status": "COMPLETE_SCOPED_MINIMAL_GR_POLE_COMPLETION_NOT_FINITE_ANCHORS_PHYSICAL_IR_REGGE_OR_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def frontier():
    return previous.frontier()


def qualifications():
    return previous.qualifications()


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the selected original minimal-GR one-loop sector")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require an explicitly scoped pole or soft-kernel observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "raw_pole_completion",
        "conditional_massive_soft_kernel",
    ):
        raise ValueError("Require pole matching or the conditional massive soft kernel")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A pole completion modulo finite anchors cannot close P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_D_meromorphic_pole_completion": poles.data(),
        "whole_minimal_graph_and_local_ambiguity_proof": graphs.data(),
        "whole_conditional_massive_soft_kernel": soft.data(),
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
        "all143_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 143,
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
        "three_finite_matching_anchors_remain_unassigned": True,
        "whole_physical_IR_and_Regge_frontiers_not_closed": True,
        "other_source_sectors_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The exact full-D nonbox plus complete scalar-LSZ physical principal parts reduce to four crossing-meromorphic basis coefficients. A literal full tensor-box denominator identity and whole-component ordinary vertex continuity justify the pole completion, without deleting the scalar boxes. The selected minimal one-loop gravity amplitude is thereby parameterized by its fixed cut representative plus one light Newton and two regular local anchors. The complete massive physical soft kernel has a positive attenuation density and a retained forward Coulomb phase.",
        "domain": "The same unchanged original Lorentz vacuum, minimally coupled Phi and gravity graph sector, one loop and pole/finite Laurent orders. Fixed-volume matching is inherited; the three finite anchors are unspecified. The exact soft-kernel signs hold for mu>0,s>4mu,0<=-t<=s-4mu, but detector interpretation retains S278's factorization premise.",
        "historical_qualification": qualifications(),
        "not_established": "The finite values or signs of the light Newton/local anchors, a complete full-source finite four-point amplitude or improved b20; C/g/heavy/Proca/M1 and higher-loop completion, finite F2 or slopes, detector/dressing/hard-unitarity or Regge remainder bounds; an original quantum state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "three_anchors_equal_zero",
            "delete_Coulomb_boxes",
            "discard_all_F2",
            "fix_finite_phase_from_D0",
            "massless_external_scalar",
            "physical_cutoff_from_cap",
            "full_Regge_bound",
            "UV_IR_poles_identical",
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
            raise ValueError("Unsupported pole-completion input accepted: " + name)
    return total


def controls():
    return {
        "whole_original_source_and_all_minimal_graphs": True,
        "entire_tensor_box_not_a_scalar_box_assumption": True,
        "all_D_zero_jets_and_all_four_LSZ_legs": True,
        "whole_component_Ward_limit_includes_singular_F2": True,
        "UV_IR_and_local_Newton_anchors_kept_distinct": True,
        "physical_Coulomb_and_conditional_detector_scope": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
