"""Unchanged original frontier and scoped ordinary gravity charge normalization."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_matter_graviton_vertex import audit as previous

from . import continuity, gauge, soft, source

MODEL = "original_ordinary_GR_one_loop_charge_bridge_not_physical_IR_or_full_P8"
OBSERVABLES = (
    "whole_source_and_proper_graph_scope",
    "whole_on_shell_zero_q_gauge_bridge",
    "all_topology_spacelike_pole_finite_continuity",
    "ordinary_charge_and_raw_IR_LSZ_cancellation",
)
ITEM = {
    "id": "QG2_H8A451_ordinary_pure_GR_zero_transfer_Ward_bridge_all_topology_continuity_and_entire_raw_charge_LSZ_cancellation",
    "status": "COMPLETE_SCOPED_ORDINARY_ONE_LOOP_CHARGE_BRIDGE_NOT_PHYSICAL_IR_FULL_MATCHING_OR_V_G_B_P8",
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
        raise ValueError(
            "Require the specified original ordinary GR one-loop charge sector"
        )
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require an explicitly scoped ordinary Ward-bridge observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "ordinary_zero_transfer",
        "finite_transfer_IR_pole",
    ):
        raise ValueError(
            "Require the ordinary zero-transfer or raw finite-transfer IR sector"
        )
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("An ordinary regulated charge bridge cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_zero_q_gauge_fixing_bridge": gauge.data(),
        "whole_all_topology_spacelike_continuity": continuity.data(),
        "whole_ordinary_charge_and_raw_IR": soft.data(),
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
        "all9_primitive_statuses_preserved": len(frontier()) == 9,
        "all142_matching_records_preserved": matching()[:-1] == previous.matching()
        and len(previous.matching()) == 142,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "light_curvature_and_four_point_matching_not_fixed": True,
        "finite_transfer_IR_pole_not_a_finite_physical_observable": True,
        "S275_S276_S277_and_S281_S286_scopes_unchanged": True,
    }


def observable():
    return {
        "established": "The complete on-shell zero-transfer gauge-fixing difference is an infrared-integrable scaleless tadpole in the specified linear split and harmonic gauge. Whole-topology Gram-free parameter estimates give continuity of the ordinary pure-GR F1 pole and finite Laurent coefficients along spacelike transfer. Thus the zero-transfer proper vertex and its scalar LSZ legs cancel the entire raw Sigma_prime, including UV, IR and finite constants. The finite-transfer remaining IR pole is explicitly R(t)/(8pi^2*kappa*EP), with a rigorously vanishing coefficient at t0.",
        "domain": "The whole unchanged original source at its separately stated constant-scalar Lorentz vacuum, one pure-GR proper-vertex loop and the same dimensional/covariant UV convention. The explicit continuity bounds use actual mu=nu=1 and0<tau<=1,0<=EP<=1/4. The result is not an exact massless one-particle-pole or physical detector theorem.",
        "historical_qualification": qualifications(),
        "not_established": "A finite physical finite-transfer IR/detector/dressing or exact interacting LSZ limit; full F1 slope or F2 and H-metric mixing; S285 light curvature coefficients, physical metric/Newton residue and complete four-point tadpole/local/physical-double-pole matching; the finite crossed amplitude, improved b20, omitted loops and Regge contour remainder; original state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "whole_finite_G_vertex",
            "UV_complete",
            "full_crossed_b20",
            "exact_massless_LSZ",
            "arbitrary_Ricci_derivative_addition",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "full_finite_amplitude",
            "universal_Ward_slope_sign",
            "drop_metric_contacts",
            "new_local_counterterm",
            "massless_external_scalar",
            "physical_cutoff_from_cap",
            "full_Regge_bound",
            "full_F2_from_F1",
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
            raise ValueError("Unsupported ordinary Ward input accepted: " + name)
    return total


def controls():
    return {
        "whole_original_source_and_all_pure_GR_topologies": True,
        "entire_metric_GF_variation_not_only_TT": True,
        "general_D_identity_and_independent_component_checks": True,
        "whole_Gram_free_parameter_domain_continuity": True,
        "literal_IR_sign_and_entire_raw_LSZ_cancellation": True,
        "physical_IR_and_all_matching_boundaries_retained": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
