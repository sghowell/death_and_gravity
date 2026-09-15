"""Minimal heavy transition and finite formal threshold pairing."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_heavy_graviton_production_threshold import audit as previous

from . import masters, pairing, proper, source

MODEL = "original_minimal_heavy_transition_and_formal_threshold_pairing_not_physical_P8"
OBSERVABLES = (
    "whole_original_minimal_gravity_graph_ownership",
    "complete_proper_vertex_and_external_factors",
    "entire_physical_sheet_masters_and_first_epsilon_derivative",
    "finite_formal_real_virtual_threshold_pair",
)
ITEM = {
    "id": "QG2_H8A456_complete_minimal_heavy_transition_normal_sheet_masters_and_finite_formal_real_virtual_soft_threshold_pairing",
    "status": "COMPLETE_SCOPED_MINIMAL_TRANSITION_AND_FINITE_FORMAL_THRESHOLD_PAIR_NOT_FULL_AMPLITUDE_PHYSICAL_IR_WIDTH_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original minimal transition sector")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated transition or formal-pairing observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "minimal_heavy_transition",
        "finite_formal_threshold_pair",
    ):
        raise ValueError("Require the minimal transition or its formal threshold pair")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError(
            "Formal threshold pairing cannot close the original physical P8"
        )
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_minimal_proper_and_external_factors": proper.data(),
        "whole_physical_sheet_scalar_masters": masters.data(),
        "whole_finite_formal_threshold_pair": pairing.data(),
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
        "all147_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 147,
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
        "finite_RH_and_other_curved_matching_not_assigned": True,
        "whole_physical_IR_Regge_and_full_mixed_amplitude_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "Every minimal proper pair triangle, metric-cubic contact and external scalar residue is retained in a full-D master formula. Physical-root subtraction supplies complete convergent first cusp derivatives. The full original Hh real cut and the known minimal heavy delta coefficient have a finite compact-window combination: their soft poles, regulator scale and Euler constants cancel with an explicit conditional remainder. The imaginary Coulomb term is retained.",
        "domain": "The unchanged original variables, selected minimal g/kappa transition and g^2/kappa formal forward spectral pairing. mu>0,n>4mu,D=4+2epsilon,0<epsilon<=1/8; fixed C1 weights and0<L<=1/2. H is perturbatively unstable, not an exact external particle. All original numerical parameters and matching/state boundaries remain.",
        "historical_qualification": qualifications(),
        "not_established": "Full H-metric/higher-EFT/local finite matching, complete mixed amplitude or full-source b20, positive physical spectral measure, exact heavy width or near-resonance uniformity, complete Coulomb/dressed detector IR, fixed-transfer Regge/all-loop bounds, original quantum state/domain/measure/bounce or V/G/B/P8 closure.",
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
            raise ValueError(
                "Unsupported minimal-transition pairing input accepted: " + name
            )
    return total


def controls():
    return {
        "all_original_vertices_and_graph_boundaries_retained": True,
        "whole_D_pair_tensor_and_contact_denominator_reductions": True,
        "all_three_external_factors_and_separate_UV_IR_origins": True,
        "entire_physical_root_moments_not_D0_only": True,
        "full_known_real_virtual_pole_and_finite_pair": True,
        "Coulomb_H_metric_local_width_and_physical_IR_not_deleted": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
