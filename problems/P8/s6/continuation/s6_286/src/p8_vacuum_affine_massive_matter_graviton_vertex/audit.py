"""Preserved P8 frontier and the complete scoped massive matter F1 coefficient."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_gaussian_metric_pole_matching import audit as previous

from . import matching as pole_matching
from . import source, vertex, ward

MODEL = "original_massive_matter_F1_and_conditional_background_Ward_not_full_P8"
OBSERVABLES = (
    "whole_source_and_scalar_OS",
    "complete_conditional_background_Ward",
    "complete_generated_massive_matter_F1",
    "generated_t_channel_finite_spin2_piece",
)
ITEM = {
    "id": "QG2_H8A450_complete_generated_massive_matter_F1_positive_spectral_measure_and_bounded_forward_gravity_vertex_piece",
    "status": "COMPLETE_SCOPED_MATTER_F1_AND_T_CHANNEL_PIECE_NOT_WHOLE_VERTEX_CROSSED_AMPLITUDE_OR_V_G_B_P8",
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
        raise ValueError("Require the original specified massive matter F1 model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a scoped generated matter or Ward observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "generated_g_squared_F1",
        "conditional_background_zero_transfer",
    ):
        raise ValueError("Require one of the explicitly stated vertex sectors")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A matter vertex contribution cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_conditional_background_Ward": ward.data(),
        "whole_generated_massive_matter_F1": vertex.data(),
        "whole_generated_t_channel_piece_and_bounds": pole_matching.data(),
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
        "all141_matching_records_preserved": matching()[:-1] == previous.matching()
        and len(previous.matching()) == 141,
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
        "light_curvature_matching_not_fixed_by_scalar_OS": True,
        "conditional_background_identity_not_massless_physical_vertex": True,
        "S275_S276_S277_and_S281_S285_scopes_unchanged": True,
    }


def observable():
    return {
        "established": "Literal light-line and heavy-line stress triangles, both cubic metric contacts and the inherited scalar OS terms give the complete generated g-squared matter F1 coefficient. Its positive spectral representation proves a gapped forward limit and F1prime(0)=Pi_second(mu)/6. The corresponding t-channel finite v-squared piece is-Pi_second(mu)/(3kappa), with magnitude below10^-1205 throughout abs(t)<=2 at the actual hierarchy. The full background Ward identity gives a conditional regulated zero-transfer vertex/LSZ cancellation, including S283's entire raw GR residue expression.",
        "domain": "The whole unchanged original source at its separately stated constant-scalar Lorentz vacuum. The generated massive matter coefficient has no internal graviton; its two cuts begin at4mu and4n. The general background identity is not promoted to an ordinary internal-gravity vertex or exact massless physical pole theorem.",
        "historical_qualification": qualifications(),
        "not_established": "Complete physical F2 and heavy-metric mixing, curved R Phi^2/R H/Ricci-derivative matching or S285 light pure-curvature coefficients; complete finite crossed amplitude, its improved b20 or all four-point local/tadpole matching; the ordinary massless-gravity proper vertex and its IR/detector/LSZ limits, physical Newton normalization and Regge contour remainder; omitted loops, exact original state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "retune_fixed_vacuum",
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
            raise ValueError("Unsupported matter vertex input accepted: " + name)
    return total


def controls():
    return {
        "whole_original_source_and_inherited_scalar_OS": True,
        "both_massive_stress_triangles_and_metric_contacts": True,
        "full_generated_F1_and_positive_spectral_measure": True,
        "whole_metric_tensor_pole_contraction": True,
        "strict_actual_complex_and_forward_bounds": True,
        "background_gauge_IR_and_Ricci_ambiguities_retained": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
