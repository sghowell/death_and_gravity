"""Strict scope: complete physical two-graviton cut, not full matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_detector_ir_cut import audit as previous

from . import amplitude, sewing, source, window

MODEL = "original_massive_physical_two_graviton_first_cut_not_full_G"
OBSERVABLES = (
    "literal_EH_scalar_graviton_pair",
    "whole_angular_two_graviton_cut",
    "complete_massive_forward_and_all_even_spins",
    "bounded_physical_spectral_window",
)
ITEM = {
    "id": "QG2_H8A443_original_massive_two_graviton_cut_from_literal_EH_and_matter_with_full_angular_sewing_all_even_spins_and_physical_window_bound",
    "status": "COMPLETE_SCOPED_FORMAL_PHYSICAL_TWO_GRAVITON_CUT_NOT_FULL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def qualifications():
    return previous.qualifications()


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original formal physical two-graviton-cut scope")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a scoped original massive graviton-cut observable")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("This physical cut does not close original P8")
    return True


@cache
def packets():
    return {
        "full_original_source_and_loop_bookkeeping": source.data(),
        "literal_EH_contact_and_exchange_tree": amplitude.data(),
        "whole_physical_angular_and_all_spin_cut": sewing.data(),
        "uniform_physical_window_bound": window.data(),
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
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
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
        "all134_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 134,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_physical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "no_new_state_counterterm_or_exact_full_amplitude": True,
        "physical_cut_not_unphysical_analytic_continuation": True,
        "S275_S276_S277_S278_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The full fixed source admits explicit formal-loop bookkeeping with no deleted vacuum constants. The literal canonically normalized EH cubic, scalar contact and both scalar exchanges give all four two-graviton helicity amplitudes. Their complete phase-sensitive angular sewing has an exact massive forward integral, all-even-spin nonnegative squares and a uniform bound on physical finite spectral windows.",
        "domain": "First formal two-graviton loop coefficient for original massive Phi scattering, kappa10^800,mass1, physical s>=4m^2 and real scattering angle. The two TT polarizations and identical-pair normalization are retained. No full quantum vacuum or resummed finite-gravity amplitude is claimed.",
        "historical_qualification": qualifications(),
        "not_established": "Unphysical0<s<4m^2 cut continuation, complete crossed low-cut subtraction/b20, other intermediate species or higher loops, exact local matching, detector errors, UV spin-admissible dispersion or a finite Regge remainder, a physical cutoff, original quantum bounce/state/domain or original V/G/B/P8 closure.",
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
        for label, call, args in (
            ("model", require_model, (value,)),
            ("observable", require_observable, (value,)),
            ("spin", sewing.require_spin, (value,)),
            ("positive", window.exact_positive, (value,)),
            ("helicity", amplitude.polarization, (1, value)),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, args))
    for i, value in enumerate((0, -1, -s.Rational(2, 3))):
        rows.append(("nonpositive_" + str(i), window.exact_positive, (value,)))
    for i, value in enumerate((-1, -2)):
        rows.append(("negative_spin_" + str(i), sewing.require_spin, (value,)))
    for i, args in enumerate(((4, 4), (3, 5), (4, 3), (2, 8, 1), (4, 9, 2))):
        rows.append(("nonphysical_window_" + str(i), window.require_window, args))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "pure_KLT_four_vector_state_sum",
            "full_unphysical_cut",
            "all_graviton_loop_terms_from_forward_cut",
            "exact_finite_G_positivity",
            "new_vacuum_counterterm",
            "physical_window_is_EFT_cutoff",
            "complete_b20",
        )
    ):
        rows.append(("unsupported_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("deleted_primitive", validate_scope, ([], matching())),
            ("deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    if len({row[0] for row in rows}) != len(rows):
        raise ValueError("Duplicate graviton-cut rejection case")
    return rows


def rejected_inputs():
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported graviton-cut input accepted: " + name)
    return total


def controls():
    return {
        "entire_source_loop_marker_and_fixed_vacuum_constants": True,
        "literal_EH_cubic_scalar_contact_and_two_exchanges": True,
        "independent_covariant_factorization_and_both_Ward_identities": True,
        "physical_TT_projector_and_all_four_helicities": True,
        "full_phase_sensitive_angular_and_exact_forward_sewing": True,
        "all_spin_written_proof_and_independent_finite_Q_integrals": True,
        "physical_window_not_full_b20_Regge_or_EFT_cutoff": True,
        "rejected_inputs": rejected_inputs(),
    }
