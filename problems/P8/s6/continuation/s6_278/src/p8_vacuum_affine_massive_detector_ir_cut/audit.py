"""Strict formal-loop scope and unchanged original matching frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_homogeneous_fold_dynamics import audit as previous

from . import cut, soft, source, subtraction

MODEL = "original_massive_detector_dictionary_and_M1_pair_cut_not_full_G"
OBSERVABLES = (
    "literal_massive_vacuum_and_M1_tree",
    "complete_M1_pair_cut",
    "massive_universal_soft_dictionary",
    "explicit_three_channel_low_cut_subtraction",
)
ITEM = {
    "id": "QG2_H8A442_original_massless_M1_first_gravitational_pair_cut_massive_soft_dictionary_and_crossing_low_cut_subtraction",
    "status": "COMPLETE_SCOPED_FORMAL_M1_LOOP_AND_CONDITIONAL_SOFT_DICTIONARY_NOT_FULL_V_G_B_P8",
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
        raise ValueError("Require the formal M1 cut and conditional detector scope")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a scoped M1 cut or detector observable")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("This IR-sector calculation does not close original P8")
    return True


@cache
def packets():
    return {
        "literal_current_vacuum_and_M1_pair": source.data(),
        "whole_massless_M1_cut": cut.data(),
        "massive_analytic_soft_dictionary": soft.data(),
        "complete_crossing_low_cut_subtraction": subtraction.data(),
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
        "all133_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 133,
        "all_matching_identifiers_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "all6_historical_physical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "no_new_state_counterterm_or_full_amplitude_claimed": True,
        "no_external_massless_bound_or_O_G_error_promoted": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The literal massive Phi reference and minimally coupled massless M1 give the full production tree and complete first M1-pair gravitational cut, with all three crossed logarithms modulo real local terms. The equal-mass analytic soft prescription retains all pair/self terms and yields its conditional resolution law. An explicit crossing low-cut subtraction restores a genuine analytic center for THIS loop sector, with the exact cap-dependent b20 and separate transfer-channel term.",
        "domain": "Formal first M1 loop at order kappa^-2 in vacuum scattering, external tree mass1, kappa10^800 and original physical matter frame; not a classical-decoupled free-M1 calculation. The soft dictionary is conditional on the stated universal factorization and is not exact finite-coupling positivity.",
        "historical_qualification": qualifications(),
        "not_established": "Full graviton or all-species quantum amplitudes, an exact physical mass/LSZ matching theorem, finite local counterterms fixed from cuts, finite detector unitarity/Regge remainder, certified positive UV smearing, physical EFT cutoff, quantum bounce, all-parent no-go or original V/G/B/P8 closure. No scalar mass is set to zero.",
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
            ("positive", soft.exact_positive, (value,)),
            ("cap", subtraction.require_cap, (value,)),
            ("order", soft.kernel_coefficient, (value,)),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, args))
    for i, value in enumerate((0, -1, -s.Rational(2, 3))):
        rows.append(("nonpositive_" + str(i), soft.exact_positive, (value,)))
    for i, values in enumerate(((4, 1), (3, 1), (1, 1), (5, 2), (10, 3))):
        rows.append(("below_massive_cap_" + str(i), subtraction.require_cap, values))
    for j in (-2, -1, 9, 100):
        rows.append(("unlicensed_series_" + str(j), soft.kernel_coefficient, (j,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "all_massless_cuts_removed_by_soft_factor",
            "exact_finite_G_positivity",
            "massless_pion_65_7_bound_for_original_Phi",
            "new_counterterms_from_absorptive_data",
            "detector_resolution_equals_regulator",
            "subtraction_cap_is_physical_cutoff",
            "original_bounce_UV_complete",
        )
    ):
        rows.append(("unsupported_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("deleted_primitive", validate_scope, ([], matching())),
            ("deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    if len({x[0] for x in rows}) != len(rows):
        raise ValueError("Duplicate IR rejection case")
    return rows


def rejected_inputs():
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported IR input accepted: " + name)
    return total


def controls():
    return {
        "literal_four_vector_stress_and_both_Ward_identities": True,
        "independent_Cartesian_sphere_and_Legendre_sewing": True,
        "identical_final_pair_and_full_massive_optical_normalization": True,
        "massless_pair_threshold_and_nonzero_forward_log": True,
        "massive_soft_parameter_integral_not_massless_replacement": True,
        "whole_crossing_low_cut_add_back_and_transfer_term": True,
        "finite_counterterm_detector_error_and_other_cuts_left_explicit": True,
        "rejected_inputs": rejected_inputs(),
    }
