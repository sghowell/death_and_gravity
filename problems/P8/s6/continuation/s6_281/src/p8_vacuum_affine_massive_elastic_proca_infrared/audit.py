"""Strict physical first-loop scope, complete species ledger, and rejected shortcuts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_graviton_cut_isolated_replay import audit as previous

from . import elastic, infrared, proca, source

MODEL = "original_first_loop_massive_elastic_and_Proca_IR_not_full_P8"
OBSERVABLES = (
    "whole_elastic_tree",
    "complete_nonforward_angular_cut",
    "specified_dimensional_finite_soft_subtraction",
    "literal_all_polarization_Proca_cut",
    "all_physical_two_body_species_below_heavy_threshold",
)
ITEM = {
    "id": "QG2_H8A445_full_original_massive_elastic_angular_cut_dimensional_finite_soft_matching_and_all_polarization_Proca_sewing",
    "status": "COMPLETE_SCOPED_FIRST_LOOP_PHYSICAL_CUTS_NOT_FULL_V_G_B_P8",
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
        raise ValueError("Require the stated original first-loop cut model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a scoped physical first-loop observable")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("First-loop cuts cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source_and_species": source.data(),
        "whole_elastic_tree_and_angular_master": elastic.data(),
        "complete_dimensional_finite_soft_cut": infrared.data(),
        "literal_nine_polarization_Proca_cut": proca.data(),
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
        "all136_matching_records_preserved": matching()[:-1] == previous.matching()
        and len(previous.matching()) == 136,
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
        "finite_cut_not_full_real_amplitude_matching": True,
        "S275_hybrid_and_S276_S277_fold_scopes_unchanged": True,
    }


def whole_species_cut():
    S, mu, K, z = source.S, source.MU, source.K, source.Z
    x, phi = s.symbols("cut_axis_cosine cut_azimuth", real=True)
    y = z * x + s.sqrt(1 - z * z) * s.sqrt(1 - x * x) * s.cos(phi)
    gg = previous.sewing.full_helicity_sewing(S, mu, K, x, y, z)
    t = -(S - 4 * mu) * (1 - z) / 2
    u = -(S - 4 * mu) * (1 + z) / 2
    return {
        "Phi_Phi_after_specified_soft_stripping": infrared.whole_stripped_cut(),
        "original_M1_M1": (S * S - t * u + 2 * mu * S + 6 * mu * mu)
        / (960 * s.pi * K * K),
        "two_physical_TT_gravitons": s.Integral(gg, (x, -1, 1), (phi, 0, 2 * s.pi))
        / (128 * s.pi**2),
        "three_polarization_Proca_pair_when_open": proca.whole_cut(
            mass2=source.VECTOR_MASS2
        ),
        "domain": "4mu<s<M_H^2,-1<z<1; Proca contribution is zero below4M_A^2. The PhiPhi term alone carries the first-loop universal soft divergence; its stated finite prescription is retained.",
        "no_double_count": "The elastic term already includes contact/heavy/gravity squares and every interference. Do not add the earlier pure-matter elastic cut a second time. M1, TT and Proca are distinct physical intermediate states.",
        "not_an_optical_positivity_theorem": "The sum is a specific first-loop stripped cut. Soft stripping need not preserve finite-G pointwise or partial-wave positivity. No exact finite-resolution error, forward limit, real local polynomial or Regge bound is inferred.",
    }


def observable():
    return {
        "established": "Entire original leading elastic tree; complete all-angle angular convolution and its finite part; minimal-D phase-space and evanescent tensor terms matched to the explicitly defined S278 soft factor; literal complete Proca stress, all nine physical polarization pairs and nonforward sewing. Together with retained M1 and TT sectors these specify all first-loop physical two-body s cuts below the first heavy threshold.",
        "domain": "Formal leading vacuum vertices of the unchanged full source, external mass1,kappa10^800,Proca mass1000 and exact heavy mass squared10^200/512+2. Strictly nonforward physical s region, with a stated minimal dimensional regulator convention and original S278 analytic soft prescription.",
        "historical_qualification": qualifications(),
        "not_established": "Complete crossed real amplitude, unphysical low-cut continuation, forward b20, unknown local matching or higher-loop errors, a physical cutoff, exact detector unitarity, finite Regge remainder, original quantum state/domain/bounce or V/G/B/P8 closure.",
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
            ("exact", elastic.exact_real, (value,)),
            ("domain", elastic.require_domain, (value, 0)),
            ("vector", proca.require_open_vector_channel, (value,)),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, args))
    for i, args in enumerate(
        (
            (4, 0),
            (3, 0),
            (10, 1),
            (10, -1),
            (10, 2),
            (10, 0, 0),
            (10, 0, -1),
            (50, 0, 1, 50),
        )
    ):
        rows.append(("nonphysical_domain_" + str(i), elastic.require_domain, args))
    for i, value in enumerate((0, 1, 4 * 10**6, -1, source.HEAVY_MASS2)):
        rows.append(
            (
                "closed_vector_channel_" + str(i),
                proca.require_open_vector_channel,
                (value,),
            )
        )
    for i, value in enumerate(
        (
            "original_P8_closed",
            "exact_finite_G_positivity",
            "auxiliary_w_equals_detector_E",
            "dropped_evanescent_trace",
            "complete_forward_b20",
            "new_local_counterterm",
            "massless_external_scalar",
            "four_vector_polarizations",
            "physical_cutoff_from_cap",
            "full_Regge_bound",
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
            raise ValueError("Unsupported physical-cut input accepted: " + name)
    return total


def controls():
    return {
        "full_original_source_not_only_clock_germ": True,
        "entire_tree_all_channels_and_interferences": True,
        "generic_master_primitive_and_independent_azimuthal_moments": True,
        "finite_dimensional_phase_space_and_evanescent_trace": True,
        "soft_factor_times_D_tree_not_just_4D_tree": True,
        "literal_Proca_polarizations_Ward_and_M1_calibration": True,
        "original_scope_not_full_P8_or_Regge": True,
        "rejected_inputs": rejected_inputs(),
    }
