"""Read-only quantitative compact classical response and conditional matching certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as original
from p8_proca_nearby_retarded import verify as parent

from . import audit, band, bumps, matching, probes

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"quantitative-compact-classical-response.json"
PARENT_SHA="1167ad429103423b492f1a53cfd6ee8b8573bb6d2351e2ef8e6e65d6e1ab2f65"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_compact_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen nearby classical retarded relational report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_91_fully_rebuilt":PARENT_SHA,
            "same_actual_classical_background_source_and_full_retarded_response":True,
            "no_temporal_EFT_cutoff_or_quantum_parent_tail_assumed":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A quantitative compact classical response gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.92.QUANTITATIVE_COMPACT_RESPONSE","date":"2026-09-09",
            "status":"EXPLICIT_COMPACT_CLASSICAL_PROBES_HAVE_POSITIVE_FULL_AND_FINITE_SPATIAL_BAND_RETARDED_PAIRINGS_WITH_COMPLETE_REMAINDER_CONTROL; ACTUAL_PARENT_MATCHING_TAIL_TEMPORAL_EFT_QUANTUM_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/probes.md","notes/signal.md","notes/band.md","notes/matching.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":len(identities),"proof_checks":gates,
            "literal_compact_bumps_normalizers_and_boundary_recurrence":serialize(bumps.data()),
            "explicit_compact_probe_geometry_norms_and_actual_widths":
                serialize({"formal_constants":probes.formal(),"actual_remainder_and_widths":probes.actual()}),
            "finite_spatial_frequency_contribution_and_complete_tail_budget":serialize(band.data()),
            "conditional_common_parent_matching_and_tail_requirements":serialize(matching.data()),
            "controls":audit.controls(),
            "verdict":"On the same actual nearby classical bounce, explicit nonnegative smooth compact source and detector probes have every support pair matter-spacelike and an exact real retarded pairing at least 3D/4, with D a declared positive rational. Both source normalizers, all four detector bumps, the actual moving clock-front cap and the complete spatial L2 remainder are retained. A finite symbolic comoving spatial cutoff leaves a contribution at least D/2. The widths retain an extremely conservative full low-frequency exponential. A matter-causal common-parent comparison would require separate low-band matching and parent high-band tail errors with total below D/2; neither is established. There is no temporal-frequency, interacting EFT, quantum or UV-completion transfer.",
            "not_established":["Numerically synthesized or physically resolved probes, a controlled nonlinear source apparatus or detector noise bound",
                "Joint temporal/spatial EFT validity, optimized carrier/width, interacting cutoff, omitted operators, loops or heavy-threshold errors",
                "Actual parent low-band matching or actual parent high-band tail bounds",
                "A quantum state, fixed profile, response or inverse on the nearby background",
                "Finite Wilson matching, healthy full-candidate vacuum, common-parent V/G/B, global bounce, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact native bump recurrence/normalization anchors, whole-support geometry and positive front-cap lower, full actual symbolic remainder and width, Plancherel operator tail factors, source/detector norm powers, strict rational conditional error budgets and negative cancellation control, with written all-orders smoothness and paired-response proofs. Every parent is natively rebuilt. No floating estimate, library patch, unspecified localization limit, discarded low frequency or assumed EFT/UV validity. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The quantitative compact classical response report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.92.QUANTITATIVE_COMPACT_RESPONSE replay passed; actual parent matching, tails, temporal EFT, quantum and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
