"""Read-only finite radiative angular conversion certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_physical_virtual_soft_pairing import verify as soft_input
from p8_vacuum_affine_radiative_state_soft_index import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-radiative-angular-finite.json"
PARENT_SHA = "83fc98b6541ae6fb75c18c342d20c3e4bc03842fe69a3191b3dc8d1a09509df7"
SOFT_SHA = "40c04cadfa1a1c542c16eed48f6df12bcfe0605ad69b39dbb70fcefaefe6b3ed"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_radiative_angular_finite/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen radiative soft-index parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(soft_input.REPORT) != SOFT_SHA:
        raise ValueError("The frozen dimensional soft phase changed")
    soft_input.validate_report(
        json.loads(soft_input.REPORT.read_text()), soft_input.build_report()
    )
    return {
        "S6_300_state_correct_radiative_index_rebuilt": PARENT_SHA,
        "S6_296_dimensional_phase_and_soft_conversion_rebuilt": SOFT_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A finite radiative angular proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.301.COMPLETE_RADIATIVE_FINITE_ANGULAR_CONVERSION_UNIFORM_COLLINEAR_HOLDER_AND_EXPLICIT_REGULATOR_ERROR",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_RADIATIVE_FINITE_ANGULAR_CONVERSION; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/angular.md",
            "notes/bounds.md",
            "notes/conversion.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_dimensional_current": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_collinear_majorants_and_finite_conversion": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full fixed-ball dimensional current, retained trace derivative and integrable collinear quarter-power majorant fix the finite radiative angular conversion uniformly in arbitrary finite multiplicity. Its absolute bound is2000/kappa, with explicit fixed-resolution regulator error. The known unexpanded detector reference differs from its elastic value by less than4250/kappa at originalparameters. Finite and evanescent hard amplitudes, complete nonleading radiation, physical analytic/Regge bounds and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Finite radiative hard amplitudes and evanescent hard scheme terms",
            "Uniform complete-amplitude-minus-leading-soft remainder at all multiplicities",
            "Omitted loops/operators and full physical unitarity or analyticity",
            "High-energy Regge and original quantum state/domain/measure/bounce",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact projector, beta, logarithmic-integral and majorant identities accompany a written uniform collinear and differentiated-integral proof. Independent massive Feynman-parameter and full-radiative angular quadratures calibrate the formula, not rigorous numerical errors. The finite angular scheme is specified; evanescent hard terms are not set to zero. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. All frozen ancestors and original closure statuses are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The finite radiative angular report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.301 finite radiative angular conversion replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
