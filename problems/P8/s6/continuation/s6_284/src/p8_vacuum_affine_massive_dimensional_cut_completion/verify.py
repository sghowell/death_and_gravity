"""Read-only whole dimensional cut and Gram-principal-part report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_common_gravity_masters import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-massive-dimensional-cut-completion.json"
)
PARENT_SHA = "affd15ac3d06770b6c46b601138212bc1ecf537e91a636447141702568b72ab3"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob("src/p8_vacuum_affine_massive_dimensional_cut_completion/*.py")
        )
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen original common cut-master and external-leg report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_283_common_cut_masters_external_legs_and_entire_ancestry_rebuilt": PARENT_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A dimensional cut completion proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.284.ORIGINAL_FULL_DIMENSIONAL_TT_SEWING_COMMON_CUT_MASTER_COEFFICIENTS_AND_FINITE_BUBBLE_RATIONAL_GRAM_PRINCIPAL_PART_COMPLETION",
        "date": "2026-09-15",
        "status": "SCOPED_DIMENSIONAL_CUT_COEFFICIENTS_AND_GRAM_PRINCIPAL_PARTS; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/threshold.md",
            "notes/polarizations.md",
            "notes/dimensional.md",
            "notes/completion.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_threshold_jets": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_dimensional_polarizations_and_finite_rational_terms": serialize(
            {
                "packets": {name: payload(packets[name]) for name in names[2:]},
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete original Einstein/matter TT tree is reconstructed in symbolic D using the written spectator-trace degree bound. The physical symmetric-traceless projector and whole D-dimensional angular reduction give both species' common scalar-master coefficients. Their ultraviolet bubble poles fix a nonzero finite evanescent rational contribution. Exact massive triangle threshold jets and crossing-symmetric rational choices remove the artificial integer Gram principal parts, including their forced tadpole-dependent part. The complete massive tadpole coefficient, threshold distributions, physical finite pole/local/vacuum matching, other sectors and full finite amplitude remain undetermined. This is not a new exact quantum vacuum, a finite detector/Regge theorem or the original bounce. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Complete massive tadpole coefficient or threshold-supported distributional matching",
            "Physical finite pole, local and vacuum normalization, full finite crossed amplitude or b20",
            "Other heavy/contact/matter local sectors, omitted loops or finite detector/Regge control",
            "A physical EFT cutoff, unique full off-background dimensional completion or exact LSZ vacuum",
            "Original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole original source, literal Einstein and matter vertices, exact symbolic-D invariant reconstruction, physical TT trace contractions, complete D-dimensional sphere moments and both-species master coefficients. Written degree, threshold and analytic-continuation arguments are not FORMALIZED. Private numerical diagnostics are not continuum certificates. Native/direct/ordinary/CLI retain original SymPy; the exact-GCD adapter is restricted to a captured complete FULL regression. Fixed sources, raw reports, vacuum constants and all historical physical qualifications remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The dimensional cut and Gram-principal-part report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.284 dimensional cut and Gram-principal-part replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
