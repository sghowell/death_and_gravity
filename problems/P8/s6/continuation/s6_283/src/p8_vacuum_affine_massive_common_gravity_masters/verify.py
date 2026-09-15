"""Read-only common gravitational cut-master and external-leg report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_graviton_threshold_sheet import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-massive-common-gravity-masters.json"
)
PARENT_SHA = "198c5f346a08b7bcba3a36a297efaf53b9a95bfe62019892512522689549af2b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_massive_common_gravity_masters/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen original graviton normal-cut sheet changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_282_complete_graviton_normal_cut_sheet_and_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A common gravitational master proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.283.ORIGINAL_COMMON_CROSSED_GRAVITY_MASTER_CUT_BASIS_AND_LITERAL_SCALAR_EXTERNAL_LEG_SOFT_MATCHING",
        "date": "2026-09-15",
        "status": "SCOPED_COMMON_GRAVITY_CUT_MASTERS_AND_LSZ_IR; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/masters.md",
            "notes/graviton.md",
            "notes/elastic.md",
            "notes/legs.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_scalar_masters": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_common_cut_basis_and_external_legs": serialize(
            {
                "packets": {name: payload(packets[name]) for name in names[2:]},
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same six ordered alternating scalar boxes reproduce the complete original two-graviton and pure-gravity elastic cuts. Exact Feynman parameter and dimensional measure calculations retain the complete finite box/triangle terms; evanescent stress traces reproduce every finite pure-gravity elastic contribution. The literal harmonic-gauge scalar self-energy has zero on-shell mass shift in D in this pure-gravity sector, while its residue has separate UV and IR contributions. Four LSZ legs restore the self term and the full crossed result matches the accepted massive analytic soft pole. A massless whole-bubble calibration is also exact. This is a cut-matched master representation, not the full finite rational amplitude, physical pole/local matching, finite detector/Regge bound, exact quantum vacuum or original bounce. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Exact interacting vacuum or LSZ existence, total mass shift or finite physical pole normalization",
            "Cut-free rational and local matching, full finite crossed amplitude, transfer subtraction or b20",
            "Other heavy/contact/matter sectors, omitted loops, finite detector errors or finite Regge allowance",
            "A physical EFT cutoff, unique off-background dimensional completion or full UV theory",
            "Original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole original source, literal off-shell stress and Ward identities, exact parameter primitives and gamma recurrences, complete both-species cut identities and finite dimensional trace matching. Written uniform integration and analytic continuation arguments are not FORMALIZED. The common basis retains all named rational-completion and renormalization boundaries. Native/direct/ordinary/CLI retain original SymPy; the exact-GCD adapter is restricted to a captured complete FULL regression. No numerical diagnostic is a continuum proof, and all frozen historical scopes remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The common gravitational master and external-leg report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.283 common gravitational master and external-leg replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
