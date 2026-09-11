"""Read-only two-loop full derivative-action physical-source transport."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_full_heavy_source import verify as parent

from . import audit, calibration, topology, transport, vertices, ward

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-physical-source-map.json"
PARENT_SHA = "7bb4498c44c3a56b8e5b7b6a8c8c608193708c94bae80ce1b09e91db4f58a9a6"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_physical_source_map/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen physical-source transport parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_159_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "full_transformed_action_and_matched_source_not_bare_truncation_or_ordinary_Psi_MS": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A physical-source transport proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.160.COMPLETE_TWO_LOOP_PHYSICAL_SOURCE_DERIVATIVE_MAP_TRANSPORT",
        "date": "2026-09-10",
        "status": "FIXED_ORDER_CANONICAL_PHI_POLE_B2_VACUUM_H_SOURCE_TRANSPORTED_WITH_FULL_DERIVATIVE_ACTION_AND_MATCHED_PHYSICAL_SOURCE; NOT_ORDINARY_PSI_MS_GLOBAL_INVERSE_PHYSICAL_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/prescription.md",
            "notes/topology.md",
            "notes/vertices.md",
            "notes/ward.md",
            "notes/transport.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_generated_action_source_and_loop_support": serialize(
            {"topology": payload(topology.data()), "vertices": payload(vertices.data())}
        ),
        "regulated_Ward_and_source_transport": serialize(
            {"ward": payload(ward.data()), "transport": payload(transport.data())}
        ),
        "actual_physical_source_two_loop_enclosures": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "With the literal fixed cubic derivative F and the full transformed parent action, counterterms, Jacobian and physical source J F(Psi), the named two-loop canonical Phi pole, b2, vacuum and stationary-H references transfer unchanged. Four-source two-loop support requires map degree three; generated octic scalar, cubic-Phi Yukawa and corresponding counterterms are retained. General Gaussian-monomial Ward identities hold through cubic map order and omission controls detect the missing terms. The continuum statement uses the regulated formal functional identity with local Jacobian loops scaleless in dimensional regularization. This does not compute an independently minimally subtracted ordinary-Psi composite operator, its second residue, a global inverse, or the physical truncation error. V/G/B and original P8 remain open.",
        "not_established": [
            "Ordinary-Psi two-loop MS composite mixing, residue positivity or a global off-shell coordinate norm",
            "Equivalence of the bare S6.109 quartic derivative truncation without the generated action and matched sources",
            "External fermion/gauge coordinate dictionaries, rolling-state transport or finite-gravity map domains",
            "Physical finite-EFT omitted-order errors, V contours/cuts, G, B or original P8 closure",
        ],
        "verification_boundary": "Exact literal jet homogeneity, weighted halfedge/Euler support, full two-site scalar/Yukawa/counterterm substitution, arbitrary-monomial Ward recurrences through map degree three and inherited rational parent enclosures. Independent tests use directly integrated nonlinear Gaussian pullbacks, a coupled two-variable map, omitted-interaction controls and independent formal inverse/LSZ expansions. Finite-dimensional diagnostics do not replace the written continuum argument. Written analysis and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science use unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The physical-source transport report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.160.COMPLETE_TWO_LOOP_PHYSICAL_SOURCE_DERIVATIVE_MAP_TRANSPORT replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
