"""Read-only overlapping-time actual-state memory/contact remainder."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_separated_response import verify as parent

from . import audit, contact, memory, response, tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-retarded-state-remainder.json"
)
PARENT_SHA = "6fec9e8992d0b6aeb13a3c95cb32e2d32c8d2478d03391ae6c94fe6f62810a6d"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_retarded_state_remainder/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_196_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A retarded actual-state remainder gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.197.ACTUAL_SPATIAL_RETARDED_CURRENT_STATE_REMAINDER_WITH_OVERLAPPING_SUPPORTS_COMPLETE_CONTACT_AND_REGULATOR_TAIL",
        "date": "2026-09-11",
        "status": "ACTUAL_SPATIAL_RETARDED_CURRENT_STATE_REMAINDER_WITH_OVERLAPPING_SUPPORTS_COMPLETE_CONTACT_AND_REGULATOR_TAIL; NOT_SINGULAR_REFERENCE_COVARIANT_MATCHING_FULL_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tail.md",
            "notes/memory.md",
            "notes/response.md",
            "notes/contact.md",
            "notes/limit.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_retarded_memory_and_contact_remainder": serialize(
            {"memory": payload(memory.data()), "contact": payload(contact.data())}
        ),
        "all_momentum_tail_and_reference_boundary": serialize(
            {"tail": payload(tail.data()), "response": payload(response.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged actual all-order Gaussian Proca state has an absolutely convergent full spatial first-order current comparison on overlapping smooth compact supports. The subtraction is the same finite-regulator mode formula with constant-alpha W8 readouts, not a new physical state or renormalization prescription. Complete pointwise pair-product differences retain both mixed terms and the error square while keeping the retarded step. The full second-metric contact retains all three physical modes, its noncommuting local product and its distinct one-momentum regulator tail. Exact all-momentum majorants give a complete comparison bound2e23 in spatial H1 and time L2 test norms, and regulator error1e27/K uniformly over all external momenta. Two canonical tensor factors give8e-777 and4e-773/K. The reference has a possible CCR defect and a nonzero equation residual; its singular memory/contact still needs matching to the original covariant retarded prescription on the diagonal. No full response inverse, finite-amplitude interacting background, physical cutoff or original V/G/B and P8 closure follows.",
        "not_established": [
            "Original covariant singular reference retarded matching on overlapping supports and the diagonal",
            "A physical reference CCR state, reference effective-action derivative or independent reference Ward identity",
            "Finite-amplitude response, full feedback inverse, interacting background or stability and physical cutoff matching",
            "Full parent state/measure/loops, vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact complete pair-product and contact differences, all-momentum radial majorants, low/high external split and canonical factors. Independent angular/radial integrals, literal overlapping-time Fock operators, noncommuting ten-field contacts and CCR-defect controls supplement continuous written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The analytic convergence and test-space arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The retarded actual-state remainder report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.197.ACTUAL_SPATIAL_RETARDED_CURRENT_STATE_REMAINDER_WITH_OVERLAPPING_SUPPORTS_COMPLETE_CONTACT_AND_REGULATOR_TAIL replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
