"""Read-only anchored full-function classical decoupling certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as parent

from . import amplitude, audit, family, gravity, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-canonical-affine-decoupling.json"
PARENT_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_canonical_affine_decoupling/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual affine parent/state source changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_176_fully_rebuilt": PARENT_SHA,
        "same_full_base_action_and_all_prior_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual classical decoupling proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.177.FULL_FIXED_CANONICAL_FUNCTION_AFFINE_FAMILY_AND_ACTUAL_CLASSICAL_GRAVITATIONAL_DECOUPLING",
        "date": "2026-09-11",
        "status": "ANCHORED_FULL_FUNCTION_CLASSICAL_DECOUPLING_AND_TREE_AMPLITUDE; NOT_QUANTUM_CONTOUR_TRUNCATION_PHYSICAL_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/family.md",
            "notes/affine.md",
            "notes/gravity.md",
            "notes/source.md",
            "notes/amplitude.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "fixed_canonical_family_and_gravitational_limit": serialize(
            {
                "family": payload(family.data()),
                "germs": payload(family.germs()),
                "gravity": payload(gravity.data()),
            }
        ),
        "full_source_and_tree_observables": serialize(
            {
                "source": payload(source.data()),
                "coefficients": payload(source.coefficients()),
                "amplitude": payload(amplitude.data()),
                "pole": payload(gravity.pole()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A continuous family through the complete S6.174/S6.176 base action holds the full independent canonical scalar functions f,r,a3 and all vacuum masses fixed. Its regular affine chart and quotient extend on a fixed canonical domain; only explicitly gravity-dependent Ia coefficients vary. After retaining the exact Einstein boundary, the full classical action tends uniformly on compact canonical jets to the complete fixed scalar action plus free normalized tensor, M1 and massive-Proca spectators. The complete canonical vector source and its contact vanish with explicit coefficient/norm rates. The exact four-scalar tree amplitude retains nonzero gamma and has b2=4lambda. Fixed-transfer graviton decoupling is not uniform in the forward limit. No quantum-continuum, full contour/cut/truncation, common interacting parent, physical cutoff, G/B or original P8 closure follows.",
        "not_established": [
            "A regulator-removed quantum decoupling limit, all-order amplitude or uniform physical omitted-order bound",
            "A full vacuum contour/cut relation and a UV completion of the complete scalar action",
            "A finite-gravity IR/Regge prescription or quantitative finite remainder",
            "A common interacting affine/light/metric state, physical cutoff, controlled quantum background/response or off-base clock solution",
            "V/G/B or original P8 closure",
        ],
        "verification_boundary": "Literal full analytic function derivatives and anchored chain rules, complete regular source scaling, exact Einstein quadratic normalization and retained boundary, labelled tensor four-point vertices and dense stress/projector controls supplement written compact-domain and finite-jet proofs. Exact algebra and written analysis are not FORMALIZED or a quantum UV theorem. Native, direct science, ordinary and CLI use original SymPy; full regression alone uses the separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual fixed-canonical-function decoupling report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.177.FULL_FIXED_CANONICAL_FUNCTION_AFFINE_FAMILY_AND_ACTUAL_CLASSICAL_GRAVITATIONAL_DECOUPLING replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
