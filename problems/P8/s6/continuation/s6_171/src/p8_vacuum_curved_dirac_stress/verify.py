"""Read-only absolute specified curved free Dirac one-loop tensor."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_curved_dirac_state import verify as parent

from . import audit, bloch, curvature, reference, stress, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-curved-dirac-stress.json"
PARENT_SHA = "68b4979bf32690ac2ebcd3ee125f003097d58c27f74cbdbfb78d9c4a958e6cc7"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_curved_dirac_stress/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual curved Hadamard-state report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_170_fully_rebuilt": PARENT_SHA,
        "same_actual_CD_free_in_out_Hadamard_states": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An absolute curved-stress proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.171.ABSOLUTE_CURVED_FREE_DIRAC_ONE_LOOP_TENSOR_WITH_EXPLICIT_EVANESCENT_AND_NEWTON_REFERENCE",
        "date": "2026-09-10",
        "status": "BOUNDED_COMPLETE_SPECIFIED_QUADRATIC_CURVED_ONE_LOOP_EC_N0_TENSOR; NOT_FULL_INTERACTING_PARENT_BACKGROUND_RESPONSE_HIGHER_LOOPS_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/projector.md",
            "notes/dimensional.md",
            "notes/curvature.md",
            "notes/reference.md",
            "notes/stress.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_projector_and_covariant_finite_parts": serialize(
            {
                mod.__name__.rsplit(".", 1)[-1]: payload(mod.data())
                for mod in (bloch, subtraction, curvature, reference)
            }
        ),
        "absolute_curved_stress": serialize(payload(stress.data())),
        "observable_and_source_context": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All42 exact curved free Dirac copies are retained. Complete twenty-frame/state remainders are below1e410 in energy and1e416 in pressure. Full dimensional local counterterms restore a common covariant finite action. The explicit E_D,C_D prescription retains the finite Euler tensor; the independently fixed quadratic Newton reference removes its zero-field second-order local contribution only. Together with the entire fixed-reference potential, both absolute rest-frame components are below1e789 uniformly in real time, with named kappa ratios below1e-11. This is not the full interacting parent tensor, controlled background or response, physical cutoff, V/G/B or original P8 closure.",
        "not_established": [
            "Complete interacting/canonically re-expressed scalar, heavy, gauge or gravitational parent tensor and Newton dictionary",
            "Higher-loop or physical-cutoff remainders, quantum target matching, scalar source or controlled background/response",
            "Finite-gravity Regge/IR remainder, full vacuum contour/truncation, V/G/B or original P8 closure",
            "Every-boost uniform stress bound or formalization of external analytic theorems",
        ],
        "verification_boundary": "Exact projector, full dimensional local/curvature variation, finite-reference and complete radial identities support explicit written estimates. Independent finite-angle/full-function projectors, nonzero-regulator quadratures, curvature lapse/scale variations and high-precision full-profile controls supplement these proofs, not validated integration. Every own payload uses the unchanged exact serializer. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The specified curved Dirac stress report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.171.ABSOLUTE_CURVED_FREE_DIRAC_ONE_LOOP_TENSOR_WITH_EXPLICIT_EVANESCENT_AND_NEWTON_REFERENCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
