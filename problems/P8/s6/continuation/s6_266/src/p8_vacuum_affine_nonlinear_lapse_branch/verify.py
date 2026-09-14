"""Read-only full nonlinear auxiliary-domain certificate and unchanged ancestry."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_nonlinear_reference_volume import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-nonlinear-lapse-branch.json"
PARENT_SHA = "7763eb35f0b2e6176bd1d5e6f25248e144c5ebc454be8c7ca706bb1365bcb349"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_nonlinear_lapse_branch/*.py"))
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
            "The frozen nonlinear-reference and coefficient-boundary report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_265_full_reference_composites_and_whole_ancestry_rebuilt": PARENT_SHA,
        "same_original_R_F_profiles_source_and_full_Hamiltonian": True,
        "original_state_and_H_Proca_preparations_unchanged": True,
        "S265_Gaussian_coefficient_boundary_and_S261_refutation_unchanged": True,
        "classical_nonlinear_branch_not_quantum_measure_or_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A whole-source nonlinear-root or response gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.266.COMPLETE_SOURCE_PINNED_QUANTITATIVE_NONLINEAR_BOUNCE_SLICE_LAPSE_NORMAL_VECTOR_BRANCH_WITH_ALL_SPATIAL_MATTER_CHANNELS_AND_FULL_IMPLICIT_PHYSICAL_RESPONSE_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "FULL_LOCAL_NONLINEAR_AUXILIARY_BRANCH_AND_PHYSICAL_CHAIN; NOT_CANONICAL_GAUSSIAN_SUPPORT_INTERACTING_MEAN_QUANTUM_REGULATOR_CUTOFF_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/domain.md",
            "notes/response.md",
            "notes/physical.md",
            "notes/measure.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_full_source_nonlinear_domain_and_auxiliary_pivots": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_full_implicit_and_physical_composite_responses": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "At u=0 the complete original Hamiltonian, including its nonzero primitive time contact and every retained matter/spatial/vector channel, admits a unique quantitative positive-lapse and normal-vector auxiliary branch on the specified twelve-invariant box. A1/30 contraction gives the full lapse within1e-245 of1, with source-only center shift below1e-375, both auxiliary pivots uniformly regular, and all complete coefficients bounded on that branch. Full first/second/third implicit responses and physical-volume chain contacts retain nonlinear canonical-invariant corrections. The invariant coordinates are not Gaussian canonical modes. The original preparations, sources, masses and physical map are unchanged. S265's distinct Gaussian coefficient-integrability obstruction and S261's refutation remain. No interacting quantum mean, ordering/continuum-gauge regulator, physical cutoff, nonlinear inhomogeneous evolution or original V/G/B/P8 closure follows.",
        "not_established": [
            "A Gaussian law, support or quantum projection for the nonlinear invariant or auxiliary coordinates",
            "A self-adjoint positive physical operator from a bounded classical Weyl symbol alone",
            "An interacting original physical quantum mean, its counterterms, state or omitted-loop bounds",
            "A continuum/spatial first-class regulator, full gauge determinant or vanishing Ward defects",
            "A physical Wilsonian matching cutoff or closure of original vacuum or finite-gravity gates",
            "A time-dependent nonlinear inhomogeneous Cauchy, stability or global geodesic-completeness theorem",
            "Global auxiliary-root uniqueness outside the explicit strip",
            "A refutation of the S265 individual Gaussian-coefficient result or restoration of S261's false map",
            "Formal verification of the written interval, implicit-function or physical-domain arguments",
        ],
        "verification_boundary": "Exact checks bind the full original Hamiltonian, complete primitive derivative, actual mass and source-C5 envelopes. Outward interval bounds certify a genuine local nonlinear root and auxiliary pivots. Independent finite-source root solves, Hamiltonian variations and mixed response checks are diagnostics, not changed actual P8 parameters. Density invariants retain their nonlinear canonical maps. Native/direct/ordinary/CLI use original SymPy; only the captured full regression uses the independently audited exact-GCD adapter. Every frozen predecessor and both historical errata remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete nonlinear auxiliary-branch report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.266 full nonlinear auxiliary branch and physical-chain replay passed; quantum domain, interacting mean, cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
