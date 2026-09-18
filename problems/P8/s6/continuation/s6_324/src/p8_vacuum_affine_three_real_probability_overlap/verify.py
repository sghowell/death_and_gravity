"""Read-only complete three-real probability-overlap certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_three_singleton_subtraction import verify as previous
from p8_vacuum_affine_two_real_soft_overlap import verify as analytic_input
from p8_vacuum_affine_uniform_all_tree_bound import verify as uniform_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-three-real-probability-overlap.json"
)
PARENT_SHA = "994a59e6ff556a08965d1c5386cd7aed1d484f017dd64a181bfa225c105ddcdc"
ANALYTIC_SHA = "79539dcf68142c2ebb62a0236293be98d2f1b6eece93d49e46cbf091fb8e1994"
UNIFORM_SHA = "6239bb0b855c48275f8598cdfdec3f8dfa2590545d5b4abb5b6ddfe765017b0b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_three_real_probability_overlap/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (previous, PARENT_SHA),
        (analytic_input, ANALYTIC_SHA),
        (uniform_input, UNIFORM_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen probability-overlap proof input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_323_complete_three_real_amplitude_and_signed_measure_rebuilt": PARENT_SHA,
        "S6_313_uniform_two_real_derivative_theorem_rebuilt": ANALYTIC_SHA,
        "S6_319_complete_tree_baseline_bound_rebuilt": UNIFORM_SHA,
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
        raise ValueError("A complete probability-overlap proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.324.COMPLETE_THREE_REAL_PROBABILITY_OVERLAP_AND_FINITE_SCHEME_TRANSFER",
        "date": "2026-09-17",
        "status": "SCOPED_THREE_REAL_PROBABILITY_RECTANGLE_AND_SIGNED_TRANSFER; NOT_INCLUSIVE_ALL_N_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/faces.md",
            "notes/kinematics.md",
            "notes/projection.md",
            "notes/measure.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_external_residues": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_proper_faces_and_probability_transfer": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The correct complete three-real soft faces include both remaining hard null legs. Literal scalar/EH3 residues and a kinematic-current-only absolute disc give uniform first and mixed proper-face bounds. Exact Boolean overlap algebra retains all12 additional proper-component interference terms. Their integrable finite transfer yields the complete anchored probability signed measure with total variation below2e-1424*x^2+1e-1340*x^4 and a common-cutoff limit. This does not match actual virtual counterterms or prove a positive normalized inclusive rate, all-N summation, state, Regge, bounce or original P8 closure.",
        "not_established": [
            "Actual real-virtual and integrated-counterterm matching",
            "A positive normalized or monotone inclusive detector measure",
            "Complete all-N subtraction/summation and finite hard or evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Generic ten-parameter soft tensor, both hard TT polarizations and general conserved-root scalar/EH3 residue identities; original434/5116 two-marked hierarchical calibrations; explicit absolute-radius kinematic recoil/Doppler margins and frozen two-real Cauchy-product bounds; all27 probability union pairs and12 proper overlaps, exact complex interferences, entropy product integrals and phase-space factors; nonzero omitted-null, omitted-overlap and radiative-state matching controls. Uniformity comes from the written proofs, not finite samples or integer gates alone. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete probability-overlap report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.324 complete three-real probability overlap replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
