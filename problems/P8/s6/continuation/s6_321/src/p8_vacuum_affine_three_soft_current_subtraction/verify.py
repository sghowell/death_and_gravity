"""Read-only three-current and one-block subtraction certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_relative_energy_complex_tube import verify as previous
from p8_vacuum_affine_two_real_soft_overlap import verify as analytic_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-three-soft-current-subtraction.json"
)
PARENT_SHA = "18292735783403d7588808de71e7620238ffbbe2f70d82a5a1012ef431f03108"
ANALYTIC_SHA = "79539dcf68142c2ebb62a0236293be98d2f1b6eece93d49e46cbf091fb8e1994"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_three_soft_current_subtraction/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen relative-energy complex-tube parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(analytic_input.REPORT) != ANALYTIC_SHA:
        raise ValueError("The frozen complex-recoil and two-real overlap input changed")
    analytic_input.validate_report(
        json.loads(analytic_input.REPORT.read_text()), analytic_input.build_report()
    )
    return {
        "S6_320_relative_energy_complex_tube_rebuilt": PARENT_SHA,
        "S6_313_complex_recoil_and_two_real_overlap_rebuilt": ANALYTIC_SHA,
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
        raise ValueError("A three-current subtraction proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.321.UNIFORM_THREE_SOFT_CURRENT_FACES_AND_ONE_BLOCK_SUBTRACTION",
        "date": "2026-09-17",
        "status": "SCOPED_THREE_CURRENT_AND188_CLASS_SUBTRACTION; NOT_FULL_INCLUSIVE_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/current.md",
            "notes/geometry.md",
            "notes/derivatives.md",
            "notes/faces.md",
            "notes/core.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_angular_coverage": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_three_current_and_one_block_subtraction": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete four-tree temporal soft current has an all-angle squarefree energy derivative hierarchy, compatible continuous energy faces and rectangle bound721*J/kappa. A global total-energy tube is proved only for the amputated47 hard functional, giving a below1e-738*J rectangle for the188-term maximal-soft-block class. The positive-coefficient certificates reconstruct all16 sector/TT cases and768 independent Taylor coefficients. This does not subtract the other4928 temporal terms or establish an inclusive rate, quantum state, Regge or original V/G/B/P8 closure.",
        "not_established": [
            "The subtracted remainder for the other4928 temporal terms or the complete5116-tree probability",
            "The full overlapping all-N soft subtraction and real-virtual completion",
            "Complete finite hard real-virtual and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Literal polynomial EH3/EH4 action, complete source Ward identities, exact known-factor division and factored assembly, positive denominator/AM-GM coefficient majorants, independent nilpotent Taylor comparisons, exact joint-axis compatibility and uniform origin control. Symbolic angular coverage plus complex-TT multilinearity; one-block hard-core topology and closed-disc gaps, original-parameter47/188 calibrations, operator Cauchy, eight-term product rule and entropy-kernel integrability. Finite points supplement the written uniform proof; no private current cache or unproved full-amplitude global-W holomorphy is used. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The three-current subtraction report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.321 three-current subtraction replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
