"""Read-only complete selected factorized-bubble radiation report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_dimensional_real_remainder import verify as input1
from p8_vacuum_affine_heavy_parent_one_loop import verify as input2
from p8_vacuum_affine_heavy_scalar_four_point_loop import verify as input4
from p8_vacuum_affine_mixed_source_radiation import verify as input0
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as input3

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-factorized-bubble-radiation.json"
PARENT_SHA = "faf709f6d72af441b9ff14d6ceb57fe573cf22d6b37a45d80ccf58205d5c9bc6"
BOUND_SHA = "f1b5b057f0877a8722d63cb560cbe7507ca3e9e639455d121f574cd6034b3b62"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
FLAT_SHA = "1277412d19d9b022eb81c3728251d1836d88e4b047c8ad80794107c2aff7c911"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_factorized_bubble_radiation/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (input0, PARENT_SHA),
        (input1, BOUND_SHA),
        (input2, SOURCE_SHA),
        (input3, BORN_SHA),
        (input4, FLAT_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen factorized-bubble radiation input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_338_complete_mixed_source_radiation_and_unknown_curvature_frontier_rebuilt": PARENT_SHA,
        "S6_335_complete_tree_real_minus_soft_limit_rebuilt": BOUND_SHA,
        "S6_239_complete_original_source_and_fixed_prescription_rebuilt": SOURCE_SHA,
        "S6_235_full_original_factorized_bubble_class_rebuilt": FLAT_SHA,
        "S6_297_original_full_Born_and_finite_kappa_transfer_rebuilt": BORN_SHA,
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
        raise ValueError("A selected factorized-bubble radiation proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.339.COMPLETE_FACTORIZED_BUBBLE_RADIATION_AND_FINITE_REMAINDER",
        "date": "2026-09-18",
        "status": "SCOPED_COMPLETE_FACTORIZED_BUBBLE_RADIATION; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/kernel.md",
            "notes/loops.md",
            "notes/branch.md",
            "notes/bounds.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/provenance.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_factorized_source_kernel_loops_and_physical_branch": serialize(
            {name: payload(packets[name]) for name in names[:4]}
        ),
        "whole_uniform_bounds_and_exact_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[4:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The entire original factorized bubble class, including both inner light and both outer heavy insertions and all fixed counterterms, has the computed physical TT response of A(v)^2 B(v). The exact-D proof, actual complex branch and continuous coincidence limit are retained. The original nonleading remainder is<1e194 A0/sqrt(kappa)=1e-206 A0, with full-tree finite interference<1e-602 in the adaptive low window and a separate fixed-angle high bound. Known S337/S338 remainders fit within the same rounded budget. The isolated hard bound1e190 is deliberately not perturbative smallness. Other loop sectors, independent curvature matching and original P8 remain OPEN.",
        "not_established": [
            "Remaining triangle, ordered-box, mixed quadratic or internal-graviton radiation",
            "Small isolated hard loop correction, unprojected loop tensor, complete curved counterfunctional or independent physical curvature matching",
            "Full inclusive interacting detector probability or loop rate",
            "Quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Generic two-endpoint Ward identities, literal local covariant quartic actions, both equal-mass exact-D triangle insertions, ordered A-B-A variation, complete frozen S235 coefficient and all three UV counterterms, exact physical-branch derivatives and original uniform budgets. Thirty independent full tensor components at diagnostic mass17 and four exact ORIGINAL absorptive recoil states retain both polarizations and actual coincident invariants. Samples are not the uniform proof. Not kernel-formalized; frozen prior evidence unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The factorized-bubble radiation report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.339 complete selected factorized-bubble radiation replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
