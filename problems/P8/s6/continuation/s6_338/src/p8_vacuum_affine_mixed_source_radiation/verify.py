"""Read-only complete selected mixed-source radiation report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_dimensional_real_remainder import verify as input1
from p8_vacuum_affine_heavy_parent_one_loop import verify as input2
from p8_vacuum_affine_local_tadpole_radiation import verify as input0
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as input3

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-mixed-source-radiation.json"
PARENT_SHA = "602dd0e80cc5fec993a20219f3423195439b0e3dba03e468794bc939b2627645"
BOUND_SHA = "f1b5b057f0877a8722d63cb560cbe7507ca3e9e639455d121f574cd6034b3b62"
FLAT_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_mixed_source_radiation/*.py"))
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
        (input2, FLAT_SHA),
        (input3, BORN_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen mixed-source radiation input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_337_complete_local_tadpoles_and_unassigned_curvature_frontier_rebuilt": PARENT_SHA,
        "S6_335_complete_tree_real_minus_soft_limit_rebuilt": BOUND_SHA,
        "S6_239_complete_flat_source_loop_prescription_rebuilt": FLAT_SHA,
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
        raise ValueError("A selected mixed-source radiation proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.338.COMPLETE_SELECTED_MIXED_SOURCE_RADIATION_AND_FINITE_REMAINDER",
        "date": "2026-09-18",
        "status": "SCOPED_COMPLETE_MIXED_SOURCE_RADIATION; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/within.md",
            "notes/cross.md",
            "notes/kernel.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/provenance.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_contractions_and_exact_D_kernel": serialize(
            {name: payload(packets[name]) for name in names[:4]}
        ),
        "whole_complete_radiation_bounds_and_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[4:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All three selected J2 G_H J4 loop-radiation classes are retained. The within-J2 and across-source physical TT classes cancel under the existing conditions, as shown by full source/onepoint, both exact-D triangle and labelled EOM calculations. The within-J4 sector gives Q times centered complete matter radiation,Q=8(4-n)/(16pi^2 kappa), including all42 source graphs and the fixed OS4 contact. Its nonleading original-source remainder is<1e-598 A0/sqrt(kappa), with finite full-tree interference<1e-1394 in the adaptive low window and a separate fixed-angle high bound. Adding S337 gives known local-plus-source bounds1e-597 and1e-1393. Other hard loops, independent curvature matching and original P8 remain OPEN.",
        "not_established": [
            "Remaining polynomial/heavy or internal-graviton radiative loop sectors",
            "Complete curved counterfunctional or independent physical curvature matching",
            "Full inclusive interacting detector probability or loop rate",
            "Quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Literal24 source labels, exact-D tadpole and both massive kernel triangles, covariant Green-function delta/EOM terms, full24 EOM emissions, generic traceful Ward and physical42-graph reductions, exact original contact and convexity identities, both polarizations at three original recoil states, and written uniform budgets. Finite samples are not the uniform proof. Not kernel-formalized; all frozen prior scientific evidence unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The mixed-source radiation report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.338 complete selected mixed-source radiation replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
