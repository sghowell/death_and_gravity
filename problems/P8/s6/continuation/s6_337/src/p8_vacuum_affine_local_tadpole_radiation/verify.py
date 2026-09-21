"""Read-only complete selected local-tadpole radiation report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_dimensional_real_remainder import verify as input1
from p8_vacuum_affine_heavy_parent_one_loop import verify as input2
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as input3
from p8_vacuum_affine_radiative_curvature_matching import verify as input0

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-local-tadpole-radiation.json"
PARENT_SHA = "f8f9cd57c0240fe6008efd743d399c27dfc393a14d50501e8362289aae204db0"
BOUND_SHA = "f1b5b057f0877a8722d63cb560cbe7507ca3e9e639455d121f574cd6034b3b62"
FLAT_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_local_tadpole_radiation/*.py"))
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
            raise ValueError("A frozen local-tadpole radiation input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_336_independent_curvature_matching_remains_unassigned_rebuilt": PARENT_SHA,
        "S6_335_complete_tree_real_minus_soft_limit_rebuilt": BOUND_SHA,
        "S6_239_complete_flat_one_loop_prescription_not_curved_matching_rebuilt": FLAT_SHA,
        "S6_297_original_full_Born_and_flat_graph_transfer_rebuilt": BORN_SHA,
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
        raise ValueError("A selected local-tadpole radiation proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.337.COMPLETE_SELECTED_LOCAL_TADPOLE_RADIATION_AND_FINITE_REMAINDER",
        "date": "2026-09-18",
        "status": "SCOPED_COMPLETE_LOCAL_TADPOLE_RADIATION; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/contractions.md",
            "notes/radiation.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/provenance.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_D_contractions_and_radiation": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_uniform_bounds_and_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete six-local-degree-six tadpole radiation sector is rebuilt from all720 assignments, exact-D loop insertion and metric tadpole, including covariant connection terms. Finite evanescence, all four external emissions and generic Ward/reduction identities yield the known OS4 local TT representative. Its original nonleading remainder is <1e-791 A0/sqrt(kappa); the finite full-tree interference is <1e-1587 in the adaptive low window, with a separate fixed-angle high-window bound. Other loop sectors, independent curvature matching and original P8 remain OPEN.",
        "not_established": [
            "Other source/matter, polynomial/heavy or internal-graviton loop radiation",
            "Complete curved counterfunctional or physical curvature matching",
            "Full inclusive interacting loop rate or detector probability",
            "Quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Literal720 six-field and24 quartic assignments; off-shell flat Wick factors; general on-shell Gram exact-D 1PI equality and independent connection comparison; eight unrestricted Ward identities and five generic physical reductions; exact original couplings, both polarizations at three recoil states and written uniform budgets. Finite calibrations are not a substitute for the analytic bounds. Not kernel-formalized; all frozen prior evidence unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The local-tadpole radiation report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.337 complete selected local-tadpole radiation replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
