"""Read-only independent hard-radiative curvature matching report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_dimensional_gravity_radiation import verify as input1
from p8_vacuum_affine_dimensional_real_remainder import verify as input0
from p8_vacuum_affine_heavy_parent_one_loop import verify as input2
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as input3

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-radiative-curvature-matching.json"
)
PARENT_SHA = "f1b5b057f0877a8722d63cb560cbe7507ca3e9e639455d121f574cd6034b3b62"
TREE_SHA = "ce0ed27669fe9fefeb0064deae3a8565c2790ce517dae44f0737077f2dfc810b"
FLAT_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_radiative_curvature_matching/*.py"))
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
        (input1, TREE_SHA),
        (input2, FLAT_SHA),
        (input3, BORN_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen radiative curvature-matching input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_335_complete_tree_real_minus_soft_limit_rebuilt": PARENT_SHA,
        "S6_334_full_original_tree_and_polarization_tensor_rebuilt": TREE_SHA,
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
        raise ValueError("A radiative curvature-matching proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.336.INDEPENDENT_HARD_RADIATIVE_CURVATURE_MATCHING_DIRECTION",
        "date": "2026-09-18",
        "status": "SCOPED_RADIATIVE_MATCHING_NONIDENTIFIABILITY; NOT_MATCHED_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/operator.md",
            "notes/matching.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_covariant_contact_and_calibrations": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_flat_matching_obstruction_and_conditional_bounds": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "An explicit Weyl-curvature counterterm direction is invisible to retained flat matching through one loop and unchanged leading/subleading/sub-subleading soft data, but gives a nonzero Bose-symmetric conserved physical four-scalar one-graviton contact proportional to chi*omega^2. Exact original47-tree interference is nonzero on three physical states. Conditional finite-rate bounds do not fix chi. No physical coefficient is selected; the complete hard radiative functional and quantum matching remain open.",
        "not_established": [
            "Physical value, sign or bound of the curvature matching coefficient",
            "Complete covariant one-loop hard radiation and curved counterfunctional",
            "An exhaustive matching basis or arbitrary-coefficient UV completion",
            "Full interacting detector probability, quantum unitarity and complex Regge",
            "Same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Literal general curvature and all256 pure-gauge components, all24 labelled scalar assignments, exact recoil and both polarizations at three original-source states, topological loop counting and written uniform conditional bounds. The known generic soft obstruction is credited to primary literature without importing its D>=5 loop theorem into D4. Not kernel-formalized; frozen prior evidence unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The radiative curvature-matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.336 independent hard radiative curvature-matching replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
