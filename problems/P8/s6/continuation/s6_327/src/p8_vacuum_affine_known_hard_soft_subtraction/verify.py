"""Read-only selected known-hard soft and sharp logarithmic report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_finite_forward_phase import verify as hard_input
from p8_vacuum_affine_minimal_gravity_radiation import verify as real_input
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-known-hard-soft-subtraction.json"
PARENT_SHA = "cb310adf6525b22e6f89eadd2d9badbbfdf1522016831c27d2fedcea9bf806fe"
HARD_SHA = "d02d6507a596c20adad3d3fd5a654f1a8abdca1857a1a9b49bc4c91bc5f2457a"
REAL_SHA = "a34d5160c55891304eda5f83bf90970854774e8aae327d9375b89e5fa797063e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_known_hard_soft_subtraction/*.py"))
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
        (hard_input, HARD_SHA),
        (real_input, REAL_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen selected hard-soft input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_326_full_radiative_logarithmic_coefficient_rebuilt": PARENT_SHA,
        "S6_303_known_hard_loop_real_interference_bound_rebuilt": HARD_SHA,
        "S6_304_complete_original47_tree_and_physical_gap_bounds_rebuilt": REAL_SHA,
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
        raise ValueError("A selected hard-soft interference proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.327.UNIFORM_KNOWN_HARD_SOFT_SUBTRACTION_AND_SHARP_LOG_INTERFERENCE",
        "date": "2026-09-18",
        "status": "SCOPED_SELECTED_KNOWN_HARD_SOFT_AND_LOG_INTERFERENCE; NOT_FULL_LOOP_ALL_N_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/phase.md",
            "notes/hard-soft.md",
            "notes/logarithm.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_selected_interference": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_uniform_signed_soft_integrals": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The known S303 hard-loop Born-soft extension has an exactly unpolarized phase cancellation in its linear interference with the complete original47-tree amplitude. The leading-pole-subtracted signed measure has TV<min(1e30,2e32*x)/kappa^2, uniformly in every nonforward hard angle. The full imaginary amplitude remains, and actual single-helicity counterexamples prevent a per-helicity deletion. Sharper original N1 estimates strengthen the named S326 logarithmic addition to TV<1e18*x*(1+ln(1/x))/kappa^2 plus its explicitly higher-order logarithmic square. These selected pieces do not supply the finite radiative hard remainder, unknown matching, common regulator, full loop square, positive detector measure, all-N quantum completion or original P8 closure.",
        "not_established": [
            "Finite radiative hard-loop remainder and full matched real-virtual rate",
            "Unknown hard coordinates and evanescent/common-regulator completion",
            "A bound on the full hard-loop modulus, square or combined squared amplitude",
            "Positive normalized detector measure, all-N hard summation or quantum unitarity",
            "Absolute complex Regge, original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact unitary helicity algebra; actual full26-matter plus21-Einstein graphs at two original-parameter kinematics with nonzero single-helicity imaginary interference; fixed original303/304 constants; analytic low/high energy estimates and exact integral/budget checks. The explicit written proof supplies uniformity and integrability, not finite samples. Only the selected fixed-Born hard-soft term and named logarithmic addition are bounded, not the whole five-point loop or regulator assembly. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The selected hard-soft report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.327 selected hard-soft and sharp logarithmic interference replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
