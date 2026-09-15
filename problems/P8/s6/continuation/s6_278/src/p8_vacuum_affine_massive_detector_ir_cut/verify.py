"""Read-only original massless-M1 cut and massive detector dictionary report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_homogeneous_fold_dynamics import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-massive-detector-ir-cut.json"
PARENT_SHA = "5fa5971b4a135548020426102c8cce3a7edbc56edb04a019a2a6c8f50e46b785"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_massive_detector_ir_cut/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen original homogeneous-fold report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_277_full_homogeneous_fold_and_entire_ancestry_rebuilt": PARENT_SHA,
        "S275_corrected_finite_hybrid_not_refuted": True,
        "S276_wave_path_not_reidentified_as_this_classical_solution": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A massive detector IR cut proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.278.ORIGINAL_MASSLESS_M1_COMPLETE_FIRST_GRAVITATIONAL_PAIR_CUT_AND_MASSIVE_SOFT_STRIPPING_WITH_EXPLICIT_CROSSING_SYMMETRIC_LOW_CUT_SUBTRACTION",
        "date": "2026-09-15",
        "status": "SCOPED_FORMAL_M1_LOOP_AND_CONDITIONAL_SOFT_DICTIONARY; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/cut.md",
            "notes/soft.md",
            "notes/subtraction.md",
            "notes/dispersion.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_massive_vacuum_and_massless_pair_cut": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_soft_resolution_and_crossing_subtraction": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original minimally coupled massless M1 is free only in the classical gravitational limit. At finite kappa its complete first pair cut in massive four-Phi scattering is [s^2-tu+2m^2s+6m^4]/(960pi kappa^2), with all crossed logarithms retained modulo real local terms. In particular a nonzero transfer logarithm remains after universal soft-graviton stripping. The equal-mass analytic soft dictionary gives a conditional resolution law, and an explicit three-channel low-cut subtraction gives an analytic center for this loop sector with a cap-dependent coefficient and a separate transfer term. Other channels, finite local matching, full detector unitarity and Regge errors remain unresolved. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A full interacting finite-gravity amplitude, quantum vacuum or exact physical mass and LSZ matching",
            "All graviton and other-species cuts, finite local counterterms or omitted loops",
            "An exact finite-detector positivity bound, positive UV smearing or explicit full Regge remainder",
            "A Wilsonian cutoff inferred from the subtraction cap, quantum bounce or original-state completion",
            "Original V/G/B/P8 closure, all-parent exclusion or formal verification of the written analysis",
        ],
        "verification_boundary": "Exact full-source low jets, four-vector Ward and stress contractions, independent angular and partial-wave sewing, massive parameter-integral algebra and exact crossing subtraction support written one-loop and analytic arguments. External soft-factor and detector-dispersion premises are explicitly conditional; massless-scalar numerical bounds are not imported. Native/direct/ordinary/CLI retain original SymPy; only a captured complete FULL replay uses the audited exact-GCD adapter. Prior frozen science and historical physical qualifications remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete massive detector IR cut report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.278 original M1 gravitational cut and massive detector dictionary replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
