"""Read-only global one-loop insertion and finite reference-window bounds."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_full_one_loop_matching import verify as parent

from . import audit, calibration, halfplane, routing, spacelike, window

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-global-one-loop-insertions.json"
PARENT_SHA = "39613ad0083b068662d7b14475aaf3054db2631dfe142f4ca81343398701805b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_global_one_loop_insertions/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full one-loop reference report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_133_fully_rebuilt": PARENT_SHA,
        "same_complete_scalar_fermion_one_loop_reference": True,
        "whole_prior_actions_and_reports_unchanged": True,
        "not_outer_two_loop_or_original_P8_closure": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A global one-loop insertion gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.134.GLOBAL_ONE_LOOP_INTERNAL_MOMENTUM_CONTROL",
        "date": "2026-09-10",
        "status": "GLOBAL_ONE_LOOP_INSERTION_AND_FINITE_REFERENCE_WINDOW_CONTROL; NOT_OUTER_TWO_LOOP_HIGHER_LOOP_ERROR_CUTOFF_FULL_SPECTRUM_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/halfplane.md",
            "notes/spacelike.md",
            "notes/routing.md",
            "notes/window.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "whole_kernel_halfplane_insertion": serialize(payload(halfplane.data())),
        "global_spacelike_decay_envelope": serialize(payload(spacelike.data())),
        "selected_all_momentum_routing_and_window": serialize(
            {"routing": payload(routing.data()), "window": payload(window.data())}
        ),
        "actual_reference_window_and_insertion_bounds": serialize(
            payload(calibration.data())
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete one-loop fermion kernel, with all affine reference terms and finite mass/residue anchors subtracted first, admits a uniform insertion bound on an unbounded complex half-plane. A separate exact spacelike logarithmic envelope retains high-momentum decay. Every real loop momentum in the named symmetric bubble routing lies inside the certified half-plane on the unit forward disc. Combining both scalar and fermion one-loop sectors and the already fixed canonical normalization bounds the selected Euclidean inverse defect below 10^-202 on the explicit reference-energy window through 10^400. That mathematical window is not a derived cutoff. This is an ingredient for subsequent loop integration, not a complete outer two-loop forest or a new-model higher-loop error bound, exact Lorentzian/gauge spectrum, reflection positivity, V contour, finite-gravity G, common-parent B or original P8 closure.",
        "not_established": [
            "Complete outer two-loop forests, integrals and primitive contributions",
            "New-model higher-loop and truncation error bounds",
            "Derived Wilsonian cutoff or global quantum potential and exact spectrum",
            "Reflection positivity or a full Lorentzian pole classification",
            "Justified high-energy V contour and strict full V verdict",
            "Finite-gravity regulated Regge Delta for G and common bounce-parent B",
            "Original P8 closure or exclusion of any complete old ladder row",
        ],
        "verification_boundary": "Independent exact differentiation, frozen-kernel agreement, resolvent/logarithm identities, all-real-momentum routing algebra, rational calibration and domain controls. Half-plane analyticity, parameter majorants and Taylor remainders are written analytic proofs, not proof-assistant formalization or all-orders QFT claims. Native and ordinary scientific replay uses unmodified SymPy with interpreter recursion and trusted-integer formatting allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The global one-loop insertion report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.134.GLOBAL_ONE_LOOP_INTERNAL_MOMENTUM_CONTROL replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
