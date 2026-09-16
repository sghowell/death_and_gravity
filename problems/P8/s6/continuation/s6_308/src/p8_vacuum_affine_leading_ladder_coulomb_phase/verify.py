"""Read-only full-D leading-ladder Coulomb phase certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_finite_forward_phase import verify as phase_input
from p8_vacuum_affine_uniform_radiation_soft_limit import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-leading-ladder-coulomb-phase.json"
)
PARENT_SHA = "4238712f24c9953227701d0eb77c817e83cb86f513c4c7d25a46f648931034c7"
PHASE_SHA = "d02d6507a596c20adad3d3fd5a654f1a8abdca1857a1a9b49bc4c91bc5f2457a"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_leading_ladder_coulomb_phase/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen uniform single-real soft-limit parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(phase_input.REPORT) != PHASE_SHA:
        raise ValueError("The frozen complete finite forward Coulomb phase changed")
    phase_input.validate_report(
        json.loads(phase_input.REPORT.read_text()), phase_input.build_report()
    )
    return {
        "S6_307_complete_single_real_uniform_soft_limit_rebuilt": PARENT_SHA,
        "S6_303_complete_one_loop_Coulomb_phase_rebuilt": PHASE_SHA,
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
        raise ValueError("A leading-ladder Coulomb proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.308.FULL_D_LEADING_LADDER_COULOMB_PHASE_AND_UNIFORM_GAMMA_CORRECTION",
        "date": "2026-09-16",
        "status": "SCOPED_EXACT_LEADING_LADDER_COULOMB_PHASE; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/dimensional.md",
            "notes/resummation.md",
            "notes/fourier.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_full_D_ladders": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_exact_phase_and_uniform_Gamma_bound": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged full-D leading pure-gravity ladder class has an exact fixed-order infrared subtraction and all-rung finite generating function. A finite-difference proof retains the dimensional2/V correction and reproduces S303's entire one-loop principal phase. The closed gamma-function Coulomb factor has unit physical-edge modulus and differs from its unexpanded logarithmic phase by less than2e-2400 at original parameters, uniformly in transfer. This is not a bound on the complete original amplitude minus that leading class; all physical matching, nonleading quantum/real and Regge obligations and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A quantitative remainder for the full interacting amplitude beyond the leading ladder class",
            "Finite Newton, curvature and other physical matching coefficients",
            "Complete all-N nonleading radiation, full detector rate or a finite forward cross section",
            "High-energy complex Regge control or a UV quantum-state/bounce construction",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact full-D Fourier normalization, general finite-difference cancellation at arbitrary fixed order, convergent finite leading-class generating function, independent one-loop phase and Gaussian-Abel Fourier calibration, and complete gamma-tail bounds. Neither a classical macroscopic hierarchy nor an uncontrolled infinite-sum/UV-limit interchange is imported. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. Frozen ancestors, original source, matching and scoped P8(a) remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The leading-ladder Coulomb report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.308 leading-ladder Coulomb replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
