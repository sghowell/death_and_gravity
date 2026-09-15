"""Read-only full-source homogeneous fold consistency and physical endpoint report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_global_lapse_chart_obstruction import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-homogeneous-fold-dynamics.json"
PARENT_SHA = "d24260483929473a16611060d16ee19166c17f749d6c85cb2385463ddf3dc5ee"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_homogeneous_fold_dynamics/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen original lapse-chart obstruction report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_276_full_lapse_chart_and_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A full homogeneous fold proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.277.COMPLETE_FULL_SOURCE_HOMOGENEOUS_LAPSE_FOLD_CONSTRAINT_PRESERVATION_OBSTRUCTION_AND_REGULAR_CLASSICAL_FINITE_PROPER_TIME_CURVATURE_ENDPOINT",
        "date": "2026-09-15",
        "status": "SCOPED_HOMOGENEOUS_CLASSICAL_FOLD_DYNAMICS; NOT_ORIGINAL_STATE_OR_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/canonical.md",
            "notes/consistency.md",
            "notes/endpoints.md",
            "notes/geometry.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_clock_source_and_canonical_flow": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_constraint_preservation_and_physical_endpoint": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Full original time differentiation retains the nonzero primitive/B contact. Exact homogeneous canonical brackets give first fold consistency K=p[-3/25+3(P-rho)]+3P_u/2-rho_u. If K vanishes, actual fixed C5 profiles force |p|<21 epsilon; the complete second-preservation quadratic has a>28, |b|<1, c>60 and discriminant<-6719. No C1 lapse solution passes through this entire specified homogeneous fold family. For p=-1/10 the smooth full-source desingularized constraint flow supplies regular C_N<0 classical solutions that approach the fold in finite proper time. The physical lapse remains N; normal expansion and the correctly reconstructed physical Ricci scalar diverge while lapse and volume have positive finite limits. These are separately specified homogeneous classical data, not the original prepared quantum state, the S276 finite-wave path or the corrected S275 finite hybrid. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "An interacting quantum state, rank-changing quantum measure, self-adjoint domain or quantum singularity resolution",
            "A singularity or generic instability of the unchanged S275 finite hybrid or the original prepared state",
            "A time realization of the different S276 mean-preserving wave path or classification of all inhomogeneous folds",
            "An evaluated initial interval, uniform near-bounce domain, physical cutoff or regulator removal",
            "Original physical UV matching, omitted-loop/Regge control, all-parent no-go, global P8 closure or formal proof verification",
        ],
        "verification_boundary": "Exact entire-source derivatives, actual fixed-profile enclosures, literal canonical brackets and independent physical-metric Ricci diagnostics support written local ODE and endpoint proofs. Native/direct/ordinary/CLI use original SymPy; only a captured complete FULL regression may use the audited exact-GCD adapter. No frozen source or report is rewritten, and prior physical qualifications remain explicit.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete homogeneous fold dynamics report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.277 full-source homogeneous fold consistency and physical endpoint replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
