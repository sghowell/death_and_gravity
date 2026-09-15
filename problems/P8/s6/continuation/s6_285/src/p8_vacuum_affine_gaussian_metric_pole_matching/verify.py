"""Read-only complete fixed gapped Gaussian metric matching report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_dimensional_cut_completion import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-gaussian-metric-pole-matching.json"
)
PARENT_SHA = "a855d4c9a602541cbce21f1432f0e6fec42c4d565a15922094c172d67f39971b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_gaussian_metric_pole_matching/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen original dimensional cut and Gram-principal-part report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_284_dimensional_cuts_Gram_parts_and_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A fixed Gaussian metric matching proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.285.ORIGINAL_THREE_GAUSSIAN_METRIC_CUTS_FIXED_H_PROCA_CURVATURE_MATCHING_AND_EXPLICIT_UNMATCHED_LIGHT_GRAVITATIONAL_POLYNOMIAL",
        "date": "2026-09-15",
        "status": "SCOPED_FIXED_H_PROCA_METRIC_MATCHING_WITH_LIGHT_CURVATURE_OPEN; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/spectral.md",
            "notes/matching.md",
            "notes/volume.md",
            "notes/poles.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_Gaussian_metric_cuts": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_fixed_Gaussian_polynomial_and_metric_poles": serialize(
            {
                "packets": {name: payload(packets[name]) for name in names[2:]},
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All three fixed Gaussian volume constants cancel the whole cosmological metric density. Literal scalar and Proca stress cuts give the complete conserved nonlocal metric sectors. The prescribed H and Proca covariance and covariant finite matching determine their curvature polynomial; their signature-correct finite Newton shift is positive and relatively below10^-602, with a positive massless residue and no additional low-disk pole in this prescribed subkernel. The light scalar finite gravitational curvature polynomial is not fixed by its vacuum constant or by the H prescription. Its three unresolved coefficients remain explicit in the full Phi-inclusive kernel and massless residue. Complete four-point local/tadpole, proper-vertex/LSZ, massless IR/Regge, exact original-state and bounce obligations remain. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Light-scalar finite gravitational curvature matching and full Phi-inclusive pole or inverse bounds; complete four-scalar tadpole/local/threshold matching",
            "Full proper scalar-graviton vertex, external LSZ, physical mass and exact interacting Newton normalization",
            "Full finite crossed amplitude or b20, massless and omitted interacting loops, finite IR detector/Regge control",
            "A physical EFT cutoff, exact interacting vacuum or original curved Gaussian/state matching",
            "Original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole original source, both fixed covariant finite prescriptions, literal scalar and Proca conserved stress cuts, arbitrary canonical covariance and second vertex, complete finite metric variations, exact convergent moments and low-disk inequalities. Written analytic/covariant matching and positivity arguments are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; the exact-GCD adapter is restricted to a captured complete FULL regression. Frozen sources, raw reports, all fixed constants and historical qualifications are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The fixed Gaussian metric matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.285 fixed H/Proca metric matching and explicit light-curvature frontier replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
