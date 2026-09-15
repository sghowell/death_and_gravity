"""Read-only covariant local and Gaussian four-point matching report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_gravity_pole_completion import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-covariant-gaussian-four-point-matching.json"
)
PARENT_SHA = "658ed46e1d638540bc34c2e2259c2ed3a05ec2fe539005904e2f8342d2f4a543"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob(
                "src/p8_vacuum_affine_covariant_gaussian_four_point_matching/*.py"
            )
        )
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen minimal-gravity pole-completion report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_288_whole_minimal_pole_completion_and_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A covariant local or Gaussian matching proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.289.ORIGINAL_ALL_D_COVARIANT_LOCAL_MATCHING_AND_COMPLETE_GAUSSIAN_METRIC_FOUR_POINT_WITH_EXACT_LIGHT_AND_BOUNDED_FIXED_HEAVY_FORWARD_COEFFICIENTS",
        "date": "2026-09-15",
        "status": "SCOPED_COVARIANT_GAUSSIAN_FOUR_POINT_MATCHING; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/local.md",
            "notes/gaussian.md",
            "notes/forward.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_covariant_local_map": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_Gaussian_four_point_and_forward_coefficients": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full dimension-dependent covariant local map has rank two with six exact null directions, validated by literal Fourier vertices and complete conserved metric response. The entire quadratic Euler kernel vanishes before the D limit, while evanescent local projection terms remain explicit. Complete H/Proca Gaussian metric four-point functions match their fixed prescriptions and all-angle physical cuts. Their full known-pole-subtracted b20 is strictly negative and below10^-1600 in magnitude; the exact light nonlocal constant is supplied without choosing its local curvature terms. No full-source positivity or finite IR/Regge conclusion follows. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Finite light or other non-Gaussian local matching coefficients or full physical Newton sign",
            "Complete full-source b20, all C/g/heavy/non-Gaussian sectors or higher-loop errors",
            "Finite physical detector/dressing/IR or fixed-transfer Regge contour remainder",
            "An exact all-order quantum EOM or field-redefinition equivalence, or a changed physical matter frame",
            "Original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact full-D local projection and nullspace, literal four-scalar Fourier vertices, independent full metric-response and quadratic Euler components in D4/D5/D6, frozen radial kernel and whole physical cut comparisons, exact rational integration and all-domain analytic bounds support the written first-insertion and gapped matching proofs. These are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; exact-GCD is FULL-only. No frozen source/raw, finite light coefficient or physical frame is changed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The covariant Gaussian four-point report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.289 covariant Gaussian four-point matching replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
