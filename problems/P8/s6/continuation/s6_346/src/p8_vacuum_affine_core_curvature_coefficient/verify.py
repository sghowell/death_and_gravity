"""Read-only known polynomial core curvature matching in one fixed jet convention."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_box_curvature_coefficient import verify as input0
from p8_vacuum_affine_curvature_contact_basis import verify as input3
from p8_vacuum_affine_factorized_bubble_radiation import verify as input2
from p8_vacuum_affine_heavy_parent_one_loop import verify as input4
from p8_vacuum_affine_heavy_scalar_four_point_loop import verify as input5
from p8_vacuum_affine_triangle_curvature_coefficient import verify as input1

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-core-curvature-coefficient.json"
PARENT_SHA = "bcc2ef4ca3ccd84c523a0159e57a34289e99a35976a7f99201fd6384735c52bd"
TRIANGLE_SHA = "227850342699960004105744b9eaad57b20a626b19d2e90deb353512c847a036"
BUBBLE_SHA = "174317dc92bc5b1c2956e6c5bebee4504f582325ae7e28b16f03d5de01f70206"
BASIS_SHA = "dd49e3ecf06b03e905779e22625ad57147f3ab1d2b69f5d49896cea55338a89e"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
FLAT_SHA = "1277412d19d9b022eb81c3728251d1836d88e4b047c8ad80794107c2aff7c911"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_core_curvature_coefficient/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (input0, PARENT_SHA),
        (input1, TRIANGLE_SHA),
        (input2, BUBBLE_SHA),
        (input3, BASIS_SHA),
        (input4, SOURCE_SHA),
        (input5, FLAT_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen common-core matching input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_345_complete_box_in_fixed_six_word_lift_rebuilt": PARENT_SHA,
        "S6_344_full_labeled_triangle_coefficient_and_dressing_rebuilt": TRIANGLE_SHA,
        "S6_339_complete_original_bubble_and_three_fixed_counterterms_rebuilt": BUBBLE_SHA,
        "S6_343_complete_first_real_TT_curvature_basis_rebuilt": BASIS_SHA,
        "S6_239_original_source_and_pole_only_derivative_scheme_rebuilt": SOURCE_SHA,
        "S6_235_complete_original_polynomial_loop_normalization_rebuilt": FLAT_SHA,
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
        raise ValueError("A known core coefficient proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.346.KNOWN_POLYNOMIAL_CORE_CURVATURE_IN_ONE_COVARIANT_JET_BASIS",
        "date": "2026-09-18",
        "status": "SCOPED_KNOWN_POLYNOMIAL_CORE_DEGREE6_CURVATURE_MATCHING; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/conversion.md",
            "notes/moments.md",
            "notes/core.md",
            "notes/bounds.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_common_basis_conversion_and_exact_moments": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_core_sum_bounds_and_original_calibration": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete known old-polynomial bubble, triangle and alternating-mass box degree6 real-TT curvature coefficient is computed in one explicit six-word covariant jet lift. Full off-shell conversion identities and exact mass moments precede the sum. For the unchanged tuned contact, the n^-2 and n^-3 terms cancel exactly and n^4*c_core tends to -23/105. An exact finite rational-log bound, not this limit alone, proves -1/(4n^4)<c_core<-1/(5n^4) and the original |chi_core|/A0<10^-207. This does not include the new local/mixed-source common-basis matching or independent extra parent chi; no physical above-threshold Taylor replacement, second copy of known loops, or original V/G/B/P8 closure is asserted.",
        "not_established": [
            "The entire known matter coefficient before off-shell conversion of the new local tadpole and mixed-source classes",
            "The independent added parent-theory curvature coefficient chi or complete physical matching",
            "A physical above-threshold amplitude approximation or local Taylor truncation-error bound",
            "Internal-graviton loops, a full curved-background counterfunctional or finite-gravity quantum decoupling",
            "Full virtual/inclusive matching, exact LSZ, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "All10 labeled-word flat projections and complete generic covariant response differences are recomputed in exact symbolic arithmetic. All four heavy-dressing orders, exact triangle primitives/equal-mass limits, original box moment, bubble integral normalization and tuned sum are retained. Every term of the original rational-log numerator enters a finite rational tail bound; the log bound uses a positive exact exponential partial sum. Sixty original-vector word checks and six actual frozen bubble kernels independently test normalization. All direct raw reports and source manifests remain immutable. Not kernel-formalized and not full nonlocal parent matching.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete common-core curvature report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.346 known polynomial core common-basis curvature replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
